#!/usr/bin/env python3
"""Merge column-partitioned GenePy Parquet files with bounded memory.

One output Parquet per normalization (minmax_full.parquet, quantile_full.parquet).
Streams row stripes from every block, preserving the input column order and
verifying every sample ID. Never holds the 200,643 x 24,625 matrix in RAM.

Usage:
  python merge_genepy_blocks_streaming.py INPUT.parquet BLOCK_OUTPUT_DIR \
      --rows-per-batch 8192
"""
import argparse
from contextlib import ExitStack
import gc
import os
from pathlib import Path
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


class BatchCursor:
    """Return an exact number of rows, even across Parquet row-group boundaries."""
    def __init__(self, parquet_file, batch_size):
        self._batches = iter(parquet_file.iter_batches(batch_size=batch_size, use_threads=False))
        self._batch = None
        self._offset = 0

    def take(self, n):
        remaining = n
        pieces = []
        while remaining:
            if self._batch is None or self._offset >= self._batch.num_rows:
                try:
                    self._batch = next(self._batches)
                except StopIteration as exc:
                    raise ValueError('A source Parquet file ended early') from exc
                self._offset = 0
            k = min(remaining, self._batch.num_rows - self._offset)
            pieces.append(self._batch.slice(self._offset, k))
            self._offset += k
            remaining -= k
        return pa.Table.from_batches(pieces)


def check_existing(dest, names, nrows):
    pf = pq.ParquetFile(dest)
    try:
        if pf.metadata.num_rows != nrows or pf.schema_arrow.names != names:
            raise ValueError(f'Existing output has incorrect rows/columns: {dest}')
    finally:
        pf.close()
    print(f'SKIPPING complete existing file: {dest}', flush=True)


def merge_one(kind, input_file, outdir, names, source_ids, block_size, rows_per_batch):
    genes = names[1:]
    nrows = len(source_ids)
    nblocks = (len(genes) + block_size - 1) // block_size
    dest = outdir / f'{kind}_full.parquet'
    tmp = outdir / f'{kind}_full.parquet.tmp.{os.getpid()}'
    t0 = time.monotonic()
    print(f'\nMerging {kind}: {nblocks} blocks, {nrows:,} rows', flush=True)

    pfs = []
    try:
        # Keep each Parquet file open and advance its record-batch iterator
        # in lockstep; this avoids repeatedly decompressing whole block files.
        for i in range(nblocks):
            path = outdir / kind / f'block_{i:05d}.parquet'
            if not path.is_file():
                raise FileNotFoundError(path)
            pf = pq.ParquetFile(path)
            expected = [names[0]] + genes[i * block_size:(i + 1) * block_size]
            if pf.schema_arrow.names != expected:
                raise ValueError(f'Unexpected column names/order in {path}')
            if pf.metadata.num_rows != nrows:
                raise ValueError(f'Unexpected row count in {path}')
            if any(f.type != pa.float64() for f in list(pf.schema_arrow)[1:]):
                raise ValueError(f'Non-float64 gene column in {path}')
            pfs.append(pf)

        cursors = [BatchCursor(pf, rows_per_batch) for pf in pfs]
        writer = None
        try:
            for start in range(0, nrows, rows_per_batch):
                n = min(rows_per_batch, nrows - start)
                expected_ids = source_ids[start:start+n]
                arrays = []
                for i, cursor in enumerate(cursors):
                    block_table = cursor.take(n)
                    actual_ids = block_table.column(0).to_numpy(zero_copy_only=False)
                    if not np.array_equal(actual_ids, expected_ids):
                        raise ValueError(f'{kind}, block {i}, row {start}: sample_ids mismatch')
                    if i == 0:
                        arrays.append(block_table.column(0))
                    arrays.extend(block_table.columns[1:])
                stripe = pa.Table.from_arrays(arrays, names=names)
                if writer is None:
                    writer = pq.ParquetWriter(
                        str(tmp), stripe.schema,
                        compression='zstd', compression_level=3,
                        use_dictionary=False, write_statistics=False,
                    )
                writer.write_table(stripe, row_group_size=n)
                print(f'  {start+n:,}/{nrows:,} rows written', flush=True)
                del stripe, arrays, actual_ids, block_table
                gc.collect()
        finally:
            if writer is not None:
                writer.close()

        result = pq.ParquetFile(tmp)
        try:
            if result.metadata.num_rows != nrows or result.schema_arrow.names != names:
                raise AssertionError('Merged output dimensions/order do not match input')
        finally:
            result.close()
        os.replace(tmp, dest)
        print(f'COMPLETE {kind}: {dest} ({dest.stat().st_size / 1024**3:.3f} GiB), '
              f'{time.monotonic()-t0:.1f}s', flush=True)
    finally:
        for pf in pfs:
            pf.close()
        if tmp.exists():
            tmp.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Original UKB Parquet')
    parser.add_argument('output_dir', type=Path, help='Directory with minmax/ and quantile/ blocks')
    parser.add_argument('--id-column', default='sample_ids')
    parser.add_argument('--block-size', type=int, default=128)
    parser.add_argument('--rows-per-batch', type=int, default=8192,
                        help='Rows written per stripe; smaller uses less RAM')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if min(args.block_size, args.rows_per_batch, args.threads) < 1:
        parser.error('Positive block-size, rows-per-batch and threads required')
    pa.set_cpu_count(args.threads)

    inp = pq.ParquetFile(args.input)
    try:
        names = inp.schema_arrow.names
        nrows = inp.metadata.num_rows
        if not names or names[0] != args.id_column:
            raise ValueError(f'First input field must be {args.id_column!r}')
        if len(names) != len(set(names)):
            raise ValueError('Duplicate input column names')
    finally:
        inp.close()
    # Only sample IDs are read from the original wide file (~1.5 MiB).
    source_ids = pq.read_table(args.input, columns=[args.id_column]).column(0).to_numpy()
    print(f'Input: {nrows:,} rows, {len(names)-1:,} genes', flush=True)
    if len(np.unique(source_ids)) != nrows:
        raise ValueError('Input sample_ids are not unique')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for kind in ('minmax', 'quantile'):
        dest = args.output_dir / f'{kind}_full.parquet'
        if dest.exists() and not args.overwrite:
            check_existing(dest, names, nrows)
            continue
        merge_one(kind, args.input, args.output_dir, names, source_ids,
                  args.block_size, args.rows_per_batch)
    print('\nALL MERGES COMPLETE', flush=True)

if __name__ == '__main__':
    main()
