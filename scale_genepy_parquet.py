#!/usr/bin/env python3
"""Memory-bounded, resumable, parallel GenePy normalization.

Input: a wide Parquet table with sample_ids and one float64 column per gene.
Output: minmax/block_XXXXX.parquet and quantile/block_XXXXX.parquet, each
containing sample_ids followed by the corresponding original gene columns.

Transformations are the same scikit-learn estimators used in scale_matrix.py:
  MinMaxScaler()
  QuantileTransformer(n_quantiles=1000, random_state=0)

Importantly, this script validates that column-block quantile fitting agrees
with a whole-matrix fit on the installed scikit-learn version before running.
It can also test real rows of the input using --validate-only.
"""

import argparse
import concurrent.futures
import gc
import hashlib
import json
import os
from pathlib import Path
import sys
import time

# Avoid BLAS/OpenMP oversubscription in each worker. Set before NumPy import.
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
import pandas as pd
import sklearn
from sklearn.preprocessing import MinMaxScaler, QuantileTransformer


def qt():
    # Do not change subsample, output_distribution, or default tie handling.
    return QuantileTransformer(n_quantiles=1000, random_state=0)


def compare(actual, expected, label, tolerance=1e-12):
    if actual.shape != expected.shape:
        raise AssertionError(f"{label}: shape mismatch {actual.shape} vs {expected.shape}")
    if not np.array_equal(np.isnan(actual), np.isnan(expected)):
        raise AssertionError(f"{label}: NaN positions differ")
    valid = ~(np.isnan(actual) | np.isnan(expected))
    diff = np.abs(actual[valid] - expected[valid])
    max_diff = float(diff.max()) if diff.size else 0.0
    nbad = int(np.count_nonzero(diff > tolerance))
    print(f"{label}: max_diff={max_diff:.3g}; values > {tolerance:g}: {nbad}", flush=True)
    if nbad:
        raise AssertionError(f"{label}: numerical results differ")


def verify_block_equivalence():
    # Both below/above default sklearn subsampling cutoff, with many ties/NaNs.
    rng = np.random.RandomState(337)
    nrows, ncols = 12007, 12
    matrix = rng.randint(0, 20, size=(nrows, ncols)).astype(np.float64)
    matrix[:, 0] = 0  # constant all-zero gene
    matrix[:, 1] = rng.normal(size=nrows)
    matrix[::89, 2] = np.nan
    matrix[::43, 3] = 0
    matrix[::47, 7] = np.nan
    full = qt().fit_transform(matrix)
    chunked = np.column_stack([
        qt().fit_transform(matrix[:, begin:begin + 3])
        for begin in range(0, ncols, 3)
    ])
    compare(chunked, full, "Synthetic full-vs-column-block QuantileTransformer")


def get_schema(path, id_column):
    import pyarrow.parquet as pq
    parquet = pq.ParquetFile(str(path))
    fields = parquet.schema_arrow
    names = fields.names
    if names.count(id_column) != 1:
        raise ValueError(f"Expected one identifier column called {id_column!r}, got {names.count(id_column)}")
    genes = [name for name in names if name != id_column]
    if len(genes) != len(set(genes)):
        raise ValueError("Duplicate gene column names")
    if not genes:
        raise ValueError("No gene columns found")
    return parquet.metadata.num_rows, genes


def read_columns(path, id_column, gene_names):
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(str(path))
    tbl = pf.read(columns=[id_column] + list(gene_names), use_threads=False)
    df = tbl.to_pandas(ignore_metadata=True)  # Keep sample_ids as an ordinary column
    del tbl, pf
    if df.columns.tolist() != [id_column] + list(gene_names):
        df = df.loc[:, [id_column] + list(gene_names)]
    if df.shape[1] != len(gene_names) + 1:
        raise ValueError("Input columns could not be read as expected")
    data = df.loc[:, list(gene_names)].to_numpy(dtype=np.float64, copy=False)
    return df[id_column], data, df


def validate_real(path, id_column, row_limit, gene_limit, block_size):
    _, genes = get_schema(path, id_column)
    selected = genes[:min(gene_limit, len(genes))]
    print(f"Validating {row_limit} rows and {len(selected)} genes from {path}", flush=True)
    _, all_data, backing_df = read_columns(path, id_column, selected)
    data = np.array(all_data[:row_limit, :], dtype=np.float64, copy=True)
    del backing_df, all_data
    if data.shape[0] <= 10000:
        print("WARNING: validation has <=10,000 rows; it will not exercise default subsampling.", flush=True)
    t0 = time.monotonic()
    expected_mm = MinMaxScaler().fit_transform(data)
    expected_qt = qt().fit_transform(data)
    actual_mm = np.column_stack([
        MinMaxScaler().fit_transform(data[:, j:j + block_size])
        for j in range(0, data.shape[1], block_size)
    ])
    actual_qt = np.column_stack([
        qt().fit_transform(data[:, j:j + block_size])
        for j in range(0, data.shape[1], block_size)
    ])
    compare(actual_mm, expected_mm, "Real-data MinMax")
    compare(actual_qt, expected_qt, "Real-data Quantile")
    print(f"Real-data validation PASSED in {time.monotonic() - t0:.1f}s", flush=True)


def atomic_json(path, obj):
    temp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(temp, "w", encoding="utf-8") as out:
        json.dump(obj, out, indent=2, sort_keys=True)
    os.replace(temp, path)


def one_block(task):
    input_path, output_dir, id_column, block_num, cols, compression, level, nrows = task
    base = Path(output_dir)
    name = f"block_{block_num:05d}.parquet"
    mm_path = base / "minmax" / name
    qt_path = base / "quantile" / name
    done_path = base / "completed" / f"block_{block_num:05d}.json"
    if done_path.exists() and mm_path.exists() and qt_path.exists():
        with open(done_path, encoding="utf-8") as inp:
            meta = json.load(inp)
        if meta.get("genes") == list(cols) and meta.get("rows") == nrows:
            return block_num, len(cols), "skipped", 0.0

    tmp_mm = mm_path.with_name(mm_path.name + f".tmp.{os.getpid()}")
    tmp_qt = qt_path.with_name(qt_path.name + f".tmp.{os.getpid()}")
    start = time.monotonic()
    try:
        sample, x, backing = read_columns(input_path, id_column, cols)
        if len(sample) != nrows:
            raise RuntimeError(f"Block {block_num}: unexpected row count")

        # Keep a single input block in memory; write and release MinMax before Quantile.
        mm = MinMaxScaler().fit_transform(x)
        out = pd.DataFrame(mm, columns=cols, copy=False)
        out.insert(0, id_column, sample.to_numpy(copy=False))
        out.to_parquet(tmp_mm, engine="pyarrow", compression=compression,
                       compression_level=level, index=False)
        del mm, out
        gc.collect()

        transformed = qt().fit_transform(x)
        out = pd.DataFrame(transformed, columns=cols, copy=False)
        out.insert(0, id_column, sample.to_numpy(copy=False))
        out.to_parquet(tmp_qt, engine="pyarrow", compression=compression,
                       compression_level=level, index=False)
        del transformed, out, x, backing
        gc.collect()

        os.replace(tmp_mm, mm_path)
        os.replace(tmp_qt, qt_path)
        atomic_json(done_path, {"genes": list(cols), "rows": nrows, "block": block_num})
        return block_num, len(cols), "done", time.monotonic() - start
    finally:
        for p in (tmp_mm, tmp_qt):
            if p.exists():
                p.unlink()


def run(args):
    source = Path(args.input).resolve()
    output = Path(args.output).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if args.jobs < 1 or args.block_size < 1:
        raise ValueError("--jobs and --block-size must be positive")

    print(f"scikit-learn {sklearn.__version__}, pandas {pd.__version__}, NumPy {np.__version__}", flush=True)
    verify_block_equivalence()

    nrows, genes = get_schema(source, args.id_column)
    print(f"Input: {nrows:,} rows x {len(genes):,} genes; "
          f"raw float64 scores {nrows * len(genes) * 8 / 1024**3:.2f} GiB", flush=True)
    if args.validate_only:
        validate_real(source, args.id_column, args.validation_rows,
                      args.validation_genes, min(args.block_size, 32))
        return

    # Guard against accidentally resuming with a different input/configuration.
    st = source.stat()
    signature = {
        "input_path": str(source), "input_bytes": st.st_size,
        "input_mtime_ns": st.st_mtime_ns,
        "input_rows": nrows, "n_genes": len(genes),
        "columns_sha256": hashlib.sha256("\n".join(genes).encode()).hexdigest(),
        "id_column": args.id_column, "block_size": args.block_size,
        "sklearn_version": sklearn.__version__,
        "quantile_n_quantiles": 1000, "quantile_random_state": 0,
        "quantile_default_subsample": qt().subsample,
        "compression": args.compression, "compression_level": args.compression_level,
    }
    output.mkdir(parents=True, exist_ok=True)
    for name in ("minmax", "quantile", "completed"):
        (output / name).mkdir(exist_ok=True)
    manifest = output / "run_manifest.json"
    if manifest.exists():
        with open(manifest, encoding="utf-8") as fh:
            old = json.load(fh)
        if old != signature:
            raise RuntimeError("Existing output was generated from a different input or configuration; use a new output directory.")
    else:
        if any((output / "completed").glob("*.json")):
            raise RuntimeError("Existing completion markers with no run manifest")
        atomic_json(manifest, signature)

    blocks = [genes[start:start + args.block_size]
              for start in range(0, len(genes), args.block_size)]
    print(f"Scheduling {len(blocks)} column blocks, {args.jobs} processes, "
          f"{args.block_size} genes/block, compression={args.compression}", flush=True)
    jobs = [
        (str(source), str(output), args.id_column, i, chunk,
         args.compression, args.compression_level, nrows)
        for i, chunk in enumerate(blocks)
    ]
    start = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        future_jobs = {pool.submit(one_block, task): task[3] for task in jobs}
        completed = 0
        for future in concurrent.futures.as_completed(future_jobs):
            block_num, n_genes, state, duration = future.result()
            completed += 1
            print(f"[{completed:3d}/{len(blocks)}] block {block_num:05d}: "
                  f"{state:7s}, genes={n_genes:3d}, duration={duration:.1f}s", flush=True)

    elapsed = time.monotonic() - start
    atomic_json(output / "all_blocks_complete.json", {
        "blocks": len(blocks), "rows": nrows, "genes": len(genes),
        "source": str(source), "seconds": elapsed,
        "note": "One Parquet file per column block; join on sample_ids if a wide table is required."
    })
    print(f"COMPLETE: {len(genes):,} genes in {elapsed:.1f} seconds", flush=True)
    print(f"Outputs: {output / 'minmax'} and {output / 'quantile'}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", help="Wide Parquet matrix; sample IDs and numeric gene columns")
    parser.add_argument("output", nargs="?", default="scaled_genepy", help="Output directory")
    parser.add_argument("--id-column", default="sample_ids")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--block-size", type=int, default=128)
    parser.add_argument("--compression", default="zstd", choices=("zstd", "snappy", "none"))
    parser.add_argument("--compression-level", type=int, default=3)
    parser.add_argument("--validate-only", action="store_true", help="Validate exact mathematics on subset; write no output")
    parser.add_argument("--validation-rows", type=int, default=25001)
    parser.add_argument("--validation-genes", type=int, default=128)
    args = parser.parse_args()
    if args.compression == "none":
        args.compression = None
        args.compression_level = None
    run(args)


if __name__ == "__main__":
    main()
