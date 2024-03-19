# -*- coding: utf-8 -*-

import sys, re
import numpy as np
import pandas as pd
import time
from numba import njit, cuda
import math
import pyarrow.csv as pa_csv

threads_per_block = 256

@njit(nopython=True)
def zeroing_with_numba(db1,to_be_zeroed):

    for i in range(db1.shape[0]):
        for j in range(db1.shape[1]):
            if db1[i, j] in to_be_zeroed:
                db1[i, j] = 0.0


@cuda.jit
def s_af0_af_func_kernel(S, af, af0, db1, i, mask_indices):
    start = cuda.grid(1) 
    stride = cuda.gridsize(1)
     
    for r in range(start, mask_indices.shape[0], stride):
        id_r = mask_indices[r, 0]
        id_c = mask_indices[r,1]
        db1[id_r,id_c] = S[id_r, i] * (-math.log10(af0[id_r] * af[id_r, i]))
    
    mask_indices = None

def s_af0_ad_func_dum(S, af, af0, db1, i, mask_indices):
    for r in range(mask_indices.shape[0]):
        id_r = mask_indices[r, 0]
        id_c = mask_indices[r,1]
        db1[id_r,id_c] = S[id_r, i] * (-math.log10(af0[id_r] * af[id_r, i]))

def s_af0_af_func(mask, S, af, af0, db1, i):
    #print("indices",np.where(mask))
    mask_indices = np.column_stack(np.where(mask))
    if mask_indices.size == 0:
        return
    #print("stacked indices",mask_indices)
    num_rows_mask = mask_indices.shape[0]
    mask_indices = cuda.to_device(mask_indices)
    func_threads = min(threads_per_block, num_rows_mask)
    blockspergrid = (num_rows_mask + func_threads - 1) // func_threads
    print(num_rows_mask,blockspergrid, func_threads)
    s_af0_af_func_kernel[blockspergrid, func_threads](S, af, af0, db1, i, mask_indices)
    #s_af0_af_func_kernel[1,1](S, af, af0, db1, i, mask_indices)#s_af0_ad_func_dum(S, af, af0, db1, i, mask_indices)
    mask_indices = None

def S_af0_af_num(i, S, af, af0, db1, samples):
    masks = [samples == '0/'+str(i+1), samples == str(i+1)+'/0']
    for mask in masks:
        s_af0_af_func(mask, S, af, af0, db1, i)


def S_af0_af(i: int, S: np.array, af: np.array, af0: np.array, db1: np.array, samples: np.array):
    def _func(mask):
        db1[mask]=S[np.where(mask)[0],i]*(-np.log10((af0[np.where(mask)[0]])*(af[np.where(mask)[0],i])))
    masks= [samples == '0/'+str(i+1),samples == str(i+1)+'/0']
    for mask in masks:
        _func(mask)

def S_af_af(i: int, S: np.array, af: np.array, db1: np.array, samples:np.array):
    mask = samples == str(i+1)+'/'+str(i+1)
    db1[mask]=S[np.where(mask)[0],i]*(-np.log10((af[np.where(mask)[0],i])*(af[np.where(mask)[0],i])))


@cuda.jit
def s_af_af_func_kernel(S, af, db1, i, mask_indices):
    start = cuda.grid(1) 
    stride = cuda.gridsize(1)
     
    for r in range(start, mask_indices.shape[0], stride):
        id_r = mask_indices[r, 0]
        id_c = mask_indices[r,1]
        db1[id_r,id_c]=S[id_r,i]*(-math.log10((af[id_r,i])*(af[id_r,i])))
    
    mask_indices = None


def s_af_af_func(mask, S, af, db1, i):
    mask_indices = np.column_stack(np.where(mask))
    if mask_indices.size == 0:
        return
    #print("stacked indices",mask_indices)
    num_rows_mask = mask_indices.shape[0]
    mask_indices = cuda.to_device(mask_indices)
    func_threads = min(threads_per_block, num_rows_mask)
    blockspergrid = (num_rows_mask + func_threads - 1) // func_threads
    print(num_rows_mask,blockspergrid, func_threads)

    s_af_af_func_kernel[blockspergrid, func_threads](S, af, db1, i, mask_indices)
    mask_indices = None

def S_af_af_num(i, S, af, db1, samples):
    mask = samples == str(i+1)+'/'+str(i+1)
    s_af_af_func(mask, S, af, db1, i)


def S_S_af_af_2(i: int, S: np.array, af: np.array, db1: np.array,samples: np.array):

    for j in range(i+1,9):
        masks = [samples == str(i+1)+'/'+str(j+1),samples== str(i+1)+'/'+str(j+1)]

        for mask in masks:
            db1[mask]=((S[np.where(mask)[0],i]+S[np.where(mask)[0],j])/2)*(-np.log10((af[np.where(mask)[0],i])*(af[np.where(mask)[0],j])))

@cuda.jit
def s_s_af_af_2_func_kernel(S, af, db1, i,j, mask_indices):
    start = cuda.grid(1) 
    stride = cuda.gridsize(1)
    #rows_per_block = (mask_indices.shape[0] + threads_per_block - 1) // threads_per_block
     
    for r in range(start, mask_indices.shape[0], stride):
        id_r = mask_indices[r, 0]
        id_c = mask_indices[r,1]
        db1[id_r,id_c]=((S[id_r,i]+S[id_r,j])/2)*(-math.log10((af[id_r,i])*(af[id_r,j])))
    
    mask_indices = None


def s_s_af_af_2_func(mask, S, af, db1, i,j):
    mask_indices = np.column_stack(np.where(mask))
    if mask_indices.size == 0:
        return
    #print("stacked indices",mask_indices)
    num_rows_mask = mask_indices.shape[0]
    mask_indices = cuda.to_device(mask_indices)
    func_threads = min(threads_per_block, num_rows_mask)
    blockspergrid = (num_rows_mask + func_threads - 1) // func_threads
    print(num_rows_mask,blockspergrid, func_threads)

    s_s_af_af_2_func_kernel[blockspergrid, func_threads](S, af, db1, i, j, mask_indices)
    mask_indices = None

def S_S_af_af_2_num(i, S, af, db1, samples):
    for j in range(i+1,9):
        masks = [samples == str(i+1)+'/'+str(j+1),samples== str(i+1)+'/'+str(j+1)]
        for mask in masks:
            s_s_af_af_2_func(mask,S,af,db1,i,j)

def S_S_af_af_1(i: int, S: np.array, af: np.array, db1: np.array,samples: np.array):
    for j in range(i+1,9):
        mask = samples == str(i+1)+'/'+str(j+1)
        db1[mask]=((S[np.where(mask)[0],i]+S[np.where(mask)[0],j])/2) * (-math.log10((af[np.where(mask)[0],i])*(af[np.where(mask)[0],j])))


@cuda.jit
def s_s_af_af_1_func_kernel(S, af, db1, i,j, mask_indices):
    start = cuda.grid(1) 
    stride = cuda.gridsize(1)
    #rows_per_block = (mask_indices.shape[0] + threads_per_block - 1) // threads_per_block
     
    for r in range(start, mask_indices.shape[0], stride):
        id_r = mask_indices[r, 0]
        id_c = mask_indices[r,1]
        db1[id_r,id_c]=((S[id_r,i]+S[id_r,j])/2)*(-math.log10((af[id_r,i])*(af[id_r,j])))
    
    mask_indices = None


def s_s_af_af_1_func(mask, S, af, db1, i,j):
    mask_indices = np.column_stack(np.where(mask))
    if mask_indices.size == 0:
        return
    #print("stacked indices",mask_indices)
    num_rows_mask = mask_indices.shape[0]
    mask_indices = cuda.to_device(mask_indices)
    func_threads = min(threads_per_block, num_rows_mask)
    blockspergrid = (num_rows_mask + func_threads - 1) // func_threads
    print(num_rows_mask,blockspergrid, func_threads)

    s_s_af_af_1_func_kernel[blockspergrid, func_threads](S, af, db1, i, j, mask_indices)
    mask_indices = None

def S_S_af_af_1_num(i, S, af, db1, samples):
    for j in range(i+1,9):
        mask = samples == str(i+1)+'/'+str(j+1)
        s_s_af_af_1_func(mask,S,af,db1,i,j)

# IMPORT data passed through the 1st argument;followed by the Gene_ENSGid as the 2nd
#data=pd.read_csv(sys.argv[1],sep='\t',header=0,dtype='object',engine="pyarrow")
t0=time.time()
parse_options = pa_csv.ParseOptions(delimiter='\t')
read_options=pa_csv.ReadOptions(block_size=1e9)
data = pa_csv.read_csv(sys.argv[1], parse_options=parse_options, read_options = read_options)
data = data.to_pandas()
t1=time.time()
#data = data.to_numpy()
gene=sys.argv[2]
CADD=sys.argv[3]
header=data.columns.values
print(f"Time to read: {gene} {t1-t0}")

data = np.array(data)
##cadd score; reformat CADD range to 0-1
## the UKBB 200k cohort raw variant has CADD16 ranges from -18.793437 to 19.100986 including the prescore and permutated score (GRCh38_v1.6;VEP 100_GRCh38)
scores = data[:,16:25]
scores = scores.astype('float')
scores = (scores - (-19.811548))/(25.028523-(-19.811548))
scores[np.isnan(scores)] = 0


##allele frequency as it is in the UKBB cohort; this is currently based on the raw pVCF data.
af = data[:,6:15]
af = af.astype('float')
##the ref allele frequency
af0 = 1-np.nansum(af,axis=1)
af[np.isnan(af)] = 1
samples=data[:,26:]
samples_header=header[26:]
print("FInished wrangling")
threads_per_block = 256
len_samples = len(samples)
#print("sample length",len(samples))
blockspergrid = (len_samples + (threads_per_block - 1)) // threads_per_block
#print("blocks",blockspergrid, threads_per_block, len_samples)

##Clalculate GenePy score; assuming genotype 0/0 has a score=0
def nan_if(arr, value):
    return np.where(arr == value, np.nan, arr)

def score_db(samples,scores,af0,af):
    
    db1=np.zeros_like(samples, dtype=float)
    S_device =cuda.to_device(scores)
    db1_device = cuda.to_device(db1)
    af0_device = cuda.to_device(af0)
    af_device = cuda.to_device(af)

    out1=[]
    to_be_zeroed = [ '0/0', '0',  './0', './.', './1', './2', './3', 
                    './4', './5', './6', './7', './8', './9', '10/', 
                    '11/', '12/', '13/', '14/', '15/', '16/', '17/', 
                    '18/', '19/', '20/', '21/', '22/', '23/', '24/', 
                    '25/', '26/', '27/','28/', '30/']

    t0=time.time()
    zeroing_with_numba(db1,to_be_zeroed)
    t1=time.time()
    print("Time to zero {t1-t0}")
    t0_l=time.time()
    for index in range(9):
        t0 = time.time()
        S_af0_af_num(index,S_device,af_device,af0_device,db1_device,samples)
        t1 = time.time()

        print(f"S_af0_af Loop {index} took {t1-t0}")
        t0=time.time()
        S_af_af_num(index,S_device,af_device,db1_device,samples)
        t1 = time.time()
        print(f"S_af_af Loop {index} took {t1-t0}")
        t0=time.time()
        S_S_af_af_2_num(index,S_device,af_device,db1_device,samples)
        t1=time.time()
        print(f"S_S_af_af_2 Loop {index} took {t1-t0}")
    for i in range(3,9):
        t0=time.time()
        S_S_af_af_1_num(i,S_device,af_device,db1_device,samples)
        t1=time.time()
        print(f"S_S_af_af_1 Loop {i} took {t1-t0}")
    t1_l=time.time()
    print(f"total loop time {t1_l-t0_l}")
    #    S_af_af_num(index,S,af,db1,samples)





    #    S_S_af_af_2_num(index,S,af,db1,samples)
    #for i in range(3,9):
    #    S_S_af_af_1_num(i,S,af,db1,samples)
    db1 = db1_device.copy_to_host()
    db1_device = None
    out1=np.nansum(nan_if(db1,'0.0'),axis=0)
    gg = np.array([gene]*len(samples_header))
    U = np.vstack((samples_header,out1,gg)).T
    #print("U",U.shape)
    return U
#%%

if (np.isnan(scores).sum()) < (scores.shape[0]): #compute metascores if at least 1 variant
    #print("In the if")
    t0=time.time()
    U = score_db(samples,scores,af0,af)
    t1=time.time()
    print(f"Time for score_db: {gene}: {t1-t0}")
    np.savetxt('./'+gene+'_'+CADD+'_matrix',U, fmt='%s', delimiter='\t')




