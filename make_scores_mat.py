# -*- coding: utf-8 -*-

import sys, re
import numpy as np
import pandas as pd


# IMPORT data passed through the 1st argument;followed by the Gene_ENSGid as the 2nd

data=pd.read_csv(sys.argv[1],sep='\t',header=0,dtype='object')
gene=sys.argv[2]
CADD=sys.argv[3]
header=data.columns.values
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

##sample genotype data as it is in format of allele_1/allele_2
samples=data[:,26:]
samples_header=header[26:]

##Clalculate GenePy score; assuming genotype 0/0 has a score=0
def nan_if(arr, value):
    return np.where(arr == value, np.nan, arr)

def score_db(samples,scores,af0,af):
    S=np.copy(scores)
    db1=np.copy(samples)
    af0=np.copy(af0)
    af=np.copy(af)
    out1=[]
    for i in range(db1.shape[0]):
        db1[i][db1[i]=='0/0']=0
        db1[i][db1[i]=='0']=0
        db1[i][db1[i]=='./0']=0
        db1[i][db1[i]=='./.']=0
        db1[i][db1[i]=='./1']=0
        db1[i][db1[i]=='./2']=0
        db1[i][db1[i]=='./3']=0
        db1[i][db1[i]=='./4']=0
        db1[i][db1[i]=='./5']=0
        db1[i][db1[i]=='./6']=0
        db1[i][db1[i]=='./7']=0
        db1[i][db1[i]=='./8']=0
        db1[i][db1[i]=='./9']=0
        db1[i][db1[i]=='10/']=0
        db1[i][db1[i]=='11/']=0
        db1[i][db1[i]=='12/']=0
        db1[i][db1[i]=='13/']=0
        db1[i][db1[i]=='14/']=0
        db1[i][db1[i]=='15/']=0
        db1[i][db1[i]=='16/']=0
        db1[i][db1[i]=='17/']=0
        db1[i][db1[i]=='18/']=0
        db1[i][db1[i]=='19/']=0
        db1[i][db1[i]=='20/']=0
        db1[i][db1[i]=='21/']=0
        db1[i][db1[i]=='22/']=0
        db1[i][db1[i]=='23/']=0
        db1[i][db1[i]=='24/']=0
        db1[i][db1[i]=='25/']=0
        db1[i][db1[i]=='26/']=0
        db1[i][db1[i]=='27/']=0
        db1[i][db1[i]=='28/']=0
        db1[i][db1[i]=='30/']=0
        db1[i][db1[i]=='0/1']=S[i,0]*(-np.log10((af0[i])*(af[i,0])))
        db1[i][db1[i]=='1/0']=S[i,0]*(-np.log10((af0[i])*(af[i,0])))
        db1[i][db1[i]=='0/2']=S[i,1]*(-np.log10((af0[i])*(af[i,1])))
        db1[i][db1[i]=='2/0']=S[i,1]*(-np.log10((af0[i])*(af[i,1])))
        db1[i][db1[i]=='0/3']=S[i,2]*(-np.log10((af0[i])*(af[i,2])))
        db1[i][db1[i]=='3/0']=S[i,2]*(-np.log10((af0[i])*(af[i,2])))
        db1[i][db1[i]=='0/4']=S[i,3]*(-np.log10((af0[i])*(af[i,3])))
        db1[i][db1[i]=='4/0']=S[i,3]*(-np.log10((af0[i])*(af[i,3])))
        db1[i][db1[i]=='0/5']=S[i,4]*(-np.log10((af0[i])*(af[i,4])))
        db1[i][db1[i]=='5/0']=S[i,4]*(-np.log10((af0[i])*(af[i,4])))
        db1[i][db1[i]=='0/6']=S[i,5]*(-np.log10((af0[i])*(af[i,5])))
        db1[i][db1[i]=='6/0']=S[i,5]*(-np.log10((af0[i])*(af[i,5])))
        db1[i][db1[i]=='0/7']=S[i,6]*(-np.log10((af0[i])*(af[i,6])))
        db1[i][db1[i]=='7/0']=S[i,6]*(-np.log10((af0[i])*(af[i,6])))
        db1[i][db1[i]=='0/8']=S[i,7]*(-np.log10((af0[i])*(af[i,7])))
        db1[i][db1[i]=='8/0']=S[i,7]*(-np.log10((af0[i])*(af[i,7])))
        db1[i][db1[i]=='0/9']=S[i,8]*(-np.log10((af0[i])*(af[i,8])))
        db1[i][db1[i]=='9/0']=S[i,8]*(-np.log10((af0[i])*(af[i,8])))
        db1[i][db1[i]=='1/1']=S[i,0]*(-np.log10((af[i,0])*(af[i,0])))
        db1[i][db1[i]=='1/2']=((S[i,0]+S[i,1])/2)*(-np.log10((af[i,0])*(af[i,1])))
        db1[i][db1[i]=='2/1']=((S[i,0]+S[i,1])/2)*(-np.log10((af[i,0])*(af[i,1])))
        db1[i][db1[i]=='1/3']=((S[i,0]+S[i,2])/2)*(-np.log10((af[i,0])*(af[i,2])))
        db1[i][db1[i]=='3/1']=((S[i,0]+S[i,2])/2)*(-np.log10((af[i,0])*(af[i,2])))
        db1[i][db1[i]=='1/4']=((S[i,0]+S[i,3])/2)*(-np.log10((af[i,0])*(af[i,3])))
        db1[i][db1[i]=='4/1']=((S[i,0]+S[i,3])/2)*(-np.log10((af[i,0])*(af[i,3])))
        db1[i][db1[i]=='1/5']=((S[i,0]+S[i,4])/2)*(-np.log10((af[i,0])*(af[i,4])))
        db1[i][db1[i]=='5/1']=((S[i,0]+S[i,4])/2)*(-np.log10((af[i,0])*(af[i,4])))
        db1[i][db1[i]=='1/6']=((S[i,0]+S[i,5])/2)*(-np.log10((af[i,0])*(af[i,5])))
        db1[i][db1[i]=='6/1']=((S[i,0]+S[i,5])/2)*(-np.log10((af[i,0])*(af[i,5])))
        db1[i][db1[i]=='1/7']=((S[i,0]+S[i,6])/2)*(-np.log10((af[i,0])*(af[i,6])))
        db1[i][db1[i]=='7/1']=((S[i,0]+S[i,6])/2)*(-np.log10((af[i,0])*(af[i,6])))
        db1[i][db1[i]=='1/8']=((S[i,0]+S[i,7])/2)*(-np.log10((af[i,0])*(af[i,7])))
        db1[i][db1[i]=='8/1']=((S[i,0]+S[i,7])/2)*(-np.log10((af[i,0])*(af[i,7])))
        db1[i][db1[i]=='1/9']=((S[i,0]+S[i,8])/2)*(-np.log10((af[i,0])*(af[i,8])))
        db1[i][db1[i]=='9/1']=((S[i,0]+S[i,8])/2)*(-np.log10((af[i,0])*(af[i,8])))
        db1[i][db1[i]=='2/2']=S[i,1]*(-np.log10((af[i,1])*(af[i,1])))
        db1[i][db1[i]=='2/3']=((S[i,1]+S[i,2])/2)*(-np.log10((af[i,1])*(af[i,2])))
        db1[i][db1[i]=='3/2']=((S[i,1]+S[i,2])/2)*(-np.log10((af[i,1])*(af[i,2])))
        db1[i][db1[i]=='2/4']=((S[i,1]+S[i,3])/2)*(-np.log10((af[i,1])*(af[i,3])))
        db1[i][db1[i]=='4/2']=((S[i,1]+S[i,3])/2)*(-np.log10((af[i,1])*(af[i,3])))
        db1[i][db1[i]=='2/5']=((S[i,1]+S[i,4])/2)*(-np.log10((af[i,1])*(af[i,4])))
        db1[i][db1[i]=='5/2']=((S[i,1]+S[i,4])/2)*(-np.log10((af[i,1])*(af[i,4])))
        db1[i][db1[i]=='2/6']=((S[i,1]+S[i,5])/2)*(-np.log10((af[i,1])*(af[i,5])))
        db1[i][db1[i]=='6/2']=((S[i,1]+S[i,5])/2)*(-np.log10((af[i,1])*(af[i,5])))
        db1[i][db1[i]=='2/7']=((S[i,1]+S[i,6])/2)*(-np.log10((af[i,1])*(af[i,6])))
        db1[i][db1[i]=='7/2']=((S[i,1]+S[i,6])/2)*(-np.log10((af[i,1])*(af[i,6])))
        db1[i][db1[i]=='2/8']=((S[i,1]+S[i,7])/2)*(-np.log10((af[i,1])*(af[i,7])))
        db1[i][db1[i]=='8/2']=((S[i,1]+S[i,7])/2)*(-np.log10((af[i,1])*(af[i,7])))
        db1[i][db1[i]=='2/9']=((S[i,1]+S[i,8])/2)*(-np.log10((af[i,1])*(af[i,8])))
        db1[i][db1[i]=='9/2']=((S[i,1]+S[i,8])/2)*(-np.log10((af[i,1])*(af[i,8])))
        db1[i][db1[i]=='3/3']=S[i,2]*(-np.log10((af[i,2])*(af[i,2])))
        db1[i][db1[i]=='3/4']=((S[i,2]+S[i,3])/2)*(-np.log10((af[i,2])*(af[i,3])))
        db1[i][db1[i]=='4/3']=((S[i,2]+S[i,3])/2)*(-np.log10((af[i,2])*(af[i,3])))
        db1[i][db1[i]=='3/5']=((S[i,2]+S[i,4])/2)*(-np.log10((af[i,2])*(af[i,4])))
        db1[i][db1[i]=='5/3']=((S[i,2]+S[i,4])/2)*(-np.log10((af[i,2])*(af[i,4])))
        db1[i][db1[i]=='3/6']=((S[i,2]+S[i,5])/2)*(-np.log10((af[i,2])*(af[i,5])))
        db1[i][db1[i]=='6/3']=((S[i,2]+S[i,5])/2)*(-np.log10((af[i,2])*(af[i,5])))
        db1[i][db1[i]=='3/7']=((S[i,2]+S[i,6])/2)*(-np.log10((af[i,2])*(af[i,6])))
        db1[i][db1[i]=='7/3']=((S[i,2]+S[i,6])/2)*(-np.log10((af[i,2])*(af[i,6])))
        db1[i][db1[i]=='3/8']=((S[i,2]+S[i,7])/2)*(-np.log10((af[i,2])*(af[i,7])))
        db1[i][db1[i]=='8/3']=((S[i,2]+S[i,7])/2)*(-np.log10((af[i,2])*(af[i,7])))
        db1[i][db1[i]=='3/9']=((S[i,2]+S[i,8])/2)*(-np.log10((af[i,2])*(af[i,8])))
        db1[i][db1[i]=='9/3']=((S[i,2]+S[i,8])/2)*(-np.log10((af[i,2])*(af[i,8])))
        db1[i][db1[i]=='4/4']=S[i,3]*(-np.log10((af[i,3])*(af[i,3])))
        db1[i][db1[i]=='4/5']=((S[i,3]+S[i,4])/2)*(-np.log10((af[i,3])*(af[i,4])))
        db1[i][db1[i]=='4/6']=((S[i,3]+S[i,5])/2)*(-np.log10((af[i,3])*(af[i,5])))
        db1[i][db1[i]=='4/7']=((S[i,3]+S[i,6])/2)*(-np.log10((af[i,3])*(af[i,6])))
        db1[i][db1[i]=='4/8']=((S[i,3]+S[i,7])/2)*(-np.log10((af[i,3])*(af[i,7])))
        db1[i][db1[i]=='4/9']=((S[i,3]+S[i,8])/2)*(-np.log10((af[i,3])*(af[i,8])))
        db1[i][db1[i]=='5/5']=S[i,4]*(-np.log10((af[i,4])*(af[i,4])))
        db1[i][db1[i]=='5/6']=((S[i,4]+S[i,5])/2)*(-np.log10((af[i,4])*(af[i,5])))
        db1[i][db1[i]=='5/7']=((S[i,4]+S[i,6])/2)*(-np.log10((af[i,4])*(af[i,6])))
        db1[i][db1[i]=='5/8']=((S[i,4]+S[i,7])/2)*(-np.log10((af[i,4])*(af[i,7])))
        db1[i][db1[i]=='5/9']=((S[i,4]+S[i,8])/2)*(-np.log10((af[i,4])*(af[i,8])))
        db1[i][db1[i]=='6/6']=S[i,5]*(-np.log10((af[i,5])*(af[i,5])))
        db1[i][db1[i]=='6/7']=((S[i,5]+S[i,6])/2)*(-np.log10((af[i,5])*(af[i,6])))
        db1[i][db1[i]=='6/8']=((S[i,5]+S[i,7])/2)*(-np.log10((af[i,5])*(af[i,7])))
        db1[i][db1[i]=='6/9']=((S[i,5]+S[i,8])/2)*(-np.log10((af[i,5])*(af[i,8])))
        db1[i][db1[i]=='7/7']=S[i,6]*(-np.log10((af[i,6])*(af[i,6])))
        db1[i][db1[i]=='7/8']=((S[i,6]+S[i,7])/2)*(-np.log10((af[i,6])*(af[i,7])))
        db1[i][db1[i]=='7/9']=((S[i,6]+S[i,8])/2)*(-np.log10((af[i,6])*(af[i,8])))
        db1[i][db1[i]=='8/8']=S[i,7]*(-np.log10((af[i,7])*(af[i,7])))
        db1[i][db1[i]=='8/9']=((S[i,7]+S[i,8])/2)*(-np.log10((af[i,7])*(af[i,8])))
        db1[i][db1[i]=='9/9']=S[i,8]*(-np.log10((af[i,8])*(af[i,8])))
        out1.append(db1[i])
    out1=np.array(out1)
    out1=np.nansum(nan_if(out1,'0.0'),axis=0)
    gg = np.array([gene]*len(samples_header))
    U = np.vstack((samples_header,out1,gg)).T
    return U
#%%

if (np.isnan(scores).sum()) < (scores.shape[0]): #compute metascores if at least 1 variant
    U = score_db(samples,scores,af0,af)
    np.savetxt('./'+gene+'_'+CADD+'_matrix',U, fmt='%s', delimiter='\t')
