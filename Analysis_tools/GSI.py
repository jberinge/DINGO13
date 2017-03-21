# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import os
import numpy as np
from scipy.optimize import curve_fit
import pdb

def sma(data,window):
    
    weights = np.repeat(1.0, window)/window
    smas = np.convolve(data, weights, 'valid')
    return smas

def GSI(data_array,Ta_lo,Ta_hi,DL_lo,DL_hi):

    Ta_array = data_array[:,0]
    DL_array = data_array[:,1]
    VWC_array = data_array[:,2]
    snow_array = data_array[:,3]
    
    Ta_ind_array = np.where(Ta_array < Ta_lo, 0, np.where(Ta_array > Ta_hi, 1, (Ta_array - Ta_lo) / (Ta_hi - Ta_lo)))
    
    DL_ind_array = np.where(DL_array < DL_lo, 0, np.where(DL_array > DL_hi, 1, (DL_array - DL_lo) / (DL_hi - DL_lo)))

#    VWC_ind_array = np.where(VWC_array < VWC_lo, 0, np.where(VWC_array > VWC_hi, 1, (VWC_array - VWC_lo) / (VWC_hi - VWC_lo)))
    
    GSI_array = Ta_ind_array * DL_ind_array * snow_array #* VWC_ind_array 
    
    GSI_array_ext = np.tile(GSI_array, 3)
    
    GSI_array_ext_run = sma(GSI_array_ext, 21)

    GSI_array_run = GSI_array_ext_run[721:1452]
    
    return GSI_array_run

file_in=os.path.join('/home/imchugh/Temp','GSI_Dargo.xlsx')
df=pd.read_excel(file_in)

df['Snow'] = df['Snow'].astype(np.float64)
in_array = np.array(df[['Ta','Daylength','VWC','Snow']][df.Aopt_norm!=np.nan])
out_array = np.array(df['Aopt_norm'][df.Aopt_norm!=np.nan])

test_rslt, test_cov = curve_fit(GSI, in_array, out_array, p0=[-2,5,8,10])

