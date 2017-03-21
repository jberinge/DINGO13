# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 05:02:50 2015

@author: imchugh
"""
import pdb
import pandas as pd
import datetime as dt
import os
import csv

def read_MODIS_ORNL_cutouts(file_in):
    
    # User configs
    pad_days = False    
    write_to_csv = True
    file_out_path = '/media/Data/Dropbox/Data_sites non flux/Site data plots and non-flux/Sites/Whroo/Data/LAI/MODIS'
    
    names_list=['row_ID','land_product_code','acq_date','ctr_pt_coords','proc_date','data_set','data_value']
    scalars_list=[0.1,0.01]
    d = {k: [] for k in range(7)}
    
    # Create dictionary containing all data    
    with open(file_in,'r') as infile:
        for line in infile:
            line_list=line.split(',')
            for i,string in enumerate(line_list):
                d[i].append(string)
    for key in d.keys():
        d[names_list[key]]=d.pop(key)
    
    # Put metadata into separate dictionary
    (lpc,)=set(d['land_product_code'])
    (cpc,)=set(d['ctr_pt_coords'])
    header_string = 'land_product_code: '+lpc+'; centre_point_coordinates: '+cpc
    name_string=(file_in.split('/')[-1]).split('.')[0]+'_output.csv'
    file_out_name=os.path.join(file_out_path,name_string)

    # Create new dictionary with individual datasets separated
    sets_list=list(set(d['data_set']))
    split_d={dataset: [] for dataset in sets_list}
    acq_dates=list(set(d['acq_date']))
    acq_dates.sort()
    split_d['acq_date']=[dt.datetime.strptime(i[1:],'%Y%j') for i in acq_dates]
    proc_dates=list(set(d['proc_date']))
    proc_dates.sort()
    split_d['proc_date']=[dt.datetime.strptime(i,'%Y%j%H%M%S') for i in proc_dates]    
    
    # Write data to it
    for i in xrange(len(d['data_value'])):
        dataset=d['data_set'][i]
        datavalue=int(d['data_value'][i])
        split_d[dataset].append(datavalue)
        
    # Put data into dataframe (easiest way to pad timestamps)
    df=pd.DataFrame(split_d)
    df['lpc']=lpc
    df.index=df['acq_date']
    df.drop('acq_date',axis=1,inplace=True)
    new_index=pd.date_range(df.index[0],df.index[-1],freq='8D')
    pdb.set_trace()
    df=df.reindex(new_index)
    if pad_days:
        new_index=pd.date_range(df.index[0],df.index[-1],freq='D')
        df=df.reindex(new_index)
    
    # Output data to .csv    
    if write_to_csv:
        pathname = os.path.join(file_out_path,file_out_name)
        df.to_csv(pathname,sep=',',index_label = 'Acquisition date')
#        f = open(pathname,'w')
#        csvWriter = csv.writer(f, delimiter=',')
#        csvWriter.writerow(list(split_d.keys()))
#        for i in xrange(len(split_d['acq_date'])):
#            csvWriter.writerow([split_d[key][i] for key in split_d.keys()])
#        f.close()
    
    # Return data to function
    print header_string
    return split_d

def combine_MODIS_cutout_csv(file1,file2):
    
    df1=pd.read_csv(file1,parse_dates=['Acquisition date'],index_col=['Acquisition date'])
    (df1_product_string, ) = set(df1.lpc)
    df1.columns=[name+'_'+df1_product_string for name in df1.columns]    
    df2=pd.read_csv(file2,parse_dates=['Acquisition date'],index_col=['Acquisition date'])
    (df2_product_string, ) = set(df2.lpc)
    df2.columns=[name+'_'+df2_product_string for name in df2.columns]
    pdb.set_trace()
                
            
            