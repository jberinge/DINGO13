# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:52:53 2014

@author: imchugh
"""
import datetime as dt
import pandas as pd
import numpy as np
import pdb

def diurnal_mean(df,mode='mean'):
    if isinstance(df,pd.Series):
        cols=[df.name]
    else:
        cols=df.columns
    if mode=='mean':
        diurnal_df=df.groupby([lambda x: x.hour, lambda y: y.minute]).mean()
    elif mode=='std':
        diurnal_df=df.groupby([lambda x: x.hour, lambda y: y.minute]).std()
    diurnal_df.index=np.arange(48.0)/2    
    diurnal_df.columns=cols
    return diurnal_df

def daily_mean(df,expr):
    if isinstance(df,pd.Series):
        cols=[df.name]
    else:
        cols=df.columns
    daily_df=df.groupby([lambda x: x.year, lambda y: y.dayofyear]).mean()*eval(expr)
    daily_df=daily_df.reset_index()
    daily_df.index=(daily_df['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                    daily_df['level_1'].apply(lambda x: dt.timedelta(int(x)-1)))
    daily_df.drop(['level_0','level_1'],axis=1,inplace=True)
    daily_df.columns=cols
    return daily_df
    
def monthly_mean_by_year(df,mode='mean'):
    if isinstance(df,pd.Series):
        cols=[df.name]
    else:
        cols=df.columns
    if mode=='mean':
        monthly_df=df.groupby([lambda x: x.year, lambda y: y.month]).mean()
    elif mode=='sum':
        monthly_df=df.groupby([lambda x: x.year, lambda y: y.month]).sum()
    elif mode=='std':
        monthly_df=df.groupby([lambda x: x.year, lambda y: y.month]).std()
    else:
        print 'Unknown mode, returning...'
        return
    monthly_df=monthly_df.reset_index()
    monthly_df.index=pd.date_range(dt.datetime(monthly_df.level_0.iloc[0],monthly_df.level_1.iloc[0],1),
                                   periods=len(monthly_df),freq='m')
    monthly_df.drop(['level_0','level_1'],axis=1,inplace=True)
    monthly_df.columns=cols
    return monthly_df

def monthly_mean_all(df,mode='mean'):
    if isinstance(df,pd.Series):
        cols=[df.name]
    else:
        cols=df.columns
    if mode=='mean':    
        monthly_df=df.groupby([lambda x: x.month]).mean()
    elif mode=='sum':
        monthly_df=df.groupby([lambda x: x.month]).sum()
    elif mode=='std':
        monthly_df=df.groupby([lambda x: x.month]).std()
    else:
        print 'Unknown mode, returning...'
        return
    monthly_df=monthly_df.reindex(df.index.month)
    monthly_df.index=df.index
    monthly_df.columns=cols
    return monthly_df
        
def cuml_sum(df):
    if isinstance(df,pd.Series):
        cols=[df.name]
        df=pd.DataFrame(df)
    else:
        cols=df.columns
    cml_df=pd.DataFrame()
    for i in cols:
        cml_df[i]=df[i].cumsum()
    cml_df.columns=[i+'_cml' for i in cml_df.columns]
    return cml_df
    
def run_mean(df,window,center):
    if isinstance(df,pd.Series):
        cols=[df.name]
        df=pd.DataFrame(df)
    else:
        cols=df.columns
    rm_df=pd.DataFrame()
    for i in cols:
        rm_df[i]=pd.rolling_mean(df[i],window,center=center)
    rm_df.columns=[i+'_rm' for i in rm_df.columns]
    return rm_df