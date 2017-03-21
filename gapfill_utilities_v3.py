import pandas as pd
import numpy as np
import datetime as dt
import pdb

# Requires single pandas series and string for climatology averaging interval
def fill_daily(S1,groupby_string):
    
    # Store the variable name for later
    var_name=S1.name
    
    # Create a dataframe, save a copy of the datetime index in it and separately
    DF=pd.DataFrame(S1)
    DF['DT']=DF.index
    DF_index=DF.index
    
    # Set dictionary to retrieve string for climatology averaging interval
    groupby_dict={'dayofyear':'%j','week':'%W','month':'%m'}
    
    # Create a series of clim averages over whole period of record by groupby_string 
    if groupby_string=='dayofyear':
	DF_temp=DF.groupby(lambda x: x.dayofyear).mean()
    elif groupby_string=='week':
	DF_temp=DF.groupby(lambda x: x.week).mean()
	DF_temp=pd.concat([DF_temp[52:],DF_temp[:]]).reset_index()
    else:
	DF_temp=DF.groupby(lambda x: x.month).mean()
    DF_temp[groupby_string]=DF_temp.index
            
    # Create the relevant interval as a variable in the DF (eg. month, week, DOY)
    DF[groupby_string]=((DF['DT'].apply(lambda x:
	                    int(dt.datetime.strftime(x,groupby_dict[groupby_string])))))
     
    # Merge the series on the groupby_string, reset index (obliterated on merge), sort then reinstate full original DF index
    DF=pd.merge(DF,DF_temp,on=groupby_string)
    DF.index=DF['DT'] 
    DF=DF.sort_index(axis=0)
    DF=DF.reindex(DF_index)
    
    # Fill the missing data using the climatology
    DF[var_name+'_climatol']=DF[var_name+'_y']
    DF[var_name+'_climatol']=DF[var_name+'_climatol'].interpolate()
    DF[var_name+'_climatol']=DF[var_name+'_climatol'].fillna(method='bfill')
    
    return DF[var_name+'_climatol']
                                                                                  
def RMSE(obs_S,est_df):
    
    RMSE_df=pd.DataFrame(index=est_df.columns,columns=['RMSE'])
    
    for i in RMSE_df.index:
        RMSE_df.ix[i]=round(np.sqrt(((obs_S-est_df[i])**2).mean()),3)
    RMSE_df.sort('RMSE',inplace=True)
    
    return RMSE_df
        