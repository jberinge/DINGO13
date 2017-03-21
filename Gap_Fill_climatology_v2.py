import pandas as pd
import numpy as np
import datetime as dt

def climatology_monthly_diurnal(New_combined,variable_to_fill):
    print "Doing climatology gap fill for " + variable_to_fill
    #Variable names defined
    #variable_to_fill passed from function
    variable_Con= variable_to_fill+'_Con'
    variable_Con_QCFlag = variable_to_fill + '_Con_QCFlag'
    variable_Clim  = variable_to_fill + '_Clim'
    
    ###Create a hierarchically-indexed dataframe containing averages for each month, hour and minute
    by = lambda x: lambda y: getattr(y, x)
    #Groupby and select three varibales of intest
    #Add a dummy variable so that the operation produces a Df not series
    climatology_DF1=New_combined[[variable_to_fill,variable_Con]].groupby([by('month'),by('hour'),by('minute')]).mean()
    climatology_DF1.index.names=['MonthC','HourC','MinuteC']
    climatology_DF1.columns=[[variable_Clim,'dummy']] 
    
    ###Create columns matching climatology DF for merging
    New_combined['DT']=New_combined.index
    New_combined['MonthC']=New_combined['DT'].apply(lambda x:int(dt.datetime.strftime(x,'%m')))
    New_combined['HourC']=New_combined['DT'].apply(lambda x:int(dt.datetime.strftime(x,'%H')))
    New_combined['MinuteC']=New_combined['DT'].apply(lambda x:int(dt.datetime.strftime(x,'%M')))
        
    #Before Merge Convert all QC Flags for variable that are missing to nan's
    #Then afterwards if not missing anymore then fill with '97'
    New_combined[variable_Con_QCFlag][New_combined[variable_Con].isnull()]=97
    count_before_fill= len(New_combined[variable_Con][New_combined[variable_Con].isnull()])
    print "Nan count before", count_before_fill
    
    ###Merge, reset index to datetimeindex and sort index (order obliterated by merge)
    #New_combined2=pd.merge(New_combined,climatology_DF,on=['MonthC','HourC','MinuteC'],how='right')
    New_combined2=New_combined.join(climatology_DF1,on=['MonthC','HourC','MinuteC'])
    New_combined2.index=New_combined['DT']
    New_combined3=New_combined2.sort_index(axis=0)
   
    #Now fill _Con with _Clim
    New_combined3[variable_Con][New_combined3[variable_Con_QCFlag]==97.]=New_combined3[variable_Clim]

    #then check to see if there are still any null values and change flag back to nan also
    New_combined3[variable_Con_QCFlag][New_combined3[variable_Con].isnull()]=np.nan
    
    count_after_fill= len(New_combined3[variable_Con][New_combined3[variable_Con].isnull()])
    print "Nan count After", count_after_fill    
    print "Total number of rows in DF ", len(New_combined)
    print "Start date for series ", New_combined.index[0]
    print "End date for series ", New_combined.index[len(New_combined)-1]

    #Clean up
    del New_combined3['MonthC']
    del New_combined3['HourC']
    del New_combined3['MinuteC']
    del New_combined3['dummy']
    del New_combined3[variable_Clim]
    
    #New_combined3.to_csv('E:/My Dropbox/Dropbox/Data_flux_data/Site data processing/HowardSprings/Advanced/test1_'+variable_to_fill+'.csv', sep=',') 
    
    #return DF
    return New_combined3

    
