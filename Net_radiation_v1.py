import numpy as np
import pandas as pd
import datetime as dt
from scipy import stats


#Import custom code modules required for adv processing
import Gap_Fill_climatology_v2 as gap_fill_climatology

###### User-set options ######


def Fn_gapfill(FLUXDataframe,VarToProcess):
    
    print "Start Fn gapfilling"   
                         
    #Calculate Corr net radiation which is best guess Fn from other steams of radiation    

    #Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
    construct_label=str(VarToProcess+"_Con")
    FLUXDataframe[construct_label]=np.nan
    #add a column for the constructed QC flag
    #This will be 1 if valid data from the tower else 99
    construct_flag_label=str(VarToProcess+"_Con_QCFlag")
    FLUXDataframe[construct_flag_label]=np.nan
    #Create a series that is just the correlated output of CABLE adjusted to Tower
    corr_label=str(VarToProcess+"_Corr")
    FLUXDataframe[corr_label]=np.nan	
    FLUXDataframe['Fn_temp']=np.nan
    #Calculate Fn and assign to temporary variable
    #Do only where all 4 streams are available
    FLUXDataframe['Fn_temp']=(FLUXDataframe['Fsd_Con']-FLUXDataframe['Fsu_Con'])+(FLUXDataframe['Fld_Con']-FLUXDataframe['Flu_Con'])
    
    #Start by filling with valid  values from tower
    FLUXDataframe[construct_flag_label][FLUXDataframe[VarToProcess].notnull()]=1
    FLUXDataframe[construct_label][FLUXDataframe[VarToProcess].notnull()]=FLUXDataframe[VarToProcess]
    #Fill with series calculated here
    FLUXDataframe[construct_flag_label][FLUXDataframe[construct_flag_label].isnull()]=99
    
    FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==99]=FLUXDataframe['Fn_temp']
    FLUXDataframe[corr_label]=FLUXDataframe['Fn_temp']
  
    #Delete temp variable
    del FLUXDataframe['Fn_temp']
  
    ###########################################
    # Call function to do climatology gap fill#
    ###########################################
    #Do climatological gap filling if still missing data
    FLUXDataframe=gap_fill_climatology.climatology_monthly_diurnal(FLUXDataframe,VarToProcess)    

    print '\n Finished Fn gapfilling returning DF...'    

    return FLUXDataframe