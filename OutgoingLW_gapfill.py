import numpy as np
import pandas as pd
import datetime as dt
from scipy import stats


#Import custom code modules required for adv processing
import Gap_Fill_climatology_v2 as gap_fill_climatology

###### User-set options ######

num_avg=10
write_to_DF=True

###### Globals ######

LW_DF_global=pd.DataFrame()
DailyLW_global=pd.DataFrame()
VarToProcess=None

###### Functions ######

def find_nearest(dateObj):
    NearDays_DF=(abs(DailyLW_global['SRT_Day_mean']-DailyLW_global['MODISDay_pred'][dateObj])+
                 abs(DailyLW_global['SRT_Night_mean']-DailyLW_global['MODISNight_pred'][dateObj]))
    NearDays_DF.sort()
    NearDays_DF=NearDays_DF[:num_avg]    
    NearData_DF=pd.concat([LW_DF_global[VarToProcess][dt.datetime.strftime(j,'%Y-%m-%d')] for j in NearDays_DF.index])
    temp_DF=NearData_DF.groupby([lambda x: x.hour,lambda y: y.minute]).mean()
    temp_DF.index=pd.date_range(dateObj,periods=1440/rec_length,freq=freq_str)
    return temp_DF

def Flu_gapfill(FLUXDataframe,variable_to_fill,fluxfreq):
    
    ###### Set constants ######
    
    # Globals
    global VarToProcess
    VarToProcess=variable_to_fill
    global rec_length
    rec_length=int(fluxfreq[:2])
    global freq_str
    freq_str=fluxfreq
    
    # Locals
    min_propn_data = 0.8
    print "Start Flu gapfilling"
      

    #Check to see if at least 50% of the longwave data is present to make this work.  Otherwise use simple approach below.                         
    Var_total=len(FLUXDataframe[VarToProcess])
    Var_notnull=FLUXDataframe[VarToProcess][FLUXDataframe[VarToProcess].notnull()].count()
    if Var_notnull>(Var_total/2):
        #Calculate surface radiative temperature
        S=FLUXDataframe[VarToProcess].apply(lambda x: (x/(5.67*10**-8))**(1/4.0)-273.15)
                
        #Split day and night with threshold of 10Wm-2
        sub_DF=pd.DataFrame({'SRT_Day':np.where(FLUXDataframe['Fsd_Con']>10,S,np.nan),
                            'SRT_Night':np.where(FLUXDataframe['Fsd_Con']<10,S,np.nan),
                            'LST_Day_interp':np.where(FLUXDataframe['Fsd_Con']>10,FLUXDataframe['LST_Day_1km_new_interp'],np.nan),
                            'LST_Night_interp':np.where(FLUXDataframe['Fsd_Con']<10,FLUXDataframe['LST_Night_1km_new_interp'],np.nan)},index=FLUXDataframe.index)
        
        #Output daily summed surface radiative temperature from obs
        by = lambda x: lambda y: getattr(y, x)
        DailyLW_DF=pd.DataFrame({'LW_count':S.groupby([by('year'),by('dayofyear')]).count(),
                                'SRT_Day_mean':sub_DF['SRT_Day'].groupby([by('year'),by('dayofyear')]).mean(),
                                'SRT_Night_mean':sub_DF['SRT_Night'].groupby([by('year'),by('dayofyear')]).mean(),
                                'LST_Day_mean':sub_DF['LST_Day_interp'].groupby([by('year'),by('dayofyear')]).mean(),
                                'LST_Night_mean':sub_DF['LST_Night_interp'].groupby([by('year'),by('dayofyear')]).mean()}).reset_index()   

       
        #Generate dates from group names (year and day of year)
        DailyLW_DF.index=(DailyLW_DF['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                          DailyLW_DF['level_1'].apply(lambda x: dt.timedelta(int(x)-1)))        
        DailyLW_DF=DailyLW_DF.drop(['level_0','level_1'],axis=1)
               
        #Throw out unwanted data
        DailyLW_DF['SRT_Day_mean']=np.where(DailyLW_DF['LW_count']>=int(min_propn_data*(1440/rec_length)),DailyLW_DF['SRT_Day_mean'],np.nan)
        DailyLW_DF['SRT_Night_mean']=np.where(DailyLW_DF['LW_count']>=int(min_propn_data*(1440/rec_length)),DailyLW_DF['SRT_Night_mean'],np.nan)
        DailyLW_DF=DailyLW_DF.drop(['LW_count'],axis=1)
    
        #Merge obs and MODIS
        temp_DF=pd.merge(DailyLW_DF[['SRT_Day_mean']],FLUXDataframe[['LST_Day_1km_new']],left_index=True,right_index=True).dropna(axis=0,how='any')
        rslt_day=stats.linregress(temp_DF['LST_Day_1km_new'],temp_DF['SRT_Day_mean'])
        DailyLW_DF['MODISDay_pred']=DailyLW_DF['LST_Day_mean']*rslt_day[0]+rslt_day[1]
        temp_DF=pd.merge(DailyLW_DF[['SRT_Night_mean']],FLUXDataframe[['LST_Night_1km_new']],left_index=True,right_index=True).dropna(axis=0,how='any')
        rslt_night=stats.linregress(temp_DF['LST_Night_1km_new'],temp_DF['SRT_Night_mean'])
        DailyLW_DF['MODISNight_pred']=DailyLW_DF['LST_Night_mean']*rslt_night[0]+rslt_night[1]
       
        ###### Create globals to avoid iterative pass of DF to search function (slooooooow!) ######
        global LW_DF_global
        LW_DF_global=FLUXDataframe
        global DailyLW_global                              
        DailyLW_global=DailyLW_DF
        
        output_LW=pd.concat([find_nearest(dateObj) for dateObj in DailyLW_DF.index])
        
        output_LW.reindex(FLUXDataframe.index)
        
        construct_flag_number=99
        
    else:
        #if not enough data then do simple stefan Boltzman on Ta. Assume emmisvity 0.95
        #Use Ta as that has best correlation and gap fill
        #Ta convert to K.  Output LW in Wm-2
        output_LW=((FLUXDataframe['Ta_Con']+273.15)**4)*0.95*5.6703*10**-8
        construct_flag_number=93
        
    ############################################    
    # Generate output and return new variables
    ############################################
    if write_to_DF==True:
   	
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
        
        #Start by filling with valid  values from tower
        FLUXDataframe[construct_flag_label][FLUXDataframe[VarToProcess].notnull()]=1
        FLUXDataframe[construct_label][FLUXDataframe[VarToProcess].notnull()]=FLUXDataframe[VarToProcess]
        #Fill with series calculated here
        FLUXDataframe[construct_flag_label][FLUXDataframe[construct_flag_label].isnull()]=construct_flag_number
        
        FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==construct_flag_number]=output_LW
        FLUXDataframe[corr_label]=output_LW
  
    ###########################################
    # Call function to do climatology gap fill#
    ###########################################
    #Do climatological gap filling if stillmissing data
    FLUXDataframe=gap_fill_climatology.climatology_monthly_diurnal(FLUXDataframe,VarToProcess)    

    print '\n Finished Flu gapfilling returning DF...'    

    return FLUXDataframe
