import numpy as np
import pandas as pd
import datetime as dt
import os
from scipy import stats
import pdb

#Import custom code modules required for adv processing
import Gap_Fill_climatology_v2 as gap_fill_climatology
import gapfill_utilities_v3 as gf

###### User-set options ######

write_to_DF=True

###### Peripheral functions ######
def linear_trans(DF,proc,year):

    if proc==True:
        DF=DF.dropna(how='any',axis=0)
        reg_stats=np.array(stats.linregress(DF['MOD_Refl'],DF['Albedo']))
        if (reg_stats[3]<0.05)&(abs(reg_stats[4]/reg_stats[0])<0.5):
            print str(year)+': data processed, p='+str(round(reg_stats[3],3))
        else:
            reg_stats[0]=1
            reg_stats[1]=0
            print str(year)+': data processed, failed statistical QC, transform not applied... using untransformed MODIS values!'
    else:
        print str(year)+': insufficient data to calculate regression parameters... using untransformed MODIS data!'
        reg_stats=np.array([1,0,0,1,0])
    return reg_stats

def RMSE(Obs_S,pred_S):
    return round(np.sqrt(((Obs_S-pred_S)**2).mean()),3)

###### Main function ######                
def Fsu_gapfill(FLUXDataframe,VarToProcess,fluxfreq):

    rec_length=int(fluxfreq[:2])
    
    print "Starting Fsu gapfilling..."
    print "Calculate daily averages for obs and add to MODIS_DF"
    
    ###### Calculate daily averages for obs and add to MODIS_DF ######
    
    #Output daily summed incoming and outgoing solar radiation from obs
    ObsDailySum_DF=pd.DataFrame({'Fsd_mean':FLUXDataframe['Fsd'].groupby([lambda x: x.year,lambda y: y.dayofyear]).mean(),
                            'Fsd_count':FLUXDataframe['Fsd'].groupby([lambda x: x.year,lambda y: y.dayofyear]).count(),
                            VarToProcess+'_mean':FLUXDataframe[VarToProcess].groupby([lambda x: x.year,lambda y: y.dayofyear]).mean(),
                            VarToProcess+'_count':FLUXDataframe[VarToProcess].groupby([lambda x: x.year,lambda y: y.dayofyear]).count(),
                            'MOD_Refl':FLUXDataframe['sur_refl_b07_new_interp'].groupby([lambda x: x.year,lambda y: y.dayofyear]).mean()}).reset_index()
    
    # Get datetimeindex from groupby variables
    ObsDailySum_DF.index=(ObsDailySum_DF['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                          ObsDailySum_DF['level_1'].apply(lambda x: dt.timedelta(int(x)-1)))    
    ObsDailySum_DF=ObsDailySum_DF.drop(['level_0','level_1'],axis=1)
    
    # Throw out data that don't have all records then calculate albedo
    ObsDailySum_DF['Fsd_mean']=np.where(ObsDailySum_DF['Fsd_count']==1440/rec_length,ObsDailySum_DF['Fsd_mean'],np.nan)
    ObsDailySum_DF[VarToProcess+'_mean']=np.where(ObsDailySum_DF[VarToProcess+'_count']==1440/rec_length,ObsDailySum_DF[VarToProcess+'_mean'],np.nan)
    ObsDailySum_DF['Albedo']=ObsDailySum_DF[VarToProcess+'_mean']/ObsDailySum_DF['Fsd_mean']

    # Calculate climatology for observational data (use interval specified above)
    climatol_S=gf.fill_daily(ObsDailySum_DF['Albedo'],'month')
    ObsDailySum_DF[climatol_S.name]=climatol_S

    # Linear transform of MODIS albedo data - year by year #
    
    #Check the number of years with sufficient data for regression
    temp_DF=ObsDailySum_DF.groupby(lambda x: x.year).count()
    temp_DF['Test']=np.logical_and(temp_DF['MOD_Refl'].map(lambda x: x>10),temp_DF['Albedo'].map(lambda x: x>10))
    
    #Pass data to the regression function, return regression parameters and statistics and add to DF
    stats_DF=pd.DataFrame(np.vstack([linear_trans(ObsDailySum_DF[['MOD_Refl','Albedo']].ix[str(i)],temp_DF['Test'].ix[i],i) for i in temp_DF.index]),
                          index=temp_DF.index,
                          columns=['slope','intcpt','r_sq','p_val','SE'])
    
    # Check if any years pass basic QC; if so, use regression paarmeters for those years to fill others
    # Linear transform of MODIS data                     
    stats_DF['QC_pass']=temp_DF['Test']
    stats_DF['QC_pass']=(stats_DF['p_val']<0.05)&(stats_DF['SE']/stats_DF['slope']<0.5)&(stats_DF['QC_pass']==True)
    if np.any(stats_DF['QC_pass']):
        MODIS_OK=True
        print 'The following years passed basic QC (minimum data threshold and significance of fit): '+', '.join(map(str, stats_DF[stats_DF['QC_pass']].index)) 
        mean_slope=stats_DF['slope'][stats_DF['QC_pass']].mean()
        mean_intcpt=stats_DF['intcpt'][stats_DF['QC_pass']].mean()
        stats_DF['slope']=np.where(stats_DF['QC_pass'],stats_DF['slope'],mean_slope)
        stats_DF['intcpt']=np.where(stats_DF['QC_pass'],stats_DF['intcpt'],mean_intcpt)
        transform_DF=pd.concat([pd.DataFrame({'MOD_Trans':ObsDailySum_DF['MOD_Refl'].ix[str(i)]*stats_DF['slope'].ix[i]+stats_DF['intcpt'].ix[i]}) 
                                for i in stats_DF.index])
        ObsDailySum_DF['MOD_Trans']=transform_DF['MOD_Trans']
    else:
        MODIS_OK=False
            
    # Find best performing algorithm and stack the data in order of best to worst
    if MODIS_OK:
        RMSE_DF=gf.RMSE(ObsDailySum_DF['Albedo'],ObsDailySum_DF[['MOD_Trans',climatol_S.name]])
        best_var=RMSE_DF.index[0]
        ObsDailySum_DF['combined_est']=np.nan
        for i in RMSE_DF.index:
            ObsDailySum_DF['combined_est']=np.where(np.isnan(ObsDailySum_DF['combined_est']),ObsDailySum_DF[i],ObsDailySum_DF['combined_est'])
        print 'Variable with lowest RMSE = '+best_var
        if best_var==climatol_S.name:
            print 'MODIS RMSE higher than climatological average!'
            print 'Using climatology'
        else:
            print 'Using MODIS data'
    else:
        ObsDailySum_DF['combined_est']=ObsDailySum_DF[climatol_S.name]
        print 'No years with acceptable fit quality - using climatology'        


    # Apply gap-filled MODIS-derived albedo to incoming shortwave radiation
    
    #Reset both dataframes to hierarchical index (year,day of year)
    ObsDailySum_DF['DT']=ObsDailySum_DF.index
    ObsDailySum_DF.index=pd.MultiIndex.from_tuples(zip(*[np.array(ObsDailySum_DF['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%Y')))),
                                                np.array(ObsDailySum_DF['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%j'))))]),names=['Year', 'DOY'])
    
    FLUXDataframe['DT']=FLUXDataframe.index
    FLUXDataframe.index=pd.MultiIndex.from_tuples(zip(*[np.array(FLUXDataframe['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%Y')))),
                                                        np.array(FLUXDataframe['DT'].apply(lambda x: int(dt.datetime.strftime(x,'%j'))))]),names=['Year', 'DOY'])
    
    #Calculate outgoing shortwave radiation using MODIS derived albedo and incoming shortwave radiation
    output_Fsu=pd.concat([pd.DataFrame(FLUXDataframe['Fsd_Corr'].ix[i]*ObsDailySum_DF['combined_est'].ix[i]).set_index(FLUXDataframe['DT'].ix[i]) for i in ObsDailySum_DF.index])
    
    #Reset the main dataframe index to original datetimeindex
    FLUXDataframe.index=FLUXDataframe['DT']
    FLUXDataframe.drop('DT',axis=1)
    
    ###########################################    
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
        FLUXDataframe[construct_flag_label][FLUXDataframe[construct_flag_label].isnull()]=99
        
        FLUXDataframe[corr_label]=output_Fsu
        FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==99]=FLUXDataframe[corr_label]
                        
        #FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==99]=output_Fsu
        #FLUXDataframe[corr_label]=output_Fsu
	
    ###########################################
    # Call function to do climatology gap fill#
    ###########################################
    #Do climatological gap filling if stillmissing data
    FLUXDataframe=gap_fill_climatology.climatology_monthly_diurnal(FLUXDataframe,VarToProcess)    

    print '\n Finished Fsu gapfilling returning DF...'
    return FLUXDataframe
