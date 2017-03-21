import os
import pandas as pd
import numpy as np
import datetime as dt
from scipy.optimize import curve_fit


###### Constants ######

time_step_noct=5            #Number of steps in days between successive windows for nocturnal fitting
avg_window_noct=30          #Width of window for nocturnal fitting
time_step_day=5             #Number of steps in days between successive windows for daytime fitting
avg_window_day=30           #Width of window for daytime fitting
light_threshold=5           #Minimum light level for fitting
VPD_threshold=1             #Estimated (literature) threshold for VPD (kPa) stomatal response
data_avail_threshold=50     #Percentage of data required in each year for fitting of Eo
data_period='30Min'         #Set the measurement period
temp_spread=5               #Minimum acceptable temperature (C) spread across sample for fitting
min_records=20              #Minimum number of acceptable records for fitting
ind_years=False             #Set to true to fit temperature sensitivity parameter individually for each year
filter_Fc=True              #Set to true if the dataset contains a QC flag for the carbon fluxes
filter_val=1                #If filtering is set to true, this value specifies the QC flag for acceptable data

###### Variable names ######

radName='Fsd_Con'           #Input name for solar radiation
radType='S'                 #Input type of measurement ('S' for full short-wave solar radiation specturm; 'PAR' for 400-700nm)
radUnits='Wm-2'             #Input units for incoming radiation ('Wm-2' or 'umol m-2 s-1')
tempName='Ta_Con'           #Input name for temperature
VPDName='VPD_Con'           #Input name for vapour pressure deficit
VPDUnits='kPa'              #Input units for VPD
CfluxName='Fc'              #Input name for carbon flux data
CfluxUnits='umolCO2 m-2 s-1' #Input units for carbon flux data ('mgCO2 m-2 s-1' or 'umolCO2 m-2 s-1')
filterName='Fc_Con_QCFlag'  #Input name for filter variable in dataframe

###### Parameters ######

Tref=10                      #Temperature reference for basal respiration (celsius)
D0=1                         #Specify units in kPa

###### Data path and name ######

path='C:\Temp'
filename='Advanced_processed_data_HowardSprings_v11a.csv'

###### Functions ######

def convert_to_PAR(df_rad):
    return df_rad*0.45
    
def convert_rad_to_umol(df_rad):
    return df_rad*4.6

def convert_Cflux_to_umol(df_Cflux):
    return df_Cflux*(10**3/44)

def TRF_Eo(local_df,rb,Eo):
    return rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df+46.02)))

def TRF_rb(local_df,rb):
    return rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df+46.02)))

def GPP(local_df,alpha,Aopt,):
    return (alpha*local_df)/(1-(local_df/2000)+(alpha*local_df/Aopt))

def LRF(local_df,alpha,Aopt_0,k,rb):
    Aopt=Aopt_0*np.exp(-k*(local_df[VPDName]-D0))
    #Aopt=1/(1+(local_df[VPDName]-k)/k)    
    index=np.where(local_df[VPDName]<=D0)[0]
    Aopt[index]=Aopt_0
    NEE=(alpha*local_df[radName])/(1-(local_df[radName]/2000)+(alpha*local_df[radName]/Aopt))
    Reco=rb*np.exp(Eo*(1/(Tref+46.02)-1/(local_df[tempName]+46.02)))
    return NEE+Reco

# Do the optimisation (handle errors) then return parameters
def get_data(local_df,date_list):
            
    for i in date_list:
                
        # Slice data from the complete dataframe
        sub_df=local_df.ix[i-dt.timedelta(days=avg_window_noct/2)+dt.timedelta(hours=t_del):
                           i+dt.timedelta(days=avg_window_noct/2)-dt.timedelta(hours=t_del)].dropna(axis=0,how='any')
        
        # If either too few data or temperature range is less than 5C, abort optimisation
        if len(sub_df) >= min_records and sub_df[tempName].max() - sub_df[tempName].min() >= temp_spread:
            
            global Eo
            Eo=params_df['Eo'].ix[i]        
            
            # Try optimisation - if causes error return nan array    
            if nocturnal_fit==True:                  
                try:
                    params_df['rb_noct'][i]=curve_fit(TRF_rb,sub_df[tempName],sub_df[CfluxName],p0=1)[0]
                except RuntimeError:
                    params_df['rb_noct'][i]=np.nan
            else:
                try:
                    a=curve_fit(LRF,sub_df[[radName,tempName,VPDName]],sub_df[CfluxName],p0=[-0.1,-10,1,1])[0]
                except RuntimeError:
                    a=[np.nan,np.nan,np.nan,np.nan]
                params_df['alpha'][i]=a[0]
                params_df['Aopt'][i]=a[1]
                params_df['k'][i]=a[2]
                params_df['rb_day'][i]=a[3]                        

            
#------------------------#                                    
###### Main program ######

def Lasslop(original_df):
    
    #Make a copy of the original dataframe then use that because the DF is later modified
    df=original_df.copy()
    
    ###### Housekeeping ######    
            
    # Globals
    global t_del, params_d, params_df, nocturnal_fit, radName
    
    # Trim for moving window (for even-numbered averaging windows)
    if avg_window_noct%2==False:
        t_del=12
    else:
        t_del=0
    

    ###### Filtering ######
    if filter_Fc==True:    
        df[CfluxName]=np.where(df[filterName]==filter_val,df[CfluxName],np.nan)
    
    
    ###### Conversions ######
    
    # Convert the radiation data
    if radType=='S':
        df[radName]=convert_to_PAR(df[radName])
    if radUnits=='Wm-2':
        df[radName]=convert_rad_to_umol(df[radName])
        
    # Convert the carbon flux data
    if CfluxUnits=='mgCO2 m-2 s-1':
        df[CfluxName]=convert_Cflux_to_umol(df[CfluxName])
    
        
    ###### Create working arrays / dataframes ######
    
    # Date arrays for all days in time series and for centre points of moving windows
    date_array=np.array([dt.datetime.strftime(i,'%Y-%m-%d') for i in df.asfreq('D').index]) # Str array of all day dates in df
    valid_noct_dates=pd.to_datetime(pd.Series(date_array[avg_window_noct/2+1:len(date_array)-(avg_window_noct/2+1):time_step_noct])) # Datetime index of daytime dates to be sampled
    valid_day_dates=pd.to_datetime(pd.Series(date_array[avg_window_day/2+1:len(date_array)-(avg_window_day/2+1):time_step_day])) # Datetime index of nocturnal dates to be sampled
    
    # Yearly frequency dataframe with number of obs, number of cases and percentage of available data for each year
    years_df=pd.DataFrame({'N_obs':df[CfluxName].groupby([lambda x: x.year]).count()})
    years_df['N_recs']=[366 if i%4 else 365 for i in years_df.index]
    years_df['N_recs']=1440/int(data_period[:2])*years_df['N_recs']
    years_df['Data_avail_pct']=np.int8(np.float64(years_df['N_obs'])/years_df['N_recs']*100)
    
    # Separate day and night data
    #There is a column called day_night which calculates sunrise and sunset each day.  1 for daytime, 2 for 3 hours of evening and 3 for night. 
    #noct_df=df[[tempName,CfluxName]][df['day_night']!=1].dropna(axis=0,how='any')
    #day_df=df[[tempName,CfluxName,radName,VPDName]][df['day_night']==1].dropna(axis=0,how='any')
    
    noct_df=df[[tempName,CfluxName]][df[radName]<light_threshold].dropna(axis=0,how='any')
    day_df=df[[tempName,CfluxName,radName,VPDName]][df[radName]>light_threshold].dropna(axis=0,how='any')    

    # Parameters dataframe to contain fit parameters of temperature and light response functions
    params_df=pd.DataFrame({'Eo':np.nan,'rb_noct':np.nan,'rb_noct_interp':np.nan,'rb_day':np.nan,
                            'rb_day_interp':np.nan,'alpha':np.nan,'alpha_interp':np.nan,'Aopt':np.nan,
                            'Aopt_interp':np.nan,'k':np.nan},index=pd.to_datetime(date_array))
    
    # Output dataframe to contain estimates of Re
    Fre_df=pd.DataFrame(index=df.index)
    
    ###### Nocturnal optimisation (Eo) ######
    
    nocturnal_fit=True
    
    # Fix Eo using nocturnal data for entire dataset
    years_df['Eo']=curve_fit(TRF_Eo,noct_df[tempName],noct_df[CfluxName],p0=[10,200])[0][1]
    
    # Fix Eo using nocturnal data for individual years        
    if ind_years==True:
        for i in years_df.iloc[np.where(years_df['Data_avail_pct']>data_avail_threshold)].index:
            years_df['Eo'][i]=curve_fit(TRF_Eo,noct_df[tempName].ix[str(i)],noct_df[CfluxName].ix[str(i)],p0=[10,200])[0][1]
    
    # Fill the parameters dataframe Eo parameter
    for i in years_df.index:
        params_df['Eo'].ix[str(i)]=years_df['Eo'].ix[i]
    
    # Fill the parameters dataframe rb_noct parameter
    get_data(noct_df,valid_noct_dates)        

    
    ######## Daytime optimisation (Aopt,alpha,rb,k) ######
    
    nocturnal_fit=False
    
    # Fill the parameters dataframe alpha, Aopt and rb_day parameters ######
    get_data(day_df,valid_day_dates)
    
    
    ###### Clean up data ######
    
    # Remove data with wrong sign and interpolate
    params_df['alpha']=np.where(params_df['alpha']<0,params_df['alpha'],np.nan)
    params_df['alpha_interp']=params_df['alpha'].interpolate()
    params_df['Aopt']=np.where(params_df['Aopt']<0,params_df['Aopt'],np.nan)
    params_df['Aopt_interp']=params_df['Aopt'].interpolate()
    params_df['k']=np.where(params_df['k']>0,params_df['k'],np.nan)
    params_df['rb_noct']=np.where(params_df['rb_noct']>0,params_df['rb_noct'],np.nan)
    params_df['rb_noct_interp']=params_df['rb_noct'].interpolate()
    params_df['rb_day']=np.where(params_df['rb_day']>0,params_df['rb_day'],np.nan)
    params_df['rb_day_interp']=params_df['rb_day'].interpolate()
    
    
    ###### Calculate and output Re for time series ######
    
    # Calculate Re for time series for day and night estimates
    Fre_df['Fre_noct']=pd.concat(list(TRF_Eo(df[tempName].ix[i:i+dt.timedelta(hours=23.5)],
                                             params_df['rb_noct_interp'][i],params_df['Eo'][i]) for i in params_df.index))
    Fre_df['Fre_day']=pd.concat(list(TRF_Eo(df[tempName].ix[i:i+dt.timedelta(hours=23.5)],
                                            params_df['rb_day_interp'][i],params_df['Eo'][i]) for i in params_df.index))

    Fre_df['GPP_Lasslop']=pd.concat(list(GPP(df[radName].ix[i:i+dt.timedelta(hours=23.5)],
                                            params_df['alpha_interp'][i],params_df['Aopt_interp'][i]) for i in params_df.index))
    
    #Try to delete existing columns in the original DF
    try:
        del df['Fre_noct']
        del df['Fre_day']
        del df['GPP_Lasslop']
    except:
        pass
  
    try:
        del original_df['Fre_noct']
        del original_df['Fre_day']
        del original_df['GPP_Lasslop']
    except:
        pass
    
    #Then join new variables to DF
    ALL_combined=original_df.join(Fre_df,how='left')   
    
    # Force nighttime GPP to 0
    ALL_combined['GPP_Lasslop'][ALL_combined['day_night']!=1] = 0 
    
    print "Finished Lasslop returning dataframe"
    
    return ALL_combined