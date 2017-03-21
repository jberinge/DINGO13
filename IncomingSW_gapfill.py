
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 11:22:36 2013

- Downscaling routine for estimating half-hourly radiation from AWAP daily estimates
- Outputs cloudiness index to be used as input for long wave incoming estimation routine
- Needs work - at Howard requires an arbitrary subsetting to get the extinction coefficient
               optimisation to work (line 265) due to persistent wet season cloudiness!!!

@author: imchugh
Modified by JB June 2013 to work with Pandas
Further Mods by IM July 2013
Rewritten May 2014                      

Refs:
    
    DiLaura, D. L. (1984), IES Calculation Procedures Committee Recommended
    practice for the calculation of daylight availability, J. Illuminating
    Engineering Soc. of North America, 13(4), 381-392.
    
    Duffie, J. and W. Beckman (1980). Solar Engineering of Thermal Processes. 
    New York, John Wiley and Sons.
    
    Wunderlich, W. (1972), Heat and Mass Transfer between a Water Surface
    and the Atmosphere, Report No 14, Report Publication No. 0-6803,
    Water Resources Research Laboratory, TennesseeValleyAuthority,Division
    of Water Control Planning, Engineering Laboratory, Norris, TN.

"""

#------------------------------------------------------------------------------#
# Imports
#------------------------------------------------------------------------------#

# Python modules
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
import datetime as dt
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['legend.fancybox'] = True
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

# Custom modules
import Timezone_v3a as TimeZone
import gapfill_utilities_v3 as gf

#------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------#

# Scale data using AWAP daily estimates
def daily_scale(obs_daily,Kdown_clr_daily,Kdown_clr_hr):
    scaling_coeff=obs_daily/Kdown_clr_daily
    Kdown_cld_hr=Kdown_clr_hr*scaling_coeff.reshape(len(scaling_coeff),1)
    cloud_factor=1-scaling_coeff
    cloud_factor=np.where(cloud_factor<0,0,cloud_factor)
    return np.ravel(Kdown_cld_hr),cloud_factor

# Create the timestamps for output
def time_out(datetime,num_periods,freq_recs,fmt):
    if fmt==0:
        t_del=0
    elif fmt==1:
        t_del=rec_length
    return pd.date_range(datetime+dt.timedelta(minutes=t_del),periods=num_periods,freq=freq_recs)

# Estimate clear sky radiation and optimise k for site obs data
def Insol_calc(DOY,k):
    
    # For each day calculate equation of time correction, solar noon, declination and TOA radiation
    array_EqofTime=0.17*np.sin(4*np.pi*(DOY-80)/373)-0.129*np.sin(2*np.pi*(DOY-8)/355) # DiLaura (1984)
    array_solar_noon=12+(GMT_zone*15.0-long_decdeg)/360*24-array_EqofTime # Me
    array_decl=np.radians(23.4)*np.sin((DOY+284)/365.0*2*np.pi) # Oke (1987)
    array_TOArad=(1+0.034*np.cos(DOY/365.25*2*np.pi))*1367.0 # Duffie and Beckman (1980)
    
    # Create an hour angle array for each minute of day and each day of year
    array_h=np.tile(np.linspace(0,1439.0/1440*24,num=1440),(len(DOY),1))
    array_h=abs(np.radians((array_solar_noon.reshape(len(DOY),1)-array_h)*15))
    
    # Duplicate declination array for each time of day
    array_decl=np.tile(array_decl,(1440,1)).T

    # Calculate zenith angles
    array_z=np.arccos(np.sin(np.radians(lat_decdeg))*np.sin(array_decl)+
            np.cos(np.radians(lat_decdeg))*np.cos(array_decl)*np.cos(array_h))
    array_z_msk=np.ma.masked_greater_equal(array_z,np.pi/2) # Mask night values    
  
    # Calculate optical air mass term for all valid Z 
    array_m=(np.exp(-1*ALT_m/8343.5)/(np.cos(array_z_msk)+0.15*
            (np.degrees(90-array_z_msk)+3.855)**-1.253)) # Wunderlich (1972)           
    
    # Instantaneous clear sky surface radiation in Wm-2 for each minute of the day
    array_Kdown_clr_mins=array_TOArad.reshape(len(array_TOArad),1)*np.exp(-k*array_m)*np.cos(array_z_msk)
    
    # Aggregate one-minute instantaneous clear sky rad to period average
    array_Kdown_clr_hr=np.empty([len(DOY),1440/rec_length])
    for i in xrange(len(DOY)):
        array_temp=array_Kdown_clr_mins[i][:].reshape(1440/rec_length,rec_length) # Temporary bins
        array_Kdown_clr_hr[i][:]=np.ma.mean(array_temp,axis=1) # Average bin content  
    
    # Aggregate to daily
    array_Kdown_clr_daily=(array_Kdown_clr_hr*(rec_length*60.0/10**6)).sum(axis=1)
        
    if boolOutput==False:
        return array_Kdown_clr_daily # Result for optimisation
    else:
        return array_Kdown_clr_daily,array_Kdown_clr_hr    # Output of final data

def plot_linear_reg(xdata,ydata,reg_slope,reg_int,reg_r2,reg_pvalue,reg_se,xlabel,ylabel,OutputPath_DailyData,Site_ID):
    #Calculate some things
    n_datapoints=len(xdata)
    startdate= xdata.index[0]
    enddate= xdata.index[n_datapoints-1]
    #Set temp variables and lists

    tempx_line=[]
    tempy_line=[]

    #For plotting we want to find the range of variable values across all sites and tower to get entire range to plot
    scale_min= int(min(xdata.min(),ydata.min())-1)
    scale_max= int(max(xdata.max(),ydata.max())+1)
    
    #create series to plot line
    #Need to extract the linear regression stats done bygroup earlier

    #slopetemp, intercepttemp, r_valuetemp, p_valuetemp, std_errtemp = stats.linregress(group[Labels[1]],group[Labels[0]])
    
    for increment in range(scale_min,scale_max):
	    tempx_line.append(increment)
	    tempy_line.append(reg_slope*increment+reg_int)
	    
    ## Could work for later  pd.merge(df, k1_means, left_on='key1', right_index=True	
    #Produce the plot 
    IDx=12345
    plt.plot(xdata, ydata, 'go',tempx_line, tempy_line, ':b' ,label=IDx,linewidth=2) 
    #Set the scale mins and maxs
    plt.xlim(scale_min, scale_max)
    plt.ylim(scale_min, scale_max)  
    #create text for ID and r2 box
    graphtext1=     str('intercept   ' + str("{0:.2f}".format(reg_int) +'\n')
                      + 'slope       ' + str("{0:.2f}".format(reg_slope)) +'\n'
                      + 'r value     ' + str("{0:.2f}".format(reg_r2)) +'\n'
                      + 'p_value     ' + str("{0:.2f}".format(reg_pvalue)) +'\n'
                      + 'std_err     ' + str("{0:.2f}".format(reg_se)) +'\n')  
    #create text for start and end dates
    graphtext2=('Data start date: '+str(startdate)+'\n'
                +'End date: '+str(enddate)+'\n'
                +'Number records: '+str(n_datapoints))   
    units=' MJ m-2 d-1'
    plt.figtext(0.7,0.3,graphtext1, bbox=dict())
    plt.figtext(0.5,0.13,graphtext2, bbox=dict())
    plt.title('Tower vs AWAP solar  - ' +Site_ID +'\n')
    plt.xlabel(xlabel + '('+units+')')
    plt.ylabel(ylabel+ '   ' + '('+units+')')
    #plt.legend(shadow=True, fancybox=True,loc='best')
    
    #Output to PDF using PdfPages a backend for MatPlotLib
    fname_graph=OutputPath_DailyData+'/'+'Linear Plot Tower vs AWAP Solar ' +Site_ID+'.pdf'
    # Create the PdfPages object to which we will save the pages:
    pdf = PdfPages(fname_graph)
    savefig(pdf, format='pdf',facecolor='w', edgecolor='w') # note the format='pdf' argument!
    close()		
    pdf.close() 
    
def RMSE(Obs_S,pred_S):
    return round(np.sqrt(((Obs_S-pred_S)**2).mean()),3)

def linear_reg(x,y):
    temp_DF=pd.DataFrame({'x':x,'y':y}).dropna(how='any',axis=0)
    return stats.linregress(temp_DF['x'],temp_DF['y'])

# Get the AWAP data
def import_AWAP(InputPathName_AWAPData,date_list,avail_data_threshold,FLUXDataframe):
    
    # If file exists, import the data and check that dates overlap (if true, return boolean and file, else just boolean)
    if os.path.exists(InputPathName_AWAPData):
        varNames_AWAP=['#YYYY','MM','DD','solar_exposure_day']
        DF=pd.read_csv(InputPathName_AWAPData,usecols=varNames_AWAP,skiprows=[1],
                       parse_dates=[[0,1,2]],index_col=0,na_values=['-9999'],keep_default_na=True)
        DF.replace(-9999,np.nan,inplace=True)
        DF.columns=['AWAP']
        try:
            AWAP_propn_avail=round((len(DF.ix[date_list[0]:date_list[1]].dropna(how='any',axis=0))/float(len(FLUXDataframe))*100),2)
            if AWAP_propn_avail>avail_data_threshold:
                print "AWAP data availability "+ str(AWAP_propn_avail) + " exceeds minimum threshold of " + str(avail_data_threshold) + " and will be included in the analysis..."
                return DF
            else:
                print "AWAP data availability "+ str(AWAP_propn_avail) + " below minimum threshold of " + str(avail_data_threshold) + " and will be excluded in the analysis..."  
        except:
            print 'Dates in AWAP file do not overlap with dates for gap filling - excluding AWAP'    
            return
    else:
        print 'No AWAP file of specified name in specified path' 
        return
            
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# MAIN CODE
#------------------------------------------------------------------------------#

def Fsd_gapfill(variable_to_fill,myBaseforResults,mypathforAWAPdata,FileName_AWAPData,FLUXDataframe,Site_ID,latitude,longitude,altitude,fluxfreq):
    

    ###### Constants ######
    
    # Globals
    global VarToProcess
    VarToProcess=variable_to_fill
    global lat_decdeg
    lat_decdeg=latitude
    global long_decdeg
    long_decdeg=longitude
    global ALT_m
    ALT_m=float(altitude)
    global rec_length
    rec_length=int(fluxfreq[:2])
    global GMT_zone
    currentdate="2013-06-01"
    GMT_zone,InDstNow = TimeZone.get_timezone_info(lat_decdeg,long_decdeg,currentdate)
    global boolOutput
    boolOutput=False
        
    # Locals
    timestamp_format=1                 # Set record naming format (0, time at beginning of period, _
            		                 # 1, time at end of period, 2, time in centre of period)
    write_to_DF=True                   # Set to true if want data appended to existing DF
    groupby_string='week'              # Frequency of binning for climatology: dayofyear, week, month
    BOM_var='Global_Solar'             # Name of dataframe variable containing Bureau estimates of incoming solar
    default_k=0.2                      # Default value for extinction coefficient if no observational data available for optimisation
    avail_data_threshold=50            # Minimum observational data threshold acceptable for optimisation
    var_list=list()                    # Stores names of variables being used
    optimisation_subset_list=[10,40]   # Confines the optimisation to specified start and end weeks
    
    
    ###### Set I/O ######       
                     
    print "Starting Fsd gap filling"
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/ANCILLARY"):
	os.mkdir(myBaseforResults+"/ANCILLARY")
    mypathforResults=myBaseforResults+"/ANCILLARY"    
    
    #Set output path
    OutputPath_DailyData=mypathforResults

    #Check whether variable is present and whether holds sufficient data
    Fsd_avail=VarToProcess in FLUXDataframe.columns
    if Fsd_avail: 
        Fsd_propn_avail=round(len(FLUXDataframe[VarToProcess].dropna())/float(len(FLUXDataframe[VarToProcess]))*100,1)
        if Fsd_propn_avail>avail_data_threshold:
            var_list.append('Obs')
            print 'Available observational data exceeds minimum data availability threshold - optimising extinction coefficient...'
        else:
            Fsd_avail=False
            print 'Available observational data below minimum data availability threshold - using default extinction coefficient...'
    else:
        print 'No observational data - checking for ancillary data'

    #Check for AWAP ancillary data and import if available
    InputPathName_AWAPData=os.path.join(mypathforAWAPdata,FileName_AWAPData)
    date_list=[FLUXDataframe.index[0].date(),FLUXDataframe.index[-1].date()]   
    AWAP_DF=import_AWAP(InputPathName_AWAPData,date_list,avail_data_threshold,FLUXDataframe)
    if type(AWAP_DF)==pd.DataFrame: 
        AWAP_avail=True
        var_list.append('AWAP')
    else: 
        AWAP_avail=False    

    #Check for BOM ancillary data
    BOM_avail=BOM_var in FLUXDataframe.columns
    if BOM_avail:
        BOM_propn_avail=round(len(FLUXDataframe[BOM_var].dropna())/float(len(FLUXDataframe[BOM_var]))*100,1)
        if BOM_propn_avail>avail_data_threshold:
            var_list.append('BOM')
            print 'BOM data availability exceeds minimum threshold and will be included in the analysis...'
        else:
            BOM_avail=False
            print 'BOM data availability below minimum threshold and will be excluded from the analysis...'
    else:
        print 'No observational data - checking for ancillary data'

    #If no data is available, return the dataframe and inform the user
    if Fsd_avail==False and AWAP_avail==False and BOM_avail==False:
        print 'No observational or ancillary data available - returning dataframe...'
        return FLUXDataframe             

                                                                                
    ###### Calculate daily averages for obs (and BOM if required) and add to daily dataframe ######
    
    # Create an empty dataframe to take the daily data
    DaySum_DF=pd.DataFrame()    
                            
    # If the solar radiation data is available, do daily averages (drop days with missing cases). 
    if Fsd_avail:
        DaySum_DF['Obs_count']=FLUXDataframe[VarToProcess].groupby([lambda x: x.year,lambda y: y.dayofyear]).count()
        DaySum_DF['Obs']=FLUXDataframe[VarToProcess].groupby([lambda x: x.year,lambda y: y.dayofyear]).mean()*0.0864
        DaySum_DF['Obs']=np.where(DaySum_DF['Obs_count']==1440/rec_length,DaySum_DF['Obs'],np.nan)                         
    
    # If BOM data is embedded in file, do daily averages as above
    if BOM_avail:
        DaySum_DF['BOM_count']=FLUXDataframe[BOM_var].groupby([lambda x: x.year,lambda y: y.dayofyear]).count()
        DaySum_DF['BOM']=FLUXDataframe[BOM_var].groupby([lambda x: x.year,lambda y: y.dayofyear]).mean()          
        DaySum_DF['BOM']=np.where(DaySum_DF['BOM_count']==1440/rec_length,DaySum_DF['BOM'],np.nan)  
    
    # If dataframe contains data, generate dates from hierarchical index (year and day of year), otherwise just use AWAP data
    if not len(DaySum_DF)==0:
        DaySum_DF=DaySum_DF.reset_index()
        # Get datetimeindex from groupby variables
        DaySum_DF.index=(DaySum_DF['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                         DaySum_DF['level_1'].apply(lambda x: dt.timedelta(int(x)-1)))
        # Buffer beginning and end datetime of file to account for timestamp naming convention
        DaySum_DF.reindex(pd.date_range(FLUXDataframe.index[0].date()-dt.timedelta(days=1),
                                        FLUXDataframe.index[-1].date()+dt.timedelta(days=1),freq='D'))
        # Add AWAP if available
        if AWAP_avail:
            DaySum_DF=DaySum_DF.join(AWAP_DF,how='left')        
    else:
        DaySum_DF=AWAP_DF
    
    # Add day of year, which is key independent variable
    DaySum_DF['DOY']=[i.timetuple().tm_yday for i in DaySum_DF.index]
                                                                
    # Drop everything except the data and day of year variables
    DaySum_DF=DaySum_DF[['DOY']+var_list]    
    
    #Force DF to have numeric values
    if Fsd_avail: DaySum_DF.Obs = pd.to_numeric(DaySum_DF.Obs, errors='coerce')
    if AWAP_avail: DaySum_DF.AWAP = pd.to_numeric(DaySum_DF.AWAP, errors='coerce')
    DaySum_DF.BOM = pd.to_numeric(DaySum_DF.BOM, errors='coerce')
    
                                                                                                                        

    ###### Do linear transforms, get climatology for available data series, calculate RMSE, rank and then fill according to rank ######

    #Linear transform
    if len(var_list)>1: 
        reg_var_list=var_list[1:]
        reg_results_df=pd.DataFrame(index=reg_var_list,columns=['slope','intcpt','rsq','pval','se'])
        for i in reg_var_list:
            reg_results_df.ix[i]=linear_reg(DaySum_DF[i],DaySum_DF[var_list[0]])
            DaySum_DF[i]=DaySum_DF[i]*reg_results_df.ix[i][0]+reg_results_df.ix[i][1]
            print ('Obs / '+i+' Regression stats: Slope = '+str(round(reg_results_df['slope'].ix[i],3))+
                    ';intercept = '+str(round(reg_results_df['intcpt'].ix[i],3))+'; r2 = '+str(round(reg_results_df['rsq'].ix[i],3)))    
            
    # If there is observational data...
    if Fsd_avail:
                                                                                                                                                               
        # Calculate climatology for observational data (use interval specified above)
        climatol_S=gf.fill_daily(DaySum_DF['Obs'],groupby_string)
        DaySum_DF[climatol_S.name]=climatol_S
        var_list.append(climatol_S.name)
    
        # Compare fill variables (i.e. climatology and external data sources) to observations and rank
        RMSE_df=gf.RMSE(DaySum_DF[var_list[0]],DaySum_DF[var_list[1:]])
        
        # Plot regression results for lowest RMSE variable
        best_var=RMSE_df.index[0]
        if not best_var=='Obs_climatol':
            temp_DF=pd.DataFrame({'x':DaySum_DF[best_var],'y':DaySum_DF[var_list[0]]}).dropna(how='any',axis=0)
            plot_linear_reg(temp_DF['x'],temp_DF['y'],reg_results_df['slope'].ix[best_var],reg_results_df['intcpt'].ix[best_var],
                            reg_results_df['rsq'].ix[best_var],reg_results_df['pval'].ix[best_var],reg_results_df['se'].ix[best_var],
                            best_var+' solar_exposure_day','Tower Solar_obs',OutputPath_DailyData,Site_ID)
        else:
            print 'Climatology provides best performance of available fill variables - skipping regression...'
        
        # Stack the corrected series together (i.e. progressively filling gaps) in order of (minimum) RMSE rank
        DaySum_DF['combined_est']=np.nan
        for i in RMSE_df.index:
            DaySum_DF['combined_est']=np.where(np.isnan(DaySum_DF['combined_est']),DaySum_DF[i],DaySum_DF['combined_est'])
            
        # Fill all gaps in Obs
        DaySum_DF['Obs']=np.where(np.isnan(DaySum_DF['Obs']),DaySum_DF['combined_est'],DaySum_DF['Obs'])
    
    # If no observational data
    else:
        
        # Stack the available series together (i.e. progressively filling gaps)
        DaySum_DF['Obs']=np.nan
        for i in var_list:
            DaySum_DF['Obs']=np.where(np.isnan(DaySum_DF['Obs']),DaySum_DF[i],DaySum_DF['Obs'])
    

    ###### Find extinction coefficient via non-linear optimisation ######    
    
    if Fsd_avail:
                                                                                                    
        # Find maximum insolation and corresponding day of year for each week
        groupbyObj=DaySum_DF.dropna(how='any',axis=0).groupby(lambda x: x.week)
        Max_DF=(pd.DataFrame({'DOY':groupbyObj['Obs'].idxmax().
    	       apply(lambda x: int64(dt.datetime.strftime(x,'%j'))),'Max':groupbyObj['Obs'].max()}))
        Max_DF=Max_DF[optimisation_subset_list[0]:optimisation_subset_list[1]] #This is arbitrary subset that removes wet season values from optimisation!
        
        # Optimise extinction coefficient
        k_opt,k_cov=curve_fit(Insol_calc,np.array(Max_DF['DOY']),np.array(Max_DF['Max']),p0=0.1)
        print 'Site extinction coefficient: ' + str(round(k_opt[0],3))
    
    else:
        
        # Use default extinction coefficient
        k_opt=default_k
    

    ###### Calculate final fluxes ######
    
    #Switch the output flag to true for final data output
    boolOutput=True   
                      
    #Calculate clear sky Kdown for all days and times in dataset
    array_Kdown_clr_daily,array_Kdown_clr_hr=Insol_calc(np.array(DaySum_DF['DOY']),k_opt)
    
    #Calculate and apply the scaling factor for each day and find cloud factor (clf)
    array_Kdown_cld_hhr,array_CLF=daily_scale(np.array(DaySum_DF['Obs']),array_Kdown_clr_daily,array_Kdown_clr_hr)
    
    #Create df to take results of calculations
    rslt_DF=pd.DataFrame({'Kdown':array_Kdown_cld_hhr},index=time_out(DaySum_DF.index[0],len(DaySum_DF)*1440/rec_length,fluxfreq,timestamp_format))
    
    #Reindex to ensure data has same datetimeindex as FLUXdataframe
    rslt_DF=rslt_DF.ix[FLUXDataframe.index[0]:FLUXDataframe.index[-1]]
    
    #Insert daily clear sky and cloud factor estimates into daily DF
    DaySum_DF['Solar_clear_sky']=array_Kdown_clr_daily
    DaySum_DF['Cloud_factor']=array_CLF


    ###### Generate output and return new variables ######

    if write_to_DF==True:
	
	#Create a new label for pandas df column for the constructed variable (and the QC flag) and column to fill label
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
	
	FLUXDataframe[construct_label][FLUXDataframe[construct_flag_label]==99]=rslt_DF['Kdown']
	FLUXDataframe[corr_label]=rslt_DF['Kdown']
 
         
    #Output daily data
    DaySum_DF.to_csv(OutputPath_DailyData+'/Daily Solar Calculations AWAP and Obs for '+Site_ID+'.csv', sep=',',index_label='DT')
    
    
    print "Finished Fsd gapfilling returning DF"
    return FLUXDataframe
