############################################################################
# This script correlates tower data with external meteorology for a variable.  Then calulates the correlation coefficient.  Then adjusts new time series and gap fills 
# Inputs:
#    base_dir_netcdf : Path where to save the netcdf dataset (e.g. '/.')
#    latlim : Latitude limits (in degrees) of the domain to examine (e.g [-40.0, -35.0])
#    lonlim : Longitude limits (in degrees) of the domain to examine (e.g [140.0, 150.0])
#    excel_file: Path to the Excel file containing the indices to correlate with (e.g 'CIall_mon_new_cropyears.xls')
#    index_name: Name of the column containing the index to correlate with the rainfall data (e.g 'Nino1+2_ANOM')
#    months_offset : Offset (in months) between the climate indices and the rainfall data 
#    output_prefix: Prefix of the output NetCDF file that will be created (e.g. '' will create a file "$(season)_$(index_name).nc)"
#
# Notes:
#   The excel file is suppose to have the dates on the first column.
#   The index name is expected to be found on the first row
#
# Programmed by Jason (Dec 1, 2012)
############################################################################

import pandas as pd
import datetime as dt
import xlrd
import numpy as np
import netCDF4
import time
import urllib2
import string
import re
import xml.etree.ElementTree as ET
import ephem
import math
from pylab import *
from scipy import stats
import os


#Import custom code modules required for adv processing
import Timezone_v3a as TimeZone
import Solar_calcs_v3 as Solar_Calcs

def mean_numbers(frame):
    return (frame.mean())

def Doplots_diurnal(PlottingDF,PlotMonth,Site_ID,mypathforResults):
    print "Doing Weekly  plot for month "+str(PlotMonth)
    #Do Diurnal Plots
    t = arange(1, 25, 1)
    xdata1a = PlottingDF['SWdown_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata1b = PlottingDF['LWdown_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata1c = PlottingDF['Fsd'].groupby([lambda x: x.hour]).mean()
    xdata1d = PlottingDF['Fld'].groupby([lambda x: x.hour]).mean()
    
    plot_solar_noon_mean=PlottingDF['Solar_noon'].mean()
    print  plot_solar_noon_mean
    
    xdata2a = PlottingDF['Rnet_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata2b = PlottingDF['Qh_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata2c = PlottingDF['Qg_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata2d = PlottingDF['Qle_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata2e = PlottingDF['Qs_CABLE'].groupby([lambda x: x.hour]).mean()

    xdata3a = PlottingDF['NEE_CABLE'].groupby([lambda x: x.hour]).mean()
    xdata3b = (PlottingDF['GPP_CABLE'].groupby([lambda x: x.hour]).mean())*-1.0
    xdata3c = PlottingDF['Fc'].groupby([lambda x: x.hour]).mean()

    rcParams['legend.loc'] = 'best'

    figure(1)
    ax1=subplot(311)
    ax1title='CABLE ensemble diurnal average for month = '+str(PlotMonth)+' '+Site_ID
    ax1.set_title(ax1title)
    plot(t,xdata1a,'r',label='SWdown_CABLE')
    plot(t,xdata1b,'b',label='LWdown_CABLE')
    plot(t,xdata1c,'r--',label='Fsd')
    plot(t,xdata1d,'b--',label='Fld')
    axvline(x=12, ymin=0, ymax=1)
    axvline(x=plot_solar_noon_mean*24, ymin=0, ymax=1,c='r')
    ax1.set_ylabel('Radiation (W m^-2)')    
    legend()
    
    ax2=subplot(312,sharex=ax1)
    plot(t,xdata2a,'k',label='Rnet_CABLE')
    plot(t,xdata2b,'r',label='Qh_CABLE')
    plot(t,xdata2c,'g',label='Qg_CABLE')
    plot(t,xdata2d,'b',label='Qle_CABLE')
    plot(t,xdata2e,'y',label='Qs_CABLE')
    axvline(x=12, ymin=0, ymax=1)
    axvline(x=plot_solar_noon_mean*24, ymin=0, ymax=1,c='r')
    ax2.set_ylabel('Energy Balance (W m^-2)')     
    legend()
    
    ax3= subplot(313,sharex=ax1)
    plot(t,xdata3a,'b',label='NEE_CABLE')
    plot(t,xdata3b,'g',label='GPP_CABLE')
    plot(t,(xdata3a-xdata3b),'r',label='Re_CABLE')
    #plot(t,xdata3c,'b--',label='Fc')
    axvline(x=12, ymin=0, ymax=1)
    axvline(x=plot_solar_noon_mean*24, ymin=0, ymax=1,c='r')  
    ax3.set_xlabel('Hour')
    ax3.set_ylabel('Carbon flux (umol m-2 s-1)')     
    legend()
    
    figure(1)
    savefig(mypathforResults+'CABLE ensemble diurnal average for month = '+str(PlotMonth)+' '+Site_ID)
    #show()  
    close()
 
def Doplots_monthly(CABLEandTOWER_DF,Site_ID,mypathforResults):   
    #Do Montly Plots
    print "Doing Weekly  plot"
    
    plot_list=['SWdown_CABLE','LWdown_CABLE','Fsd','Fld','Rnet_CABLE','Qh_CABLE','Qg_CABLE','Qle_CABLE','Qs_CABLE','NEE_CABLE','GPP_CABLE','Fc']
    plot_grouped=CABLEandTOWER_DF[plot_list].groupby([lambda x: x.week]).mean()
    
    xdata1a = plot_grouped['SWdown_CABLE']
    xdata1b = plot_grouped['LWdown_CABLE']
    xdata1c = plot_grouped['Fsd']
    xdata1d = plot_grouped['Fld']
    
    xdata2a = plot_grouped['Rnet_CABLE']
    xdata2b = plot_grouped['Qh_CABLE']
    xdata2c = plot_grouped['Qg_CABLE']
    xdata2d = plot_grouped['Qle_CABLE']
    xdata2e = plot_grouped['Qs_CABLE']

    xdata3a = plot_grouped['NEE_CABLE']
    xdata3b = (plot_grouped['GPP_CABLE'])*-1.0
    xdata3c = plot_grouped['Fc']
    
    rcParams['legend.loc'] = 'best'
    
    ticks=shape(xdata3c)[0]
    t = arange(1, (ticks+1), 1)
    
    figure(2)
    ax1=subplot(311)
    ax1title='CABLE output ensemble average by hour ALL years for '+Site_ID
    ax1.set_title(ax1title)
    plot(t,xdata1a,'r',label='SWdown_CABLE')
    plot(t,xdata1b,'b',label='LWdown_CABLE')
    plot(t,xdata1c,'r--',label='Fsd')
    plot(t,xdata1d,'b--',label='Fld')
    ax1.set_ylabel('Radiation (W m^-2)')  
    legend()
    
    ax2=subplot(312,sharex=ax1)
    plot(t,xdata2a,'k',label='Rnet_CABLE')
    plot(t,xdata2b,'r',label='Qh_CABLE')
    plot(t,xdata2c,'g',label='Qg_CABLE')
    plot(t,xdata2d,'b',label='Qle_CABLE')
    plot(t,xdata2e,'y',label='Qs_CABLE')
    ax2.set_ylabel('Energy Balance (W m^-2)') 
    legend()
    
    ax3= subplot(313,sharex=ax1)
    plot(t,xdata3a,'b',label='NEE_CABLE')
    plot(t,xdata3b,'g',label='GPP_CABLE')
    plot(t,(xdata3a-xdata3b),'r',label='Re')
    #plot(t,xdata3c,'b--',label='Fc')
    ax3.set_xlabel('Week')
    ax3.set_ylabel('Carbon flux (umol m-2 s-1)')         
    legend()
    
    figure(2)
    savefig(mypathforResults+'CABLE ensemble montly plot all years for '+Site_ID)
    #show()  
    close()   


def regress_func(x,xlabel,ylabel):
    #Do the regression.  Start by subsetting the two columns required.
    #Then drop any NaN case wise
    #reset (get rid of the index) and drop the index rather than keeping it as a column
    #so it can be passed to linregress
    xnow=x[[xlabel,ylabel]]
    xnow=xnow.dropna(how='any')
    xdata=xnow[xlabel].dropna().reset_index(drop=True)
    ydata=xnow[ylabel].dropna().reset_index(drop=True)   

   #Check to see if minumum number of samples before writing the results.
    if xnow[xlabel].count()>10 and xnow[ylabel].count()>10:
	slope, inter, rsqu, pval, se= stats.linregress(xdata,ydata)
	print "stats:",slope, inter, rsqu, pval, se
	#Here use the original column to do apply the lin regresssion as 
	#values had been dropped previously
	#x[newcolumn]=slope*x[xlabel]+inter
    #else:
	#x[newcolumn]=np.NaN
     #Return the new column to then join to existing file
    return (slope, inter, rsqu, pval, se, xdata, ydata)

def linear_reg_plot(slope, intercept, r_value, p_value, std_err, xdata, ydata,VarToCorrelate,units):
    #Produce the plot 
    line = slope*xdata+intercept
    plt.plot(xdata, ydata, 'go',xdata, line, ':b' ,linewidth=2) 
    
    #create text for ID and r2 box
    graphtext1=str('intercept  ' + str("{0:.2f}".format(intercept) +'\n')
	                              + 'r value      ' + str("{0:.2f}".format(r_value)) +'\n'+'p_value      ' 
                                      + str("{0:.2f}".format(p_value)) +'\n'
	                              + 'slope       ' + str("{0:.2f}".format(slope)) +'\n'
                                      + 'std_err      ' + str("{0:.2f}".format(std_err)) +'\n')  
    #create text for start and end dates
    graphtext2=('Data start date: '+str(startdate)+'\n'
	                        +'End date: '+str(enddate)+'\n'
	                        +'Number records: '+str(n_datapoints))     
    plt.figtext(0.7,0.3,graphtext1, bbox=dict())
    plt.figtext(0.5,0.13,graphtext2, bbox=dict())
    plt.title('Tower vs CABLE for Variable '+VarToCorrelate+ ' at ' +Site_ID)
    plt.xlabel('Tower ' + '('+units+')')
    plt.ylabel('CABLE ' + '('+units+')')
    plt.legend(shadow=True, fancybox=True,loc='best')
    
    savefig('Tower vs CABLE for Variable '+VarToCorrelate+ ' at ' +Site_ID)
    
    show()
    
####################
##START MAIN CODE
####################

def fetch_CABLE(FLUXDataframe,myBaseforResults,CABLEfilename,Site_ID,latitude,longitude):

    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/CABLE"):
	os.mkdir(myBaseforResults+"/CABLE")
    mypathforResults=myBaseforResults+"/CABLE/"
    
    #Check to see if columns already exist.  If so delete them.  Then process and add them again
    try:
	listnames=  FLUXDataframe.columns
	regex=re.compile(".*(CABLE).*")    
	listtodelete=[m.group(0) for l in listnames for m in [regex.search(l)] if m]   
	for item in listtodelete:
	    del FLUXDataframe[item]
    except:
	pass
    
    try:
	del FLUXDataframe['Solar_noon']
    except:
	pass
    
    #set start date and time for data
    start_date=dt.datetime(1990,1,1)
    #Depths in model used later for labels (in cm)
    #derived from thickness of 0.022, 0.058, 0.07, 0.15, 0.30, 0.30, 0.30, 1.20, 3.0, and 4.5m
    depths=[3,8,15,30,60,90,120,240,540,990]
    
    #Next Fetch the Flux Tower Data
    #start by Opening the NetCDF file
    nc_file = netCDF4.Dataset(CABLEfilename) 
    
    #List all variable names
    nc_variableNames = nc_file.variables.keys() 
    
    #Create some lists of size equal to length of vnames list.
    temp=list(xrange(len(nc_variableNames)))
    vartemp=list(xrange(len(nc_variableNames)))
    
    # call the function to convert to datetime from excel. Assume datemode: 0
    # Need to know the name of variable Exceldatetime_var, defined above
    #times = [excel_to_pydate(elem) for elem in nc_file.variables[Exceldatetime_var]]
    #vartemp[0] = times    
    
    #Create series for each variable then loop through other variables in the list 
    #Start with the time variable
    #temp[0]= nc_file.variables['time']
    #vartemp[0] = pd.Series(temp[0][:],name='time')
    #Now this is in Solar time.  Need to convert
    
    #Get number of lines in the CABLE nc file then use that to create a TIME series based on startdate and 30 minute period
    #Assume that the model has PERFECT time stamps.
    n_datapoints_ncfile=len(nc_file.variables['time'])
    rng=pd.date_range(start_date,periods=n_datapoints_ncfile,freq='30min')

    index_vars=1
	
    currentdate="2013-06-01"
    #Get timezone info, call routines from external code
    timezone,InDstNow=TimeZone.get_timezone_info(latitude,longitude,currentdate)
    print "AskGEO TimZone offset (hrs):  ", timezone
    
  
    
    #Should we shift the CABLE output (reported in solar time, centred around solar noon)?
    
    #Only need to do so if using 30 minute information.  If using the daily soil moisture then no need.
    
    #calcualte the difference between local time (recordede by logger and the CABLE output)
    #Need to calculate solar noon for the given site.  Need to inout lat/long
    #Then the code will fetch TimeZone info (hours offset)
    #Then calls the solar calculator to get solar noon for the particular day of year and year.
    #Solar noon time varies across the year resulting in different offset periods (something zero somtimes an hour or more).
    #Will have to shift a single 30 minute time step when the difference in time is greater than 15  minutes.
    #Also if we shift on a daily basis  then as we switch from a 30  minute offset to 1 hour we will leave some gaps.  Is that ok? linear gap fill
    # Also need to take into account the fact that datalogger time is dataloggers report the end of the period as the time stamp. 
    # What about Cable?
    
    shifttime=True
    
    for variable in (nc_variableNames):  
    
	try:
	    number_of_subarrays = nc_file.variables[variable].shape[2]
	    #For multiple dimension variables like Ts and SWC we need to break this down
	    #Should be 11 levels in SOIL
	    for subvariable_index in range(nc_file.variables[variable].shape[1]):
		
		temps= nc_file.variables[variable]
		tempsubvariable = temps[:,subvariable_index,]
		tempsubvariable2 = tempsubvariable[:,0]            
		tempsubvar_name = variable + '_' + str(depths[subvariable_index])+"cm"+"_CABLE"
		#A[:,2] # returns the xth columm
		vartemp[index_vars] = pd.Series(tempsubvariable2[:],name=tempsubvar_name,index=rng)
		index_vars += 1
	except:
	    #Get all the cable variables that are time series.  quuickly check to see if valid.  we expect
	    #about 370992
	    try:
		a=nc_file.variables[variable].shape[0]
		b=nc_file.variables[variable].shape[1]
		if a>1000 and b==1:
		    temps2 = nc_file.variables[variable]
		    temp[index_vars] = temps2[:,0]   
		    tempname=variable+"_CABLE"
		    vartemp[index_vars] = pd.Series(temp[index_vars][:],name=tempname,index=rng)
		    #print vartemp[index_vars].shape             
		    #print vartemp[index_vars]
		    index_vars += 1   
	    except:
		pass
    
    #Concatenate all the series together
    #Specifically only 1 to 57 elements in the list.  Elements later in the list are just integers
    #Element 0 is time and we dont need thins b/c we has made our own time series.We ignore those
    print "Concatenate all the series together"
    theDataFrame=pd.concat((vartemp[1:56]),axis=1)
    
    print "length of theDataFrame",len(theDataFrame)
    
    ##Write out to CSV
    #theDataFrame.to_csv(mypathforResults+'CABLE_model_input_for_'+Site_ID+'.csv', sep=',')
    
    #Write the solar noon to the DF so we can use it later
    theDataFrame['Solar_noon']=np.nan
    
    def shifttimeDF(group):
	#Get date from data group
	startdate= group.index[0]
	basedate=dt.datetime(1900, 1, 1,0,0,0,0)    
	delta=(startdate-basedate).days
	#Calculate Approximate Solar noon, call routines from external code
	#Call Solar_Calcs.solar_calculations(Date_input_excel,latitude,longitude,timezone)
	solar_sunrise,solar_noon_for_date,solar_sunset  =Solar_Calcs.solar_calculations(delta,latitude,longitude,timezone)
	time_difference_minutes=((solar_noon_for_date-0.5)*24*60)
	intervals=round(time_difference_minutes/30)
	#Now fore each 30  minutes we shift
	shiftit=group.shift(intervals)
	shiftit['Solar_noon']=solar_noon_for_date
	#print solar_noon_for_date,intervals
	return shiftit
    
    #Offset by variable amount if variable is True
    if shifttime==True:
	print "Do time shift of CABLE time to local Tower time"
	dailyCABLE_DF=theDataFrame.groupby([lambda x: x.year, lambda x: x.month,lambda x: x.day]).apply(shifttimeDF)
	#gapfill missing numbers due to time shift.  Use backfill because there was a fwd shift
	dailyCABLE_DF=dailyCABLE_DF.fillna(method='backfill')
	
	#Write out to CSV
	#Uncomment if you need the CABLE output independantly
	#dailyCABLE_DF.to_csv(mypathforResults+'CABLE_model_input_for_'+Site_ID+'_shifted.csv', sep=',')



    
    #now join original Pandas DF to the CABLE DF
    CABLEandTOWER_DF=FLUXDataframe.join(dailyCABLE_DF,how='left')    
    #To compare with the tower we need to average CABLE layers SoilMoist_3cm_CABLE and SoilMoist_8cm_CABLE
    print "length of dailyCABLE_DF",len(dailyCABLE_DF)
    print "length of FLUXDataframe",len(FLUXDataframe)
    print "length of CABLEandTOWER_DF",len(CABLEandTOWER_DF)
    
    CABLEandTOWER_DF['Sws_CABLE']=(CABLEandTOWER_DF['SoilMoist_3cm_CABLE']+CABLEandTOWER_DF['SoilMoist_8cm_CABLE'])/2    
    #To compare with the tower we need to average CABLE layers SoilTemp_3cm_CABLE and SoilTemp_8cm_CABLE
    #and convert to oC from Kelvin -272.15
    CABLEandTOWER_DF['Ts_CABLE']=((CABLEandTOWER_DF['SoilTemp_3cm_CABLE']+CABLEandTOWER_DF['SoilTemp_8cm_CABLE'])/2)-272.15    

    #Write out to CSV
    #CABLEandTOWER_DF.to_csv(mypathforResults+'CABLEandTOWER_DF_for_'+Site_ID+'.csv', sep=',')
    # Save dataframe. 
    #CABLEandTOWER_DF.save(myBaseforResults+'CABLEandTOWER_DF_for_'+Site_ID+'.df')  
    
   
    ###########################    
    #Do additional plots if required
    ##########################
    print 'doing plots for CABLE'
    #Create temp DF that can be passed for plotting
    #Select for Jan
    PlotMonth=1
    PlottingDF=CABLEandTOWER_DF[(CABLEandTOWER_DF.index.month==PlotMonth)]
    Doplots_diurnal(PlottingDF,PlotMonth,Site_ID,mypathforResults)
    #Select for July
    PlotMonth=7
    PlottingDF=CABLEandTOWER_DF[(CABLEandTOWER_DF.index.month==PlotMonth)]
    Doplots_diurnal(PlottingDF,PlotMonth, Site_ID,mypathforResults)
    #Do Monthly plots
    Doplots_monthly(CABLEandTOWER_DF,Site_ID,mypathforResults)  
    
    #####################
    # Finish up
    ######################
    
    print "FINISHED CABLE IMPORT"
    return CABLEandTOWER_DF
    
    
