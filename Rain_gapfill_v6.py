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
#import ephem
import math
from pylab import *
from scipy import stats, optimize
import os
import pickle
import operator
import matplotlib.dates as mdates
import matplotlib.cbook as cbook


#Import custom code for processing
import Gap_Fill_climatology_v2 as gap_fill_climatology
import Timezone_v3a as TimeZone

def mean_numbers(frame):
    return (frame.mean())

def Doplots_weekly(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,sorted_Rainstats):   

    #for ID, slope,val_r2 in sorted_Rainstats:
    
    #Do Montly Plots
    print "Doing Monthly  plot"
    figure(1)
    ax1=subplot(111)
    ax1title='Monthly rainfall all sources ALL years for '+Site_ID
    ax1.set_title(ax1title)
    
    xdata=[]
    xlabel=[]
    for index,ID in enumerate(bestAWS_ID):  
	if index==3:
	    corr_label='Rainf_CABLE_mm'
	else:
	    corr_label=str(variable_to_fill+"_"+str(ID)+"_new")
	xdata.append(Rainfall_DF[corr_label].groupby([lambda x: x.year, lambda x: x.month]).mean())
	xlabel.append(ID)
    
    xdata.append( Rainfall_DF[construct_label].groupby([lambda x: x.year, lambda x: x.month]).mean())
    xlabel.append('Rainfall Construct')
  
    xdata.append( Rainfall_DF[variable_to_fill].groupby([lambda x: x.year, lambda x: x.month]).mean())
    xlabel.append('Rainfall Tower')
    
    rcParams['legend.loc'] = 'best'
    plot(xdata[0],'y',label=str(xlabel[0]))
    #plot(xdata[1],'b--',label=str(xlabel[1]))
    #plot(xdata[2],'r--',label=str(xlabel[2]))
    plot(xdata[3],'g--',label=str(xlabel[3]))
    plot(xdata[4],'m',label=str(xlabel[4]))
    if cable_rain==True:
	plot(xdata[5],linewidth=2, color='b',label=str(xlabel[5]))
    
    ax1.set_ylabel('Rainfall (mm)')  
    legend()
   
    savefig(mypathforResults+'/'+'Monthly rainfall all sources ALL years for '+Site_ID)
    #show()   
    close()

    
def Doplots_at_freq(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,sorted_Rainstats,plot_freq,versionID):   
   
    print "Doing Rainfall plots for "+Site_ID+" at freq "+plot_freq
    
    #Create new data and variables
    by = lambda x: lambda y: getattr(y, x)
    #tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    if plot_freq=='dayofyear': scale_factor= 1 * 48
    if plot_freq=='week': scale_factor=1* 48 * 7
    if plot_freq=='month': scale_factor=1* 48 * 30
    if plot_freq=='year': scale_factor=1 * 48 * 365
    
    #Sort the list based on best rsquared and then fill rain in that priority.
    sorted_Rainstats2 = sorted(sorted_Rainstats, key=operator.itemgetter(0),reverse=False)    
    
    xdata1a=Rainfall_DF['Rainf_CABLE_mm'].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x1label='Rainf_CABLE_mm' + " r2=" +  str("{0:.2f}".format(sorted_Rainstats2[3][1]))
	
    corr_label2=str(variable_to_fill+"_"+sorted_Rainstats2[0][0]+"_new")
    xdata2a=Rainfall_DF[corr_label2].groupby([by('year'),by(plot_freq)]).mean()	 * scale_factor
    x2label="AWS ID "+str(sorted_Rainstats2[0][0]) + " r2=" +  str("{0:.2f}".format(sorted_Rainstats2[0][1]))
    
    corr_label3=str(variable_to_fill+"_"+sorted_Rainstats2[1][0]+"_new")
    xdata3a=Rainfall_DF[corr_label3].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x3label="AWS ID "+str(sorted_Rainstats2[1][0])+ " r2=" +  str("{0:.2f}".format(sorted_Rainstats2[1][1]))
    
    corr_label4=str(variable_to_fill+"_"+sorted_Rainstats2[2][0]+"_new")
    xdata4a=Rainfall_DF[corr_label4].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x4label="AWS ID "+str(sorted_Rainstats2[2][0])+ " r2=" +  str("{0:.2f}".format(sorted_Rainstats2[2][1]))

    xdata5a=Rainfall_DF[construct_label].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x5label='Rainfall Construct'
    
    xdata6a=Rainfall_DF[variable_to_fill].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x6label='Rainfall Tower'
    
    startplot=dt.datetime(Rainfall_DF.index[0].year,Rainfall_DF.index[0].month,1)
    enddate1=Rainfall_DF.index[len(Rainfall_DF)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(Rainfall_DF)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
    
    plot(ydata,xdata1a, color='#848484',linewidth=1.5, label=x1label)
    plot(ydata,xdata2a, color='#33CCFF',linewidth=1.5, label=x2label)	
    plot(ydata,xdata3a, color='#9966FF',linewidth=1.5, label=x3label)	
    plot(ydata,xdata4a, color='#00FFCC',linewidth=1.5, label=x4label)	
    plot(ydata,xdata5a, color='#669900',linewidth=2.5, label=x5label)
    plot(ydata,xdata6a, color='#CC9900',linewidth=2.5, label=x6label)
 

    legend(loc='upper right')
    ylabel('Rainfall sum for freq (mm per freq)', size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    suptitle('Monthly rainfall all sources ALL years for '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID,size=24)
    savefig(mypathforResults+'/Monthly rainfall all sources ALL years for '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
    close()

def Doplots_at_freq_nostats(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,plot_freq,versionID):   
   
    print "Doing Rainfall plots for "+Site_ID+" at freq "+plot_freq
    
    #Create new data and variables
    by = lambda x: lambda y: getattr(y, x)
    #tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    if plot_freq=='dayofyear': scale_factor= 1 * 48
    if plot_freq=='week': scale_factor=1* 48 * 7
    if plot_freq=='month': scale_factor=1* 48 * 30
    if plot_freq=='year': scale_factor=1 * 48 * 365
    
    #Sort the list based on best rsquared and then fill rain in that priority.
    #sorted_Rainstats2 = sorted(sorted_Rainstats, key=operator.itemgetter(0),reverse=False)    
    
    xdata1a=Rainfall_DF['Rainf_CABLE_mm'].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x1label='Rainf_CABLE_mm' 
	
    corr_label2=str(variable_to_fill+"_"+bestAWS_ID[0]+"_new")
    xdata2a=Rainfall_DF[corr_label2].groupby([by('year'),by(plot_freq)]).mean()	 * scale_factor
    x2label="AWS ID "+str(bestAWS_ID[0])
    
    corr_label3=str(variable_to_fill+"_"+bestAWS_ID[1]+"_new")
    xdata3a=Rainfall_DF[corr_label3].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x3label="AWS ID "+str(bestAWS_ID[1])
        
    corr_label4=str(variable_to_fill+"_"+bestAWS_ID[2]+"_new")
    xdata4a=Rainfall_DF[corr_label4].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x4label="AWS ID "+str(bestAWS_ID[2])

    xdata5a=Rainfall_DF[construct_label].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x5label='Rainfall Construct'
    
    xdata6a=Rainfall_DF[variable_to_fill].groupby([by('year'),by(plot_freq)]).mean() * scale_factor
    x6label='Rainfall Tower'
    
    startplot=dt.datetime(Rainfall_DF.index[0].year,Rainfall_DF.index[0].month,1)
    enddate1=Rainfall_DF.index[len(Rainfall_DF)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(Rainfall_DF)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
    
    plot(ydata,xdata1a, color='#848484',linewidth=1.5, label=x1label)
    plot(ydata,xdata2a, color='#33CCFF',linewidth=1.5, label=x2label)	
    plot(ydata,xdata3a, color='#9966FF',linewidth=1.5, label=x3label)	
    plot(ydata,xdata4a, color='#00FFCC',linewidth=1.5, label=x4label)	
    plot(ydata,xdata5a, color='#669900',linewidth=2.5, label=x5label)
    plot(ydata,xdata6a, color='#CC9900',linewidth=2.5, label=x6label)
 

    legend(loc='upper right')
    ylabel('Rainfall sum for freq (mm per freq)', size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    suptitle('Monthly rainfall all sources ALL years for '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID,size=24)
    savefig(mypathforResults+'/Monthly rainfall all sources ALL years for '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
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

def regress_func_origin(x,xlabel,ylabel):
    #Do the regression.  Start by subsetting the two columns required.
    #Then drop any NaN case wise
    #reset (get rid of the index) and drop the index rather than keeping it as a column
    #so it can be passed to linregress
    xnow=x[[xlabel,ylabel]]
    xnow=xnow.dropna(how='any')
    x=xnow[xlabel].dropna().reset_index(drop=True)
    y=xnow[ylabel].dropna().reset_index(drop=True)      
    #Fits a linear fit of the form mx+b to the data

    fitfunc = lambda params, x: params[0] * x #+ params[1]    #create fitting function of form mx+b
    errfunc = lambda p, x, y: fitfunc(p, x) - y              #create error function for least squares fit

    init_a = 0.5                            #find initial value for a (gradient)
    init_b = min(y)                          #find initial value for b (y axis intersection)
    init_p = np.array((init_a))          #bundle initial values in initial parameters

    #calculate best fitting parameters (i.e. m and b) using the error function
    p1, success = optimize.leastsq(errfunc, init_p.copy(), args = (x, y))
    f = fitfunc(p1, x)          #create a fit with those parameters
    
    correlation = np.corrcoef(x, y)
    rsqu= (correlation[0,1])**2 
    return (p1, f ,rsqu, x, y)

def linear_reg_plot(slope, rsqu, xlabel1, xdata, ydata,variable_to_fill,units,startdate,enddate,n_datapoints,Site_ID,ID,mypathforResults):
    #Produce the plot 
    line = slope*xdata
    plt.plot(xdata, ydata, 'go',xdata, line, ':b' ,linewidth=2) 
    
    #create text for ID and r2 box
    graphtext1=str('slope  ' + str(slope)+'\n' + 'r squared '+ str(rsqu))  
    #create text for start and end dates
    graphtext2=('Data start date: '+str(startdate)+'\n'
	                        +'End date: '+str(enddate)+'\n'
	                        +'Number records: '+str(n_datapoints))     
    plt.figtext(0.7,0.3,graphtext1, bbox=dict())
    plt.figtext(0.5,0.13,graphtext2, bbox=dict())
    plt.title('Tower vs Rainfall '+xlabel1+' at ' +Site_ID)
    plt.xlabel('Tower ' + '('+units+')')
    plt.ylabel('CABLE ' + '('+units+')')
    plt.legend(shadow=True, fancybox=True,loc='best')
    
    savefig(mypathforResults+'/'+'Tower vs AWS site '+ID+' Rainfall at ' +Site_ID)
    #show()   
    close()
    
def adjust9am_rainfall(x):
    #Convert cummulative to 30 minute for both normal and daylight savings time
    currentdate=dt.date(x.index[0].year,x.index[0].month,x.index[0].day)
    #Get timezone info, call routines from external code
    timezone,InDstNow=TimeZone.get_timezone_info(latitude,longitude,currentdate)
    #print "Current date processing is "+str(currentdate)+". It is a daylight saving day: "+str(InDstNow)    
    
    x.ix[0,AWSlabelnew[0]]=float(x.ix[0,AWSlabel[0]])
    x.ix[0,AWSlabelnew[1]]=float(x.ix[0,AWSlabel[1]])
    x.ix[0,AWSlabelnew[2]]=float(x.ix[0,AWSlabel[2]])
    
    x.ix[0,AWSlabelnew_DST[0]]=float(x.ix[0,AWSlabel_DST[0]])
    x.ix[0,AWSlabelnew_DST[1]]=float(x.ix[0,AWSlabel_DST[1]])
    x.ix[0,AWSlabelnew_DST[2]]=float(x.ix[0,AWSlabel_DST[2]])
    
    #If normal
    if InDstNow == '"false"':
	for index1 in range(1,len(x)):
	    x.ix[index1,AWSlabelnew[0]]=float(x.ix[index1,AWSlabel[0]])-float(x.ix[(index1-1),AWSlabel[0]])
	    x.ix[index1,AWSlabelnew[1]]=float(x.ix[index1,AWSlabel[1]])-float(x.ix[(index1-1),AWSlabel[1]])
	    x.ix[index1,AWSlabelnew[2]]=float(x.ix[index1,AWSlabel[2]])-float(x.ix[(index1-1),AWSlabel[2]])
    
    #If daylight savings true
    if InDstNow == '"true"':
	for index1 in range(1,len(x)):
	    x.ix[index1,AWSlabelnew[0]]=float(x.ix[index1,AWSlabel_DST[0]])-float(x.ix[(index1-1),AWSlabel_DST[0]])
	    x.ix[index1,AWSlabelnew[1]]=float(x.ix[index1,AWSlabel_DST[1]])-float(x.ix[(index1-1),AWSlabel_DST[1]])
	    x.ix[index1,AWSlabelnew[2]]=float(x.ix[index1,AWSlabel_DST[2]])-float(x.ix[(index1-1),AWSlabel_DST[2]])   
    #Return a series _new that has the converted time series
    return x

####################
##START MAIN CODE
####################
def Rain_gapfill(variable_to_fill,CABLEandTOWER_DF,myBaseforResults,Site_ID,bestAWS_ID,cable_rain,FluxFreq,latitude1,longitude1,versionID):
    global latitude
    global longitude
    latitude=latitude1
    longitude=longitude1
    
    print "================================================================"
    print "DINGO: Starting Rainfall gap filling"
    print "================================================================"
    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/AWS"):
	os.mkdir(myBaseforResults+"/AWS")
    if not os.path.isdir(myBaseforResults+"/AWS/"+variable_to_fill):
	os.mkdir(myBaseforResults+"/AWS/"+variable_to_fill)
    mypathforResults=myBaseforResults+"/AWS/"+variable_to_fill
    
    #Load AWS dataframe. 
    AWS_combined=pd.read_pickle(myBaseforResults+'/AWS/AWS_combined_'+Site_ID+'.df') 
    #create and populate two lists that will be used in a new dataframe
    #Three labels each for the 3 AWS ID that come from the BEST AWS correlations
    
    
    #Try deleting the columns first.  This is necessary as if the module is run out of sequence this will throw an error as the column is already there
    try:
	del CABLEandTOWER_DF[construct_label]
    except:
	pass

    try:
	del CABLEandTOWER_DF[construct_flag_label]
    except:
	pass

    try:
	del CABLEandTOWER_DF[Newcorr_label]
    except:
	pass   
    
    #Include CABLE/AWAP data.  The Cable data has AWAP as input so they are the same
    #Convert to mm per 30 minutes from mm.s-1.  Then add the lable to the best ID list so it is pocessed later
    #Include cable yes or no?
    if cable_rain==True:
	CABLEandTOWER_DF['Rainf_CABLE_mm']=CABLEandTOWER_DF['Rainf_CABLE']*60*30
	bestAWS_ID.append('Rainf_CABLE_mm')
    else:
	#trick the following code into thinking there is Cable data. Assign zeros.  This will produce a poor correlation and therfeore will not be selected for gap filling
	CABLEandTOWER_DF['Rainf_CABLE_mm']=0
	bestAWS_ID.append('Rainf_CABLE_mm')
    
    #Here create two sets of columns.  One with normal shift and one for Daylight Saving Time shift _DST
    global AWSlabel
    global AWSlabelnew
    global AWSlabel_DST
    global AWSlabelnew_DST
    AWSlabel=[]
    AWSlabelnew=[]
    AWSlabel_DST=[]
    AWSlabelnew_DST=[]
    
    for index,ID in enumerate(bestAWS_ID):
	AWSlabel.append('Rain_'+ID)
	AWSlabelnew.append(variable_to_fill+'_'+ID+'_new')
	AWSlabel_DST.append('Rain_'+ID+'_DST')
	AWSlabelnew_DST.append(variable_to_fill+'_'+ID+'_new_DST')
	
    #Calculate some things for plots
    n_datapoints=len(CABLEandTOWER_DF)
    startdate= CABLEandTOWER_DF.index[0]
    enddate= CABLEandTOWER_DF.index[n_datapoints-1]
	
    #Combine rainfall fields together
    #Here we need to convert rainfall BoM as cummulative rainfall to 9am INTO 30 minute intervals
    #Lets create a series to manipulate then shift itx intervals so that brings it to rain till midnight
    #Then we can group by day to cut by day and apply a function
    #Now already shifted by 30 minutes when we read in the data due to data logger and BoM AWS timestamp offset
    
    #Create 3 new columns that will be used for the daylight savings shofted column 
    AWS_combined[AWSlabel_DST[0]]=AWS_combined[AWSlabel[0]]
    AWS_combined[AWSlabel_DST[1]]=AWS_combined[AWSlabel[1]]
    AWS_combined[AWSlabel_DST[2]]=AWS_combined[AWSlabel[2]]

    #The BoM AWS files typically have much earlier years than required for the Tower gap filling.
    #Sety this to True if you want to process all the AWS years.  It takes much longer
    #Otherwise get start and end dates from the Tower file and select data based on that
    process_all = False
    if process_all== True:
	rainfall_AWS=AWS_combined[[AWSlabel[0],AWSlabel[1],AWSlabel[2],AWSlabel_DST[0],AWSlabel_DST[1],AWSlabel_DST[2]  ]]
    else:
	rainfall_AWS=AWS_combined[[AWSlabel[0],AWSlabel[1],AWSlabel[2],AWSlabel_DST[0],AWSlabel_DST[1],AWSlabel_DST[2]  ]][startdate:enddate]        

    #Make sure that file is at right frequency
    #Check dataframe for duplicates, pad as necessary and sort
    rainfall_AWS.sort(inplace=True)
    rainfall_AWS["index"] = rainfall_AWS.index
    rainfall_AWS.drop_duplicates('index', take_last=True, inplace=True)
    del rainfall_AWS["index"]	
    rainfall_AWS=rainfall_AWS.asfreq(FluxFreq, method=None)      
    
    #Shift rainfall by 1 to get it back to an hourly timestamp becayuse it was previously shifted
    #due to datalogger timestamp mismatch
    rainfall_AWS=rainfall_AWS.shift(-1)
    # Convert all blank cells to nans.  Get rid of the spaces in the file
    rainfall_AWS = rainfall_AWS.applymap(lambda x: np.nan if isinstance(x, basestring) and x.isspace() else x)
    #Then fill forward for max gap of 1.  This ensures that we only fill forward from the hourly data 
    #which has every second line a gap, to half hourly
    rainfall_AWS=rainfall_AWS.fillna(method='ffill', limit=1)
    
    #Shift from 9am to midnight.
    #Rainfall shifts with daylight savings time
    #Do normal shift for normal time
    rainfall_AWS_18=rainfall_AWS[[AWSlabel[0],AWSlabel[1],AWSlabel[2]]].shift(-18)
    #Do DST shift for daylight savings time
    rainfall_AWS_16=rainfall_AWS[[AWSlabel_DST[0],AWSlabel_DST[1],AWSlabel_DST[2]]].shift(-16)  
    #Recombine columns
    rainfall_AWS=rainfall_AWS_18.join(rainfall_AWS_16,how='left')   
    
    rainfall_AWS[AWSlabelnew[0]]=np.NaN
    rainfall_AWS[AWSlabelnew[1]]=np.NaN
    rainfall_AWS[AWSlabelnew[2]]=np.NaN   

    rainfall_AWS[AWSlabelnew_DST[0]]=np.NaN
    rainfall_AWS[AWSlabelnew_DST[1]]=np.NaN
    rainfall_AWS[AWSlabelnew_DST[2]]=np.NaN      
    
    
    #Apply function to convert to 
    #Will display the day of year being processed as it  goes
    print "Apply function to convert to intervals"
    daily_DF=rainfall_AWS.groupby([lambda x: x.year, lambda x: x.month,lambda x: x.day]).apply(adjust9am_rainfall)
    
    #There may still be some times when the rainfall doesnt align and totals are negative.   Force these to 0
    print "Number of negative values in timeseries for deletion "+AWSlabelnew[0]+" is "+str(daily_DF[AWSlabelnew[0]][daily_DF[AWSlabelnew[0]]<0].count())
    print "Number of negative values in timeseries for deletion "+AWSlabelnew[1]+" is "+str(daily_DF[AWSlabelnew[1]][daily_DF[AWSlabelnew[1]]<0].count())
    print "Number of negative values in timeseries for deletion "+AWSlabelnew[2]+" is "+str(daily_DF[AWSlabelnew[2]][daily_DF[AWSlabelnew[2]]<0].count())
    
    print "Number of negative values in timeseries for deletion "+AWSlabelnew_DST[0]+" is "+str(daily_DF[AWSlabelnew_DST[0]][daily_DF[AWSlabelnew_DST[0]]<0].count())
    print "Number of negative values in timeseries for deletion "+AWSlabelnew_DST[1]+" is "+str(daily_DF[AWSlabelnew_DST[1]][daily_DF[AWSlabelnew_DST[1]]<0].count())
    print "Number of negative values in timeseries for deletion "+AWSlabelnew_DST[2]+" is "+str(daily_DF[AWSlabelnew_DST[2]][daily_DF[AWSlabelnew_DST[2]]<0].count())
    
    daily_DF[AWSlabelnew[0]][daily_DF[AWSlabelnew[0]]<0]=0
    daily_DF[AWSlabelnew[1]][daily_DF[AWSlabelnew[1]]<0]=0
    daily_DF[AWSlabelnew[2]][daily_DF[AWSlabelnew[2]]<0]=0
    
    daily_DF[AWSlabelnew_DST[0]][daily_DF[AWSlabelnew_DST[0]]<0]=0
    daily_DF[AWSlabelnew_DST[1]][daily_DF[AWSlabelnew_DST[1]]<0]=0
    daily_DF[AWSlabelnew_DST[2]][daily_DF[AWSlabelnew_DST[2]]<0]=0
    
    #now join original Rainfall AWS and relevant rainnfall from Fluxtower into a Pandas D
    Rainfall_DF=CABLEandTOWER_DF[[variable_to_fill,'Rainf_CABLE_mm']].join(daily_DF,how='inner') 
    
    #########################################################################################################
    #Check to see if at least 33% of the Rainfall data is present to make this work.  Otherwise use AWAP data.                         
    #########################################################################################################
    Var_total=len(CABLEandTOWER_DF[variable_to_fill])
    Var_notnull=CABLEandTOWER_DF[variable_to_fill][CABLEandTOWER_DF[variable_to_fill].notnull()].count()
    print "Number of values in series: ",Var_total
    print "Number of valid values in series: ", Var_notnull
    #Make an exclude list.  If siteID in the list then only do CABLE gapfill
    #excludeSiteID = ['Daintree']
    #and (Site_ID not in excludeSiteID):
    if (Var_notnull>(Var_total/5)):
	print "Number of valid values greater than a third so do AWS processing"    

	#Do this for each AWS ID.  Build table of Rain stats to use later 
	Rainstats=[]
	for index,ID in enumerate(bestAWS_ID):
	    if index==3:
		data_to_plot=Rainfall_DF[[variable_to_fill,'Rainf_CABLE_mm']].dropna(how='any')
	    else:
		data_to_plot=Rainfall_DF[[variable_to_fill,AWSlabelnew[index]]].dropna(how='any')
	    
	    data_to_plot_grouped=data_to_plot.groupby([lambda x: x.year, lambda x: x.week]).mean()
	    #set lables to send to linear regression function
	    if index==3:
		ylabel=variable_to_fill
		xlabel='Rainf_CABLE_mm'    
	    else:
		ylabel=variable_to_fill
		xlabel=variable_to_fill+'_'+ID+'_new'
	    #Call regression function on dataframe
	    print "Doing regression for AWS ID ",Site_ID
	    #slope, intercept, r_value, p_value, std_err, xdata, ydata = regress_func(data_to_plot_grouped,xlabel,ylabel)
	    #Tried the linear regression but that gives an intercept which when applied gave a small rainfall offset
	    #which meant that there was always a small rainfall amount on the gap filled rain.
	    #Try using this function whcih has no intercept.
	    p1,f, r2, xdata, ydata=regress_func_origin(data_to_plot_grouped,xlabel,ylabel)
	    p1=float(p1)
	    Rainstats.append([ID,p1,r2])    
	    print 'Correlation of Site AWS ID '+str(Rainstats[index][0]) + ' the slope is '+ str(Rainstats[index][1])
	    print 'Correlation of Site AWS ID '+str(Rainstats[index][0]) + ' the r squared is '+ str(Rainstats[index][2])        
	    
	    ###########################
	    #Plot linear regression.
	    ##########################
	    units='mm'  
	    #Calculate some things for plots
	    n_datapoints_plot=len(data_to_plot_grouped)
	    startdate_plot= Rainfall_DF.index[0]
	    enddate_plot= Rainfall_DF.index[n_datapoints_plot-1]  
	    #enddate_plot= Rainfall_DF.index[n_datapoints-1]  
	    
	    #Call linear regression plot    
	    linear_reg_plot(Rainstats[index][1], Rainstats[index][2], xlabel, xdata, ydata,variable_to_fill,units,startdate_plot,enddate_plot,n_datapoints_plot,Site_ID,ID,mypathforResults)
	    
	    #Create a series for each AWS ID that is just the correlated output of CABLE adjusted to Tower
	    corr_label=str(variable_to_fill+"_Corr"+str(ID))
	    Rainfall_DF[corr_label]=Rainfall_DF[xlabel]*Rainstats[index][1]
    
	############################################################
	# Now gapfill Actual Precip series
	# APPLY GAP FILL SWC AND TS   
	# Create variable XX_con for constrcuted variable
	#
	############################################################
	
	#Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
	construct_label=str(variable_to_fill+"_Con")
	Rainfall_DF[construct_label]=np.nan
	#add a column for the constructed QC flag
	#This will be 1 if valid data from the tower else 99
	construct_flag_label=str(variable_to_fill+"_Con_QCFlag")
	Rainfall_DF[construct_flag_label]=np.nan
	#Start by filling with valid rainfal values from tower
	Rainfall_DF[construct_flag_label][Rainfall_DF[variable_to_fill].notnull()]=1
	Rainfall_DF[construct_label][Rainfall_DF[variable_to_fill].notnull()]=Rainfall_DF[variable_to_fill]
	
	#Sort the list based on best rsquared and then fill rain in that priority.
	sorted_Rainstats = sorted(Rainstats, key=operator.itemgetter(2),reverse=True)
	print "Some stats for Rainfall.  All period means for all sources"
	Rainstats[index][1], Rainstats[index][2]
	
	#Also create a SINGLE new label for pandas df column for the contructed variable based on ALL 3 AWS sites and CABLE/AWAP
	#i.e XXX_Corr, No tower data
	Newcorr_label=str(variable_to_fill+"_Corr")
	Rainfall_DF[Newcorr_label]=np.nan
	
	for ID, slope,val_r2 in sorted_Rainstats:
	    #Gap fill  based on priority of best r2 to worst including 3 AWS sites and CABLE/AWAP
	    #Use the previously 'corrected data' i.e. AWS adjusted to Tower rainfall using rehression and slope i.e. XX_Corr
	    corr_label=str(variable_to_fill+"_Corr"+str(ID))
	    if ID=='Rainf_CABLE_mm':
		Rainfall_DF[construct_flag_label][Rainfall_DF[construct_flag_label].isnull()]=100
		Rainfall_DF[construct_label][Rainfall_DF[construct_flag_label]==100]=Rainfall_DF['Rainf_CABLE_mm']
		Rainfall_DF[Newcorr_label][Rainfall_DF[Newcorr_label].isnull()]=Rainfall_DF['Rainf_CABLE_mm']	
	    else:
		Rainfall_DF[construct_flag_label][Rainfall_DF[construct_flag_label].isnull()]=99
		Rainfall_DF[construct_label][Rainfall_DF[construct_flag_label]==99]=Rainfall_DF[corr_label]
		Rainfall_DF[Newcorr_label][Rainfall_DF[Newcorr_label].isnull()]=Rainfall_DF[corr_label]
	    print "For source ", ID, " mean rainfall ", Rainfall_DF[corr_label].mean()
	print "For source TOWER mean rainfall ", Rainfall_DF[variable_to_fill].mean()
	
	#Do plot of  rainfall at frequency
	plot_freq="month"
	Doplots_at_freq(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,sorted_Rainstats,plot_freq,versionID)

	plot_freq="week"
	Doplots_at_freq(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,sorted_Rainstats,plot_freq,versionID)

    

    else:
	##########################################################################################
	# Number of valid values NOT greater than a third so just using CABLE data for rainfall gapfilling
	##############################################################################################
	print "Number of valid values NOT greater than a third so just using CABLE or AWS data for rainfall gapfilling"
	
	#Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
	construct_label=str(variable_to_fill+"_Con")
	Rainfall_DF[construct_label]=np.nan
	#add a column for the constructed QC flag
	#This will be 1 if valid data from the tower else 99
	construct_flag_label=str(variable_to_fill+"_Con_QCFlag")
	Rainfall_DF[construct_flag_label]=np.nan
	#Start by filling with valid rainfal values from tower
	Rainfall_DF[construct_flag_label][Rainfall_DF[variable_to_fill].notnull()]=1
	Rainfall_DF[construct_label][Rainfall_DF[variable_to_fill].notnull()]=Rainfall_DF[variable_to_fill]

	#Also create a SINGLE new label for pandas df column for the contructed variable based on ALL 3 AWS sites and CABLE/AWAP
	#i.e XXX_Corr, No tower data
	Newcorr_label=str(variable_to_fill+"_Corr")
	Rainfall_DF[Newcorr_label]=np.nan

	#Gap fill  using best AWS - This will not be correlated with site based data because there is not enough 
	#to do this robustly
	ID= bestAWS_ID[0]
	print "Using AWS ID number ", str(ID)," for gap filling"
	
	original_label='Rain_'+str(ID)
	New_label='Precip_'+str(ID)+'_new'
	
	#Here if Cable data is selected then use that.  Else use nearest AWS (UNCORRELATED)
	if cable_rain==True:
	    Rainfall_DF[construct_flag_label][Rainfall_DF[construct_flag_label].isnull()]=100
	    Rainfall_DF[construct_label][Rainfall_DF[construct_flag_label]==100]=Rainfall_DF['Rainf_CABLE_mm']
	    Rainfall_DF[Newcorr_label][Rainfall_DF[Newcorr_label].isnull()]=Rainfall_DF['Rainf_CABLE_mm']	
	else:
	    Rainfall_DF[original_label]=Rainfall_DF[original_label].astype('float64')

	    Rainfall_DF[construct_flag_label][Rainfall_DF[construct_flag_label].isnull()]=93
	    Rainfall_DF[construct_label][Rainfall_DF[construct_flag_label]==93]=Rainfall_DF[New_label]
	    Rainfall_DF[Newcorr_label][Rainfall_DF[Newcorr_label].isnull()]=Rainfall_DF[New_label]
	
	#Do plot of  rainfall at frequency
	plot_freq="month"
	Doplots_at_freq_nostats(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,plot_freq,versionID)

	plot_freq="week"
	Doplots_at_freq_nostats(Rainfall_DF,Site_ID,bestAWS_ID,construct_label,variable_to_fill,cable_rain,mypathforResults,plot_freq,versionID)
	

    #Write out Rainfall file
    Rainfall_DF.to_csv(mypathforResults+'/'+'Rainfall_DF_for_' +Site_ID+'.csv')
    Rainfall_DF.to_pickle(mypathforResults+'Rainfall_DF_for_'+Site_ID+'.df')  
    
    #Delete any existing columns in DF
    try:
	CABLEandTOWER_DF=CABLEandTOWER_DF.drop([construct_label,construct_flag_label,Newcorr_label],axis=1)
    except:
	pass
    #add new columns to the Flux Tower DF
    #now join original Rainfall AWS and relevant rainnfall from Fluxtower into a Pandas DF
    CABLEandTOWER_DF=CABLEandTOWER_DF.join(Rainfall_DF[[construct_label,construct_flag_label,Newcorr_label]],how='left') 
    
    ###########################################
    # Call function to do climatology gap fill#
    ###########################################
    CABLEandTOWER_DF=gap_fill_climatology.climatology_monthly_diurnal(CABLEandTOWER_DF,variable_to_fill)    
       
    #####################
    # Finish up
    ######################
    print "================================================================"
    print "DINGO: Starting Rainfall gap filling "+ variable_to_fill
    print "================================================================"
    
    return CABLEandTOWER_DF




