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
from scipy import stats
import os

#Import custom code modules required for adv processing
import Gap_Fill_climatology_v2 as gap_fill_climatology
 

def mean_numbers(frame):
    return (frame.mean())

def Doplots_diurnal(PlottingDF,PlotMonth):
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
    savefig(mypathforResults+'/CABLE ensemble diurnal average for month = '+str(PlotMonth)+' '+Site_ID)
    close()
    #show()    
 
def Doplots_monthly(CABLEandTOWER_DF):   
    #Do Montly Plots
    print "Doing Weekly  plot"
    t = arange(1, 54, 1)

    xdata1a = CABLEandTOWER_DF['SWdown_CABLE'].groupby([lambda x: x.week]).mean()
    xdata1b = CABLEandTOWER_DF['LWdown_CABLE'].groupby([lambda x: x.week]).mean()
    xdata1c = CABLEandTOWER_DF['Fsd'].groupby([lambda x: x.week]).mean()
    xdata1d = CABLEandTOWER_DF['Fld'].groupby([lambda x: x.week]).mean()
    
    xdata2a = CABLEandTOWER_DF['Rnet_CABLE'].groupby([lambda x: x.week]).mean()
    xdata2b = CABLEandTOWER_DF['Qh_CABLE'].groupby([lambda x: x.week]).mean()
    xdata2c = CABLEandTOWER_DF['Qg_CABLE'].groupby([lambda x: x.week]).mean()
    xdata2d = CABLEandTOWER_DF['Qle_CABLE'].groupby([lambda x: x.week]).mean()
    xdata2e = CABLEandTOWER_DF['Qs_CABLE'].groupby([lambda x: x.week]).mean()

    xdata3a = CABLEandTOWER_DF['NEE_CABLE'].groupby([lambda x: x.week]).mean()
    xdata3b = (CABLEandTOWER_DF['GPP_CABLE'].groupby([lambda x: x.week]).mean())*-1.0
    xdata3c = CABLEandTOWER_DF['Fc'].groupby([lambda x: x.hour]).mean()
    
    rcParams['legend.loc'] = 'best'
    print shape(t)
    print shape(xdata3c)
    
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
    savefig(mypathforResults+'/CABLE ensemble monthly plot all years for '+Site_ID)
    close()
    #show()    


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

def linear_reg_plot(slope, intercept, r_value, p_value, std_err, xdata, ydata,VarToCorrelate,units,startdate,enddate,n_datapoints,Site_ID):
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
    
    savefig(mypathforResults+'Tower vs CABLE for Variable '+VarToCorrelate+ ' at ' +Site_ID)
    close()
    #show()  
    
####################
##START MAIN CODE
####################
def CABLE_gapfill(variable_to_fill,CABLEandTOWER_DF,myBaseforResults,Site_ID):
    print ' Start gap filling of variable '+ variable_to_fill +' using CABLE LSM data for site '+ Site_ID
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/CABLE"):
	os.mkdir(myBaseforResults+"/CABLE")
    global mypathforResults
    mypathforResults=myBaseforResults+"/CABLE/"
    
    #################################################
    ##Do comparison with Fluxtower data is required #
    ################################################
    #First check to see if there is any data.  The Cable files stop at 2010.  If there is no data skip this part
    if len(CABLEandTOWER_DF[['Ts_CABLE','Ts']].dropna(how='any'))>1000:
	CABLEData=True
    else:
	CABLEData=False
    
    ylabel=variable_to_fill
    xlabel=variable_to_fill+'_CABLE'
	
    if CABLEData==True:
	#Do corelation with Variable 
	#create dataframe of daily average
	CABLEandTOWER_daily_DF=CABLEandTOWER_DF.groupby([lambda x: x.year, lambda x: x.strftime('%j')]).mean()

	
	#Call regression function on dataframe
	slope, intercept, r_value, p_value, std_err, xdata, ydata = regress_func(CABLEandTOWER_daily_DF,xlabel,ylabel)
	Var_stats=[slope, intercept, r_value, p_value, std_err]
	
	#Calculate some things for plots
	n_datapoints=len(CABLEandTOWER_DF)
	startdate= CABLEandTOWER_DF.index[0]
	enddate= CABLEandTOWER_DF.index[n_datapoints-1]
	
	#Plot linear regression.  Try and detect assign the units based on variable to fill
	if variable_to_fill=='Ts':
	    units='oC'
	elif variable_to_fill=='Sws':
	    units='m3 m-3'
	elif variable_to_fill=='Ta':
	    units='oC'
	elif variable_to_fill=='Fsd':
	    units='W m-2'    
	elif variable_to_fill=='Fld':
	    units='W m-2' 
	elif variable_to_fill=='Ah':
	    units='g m-3'
	elif variable_to_fill=='Ws':
	    yunits='m s-1' 
	elif variable_to_fill=='P':
	    yunits='kPa'
	elif variable_to_fill=='ps':
	    yunits='kPa'      
	elif variable_to_fill=='Precip':
	    yunits='mm'
	elif variable_to_fill=='Rain':
	    yunits='mm'  
	else:
	    units=''
	
       
	#Call linear regression plot    
    
	linear_reg_plot(slope, intercept, r_value, p_value, std_err, xdata, ydata,variable_to_fill,units,startdate,enddate,n_datapoints,Site_ID)
	
    
    ############################################################
    # APPLY GAP FILL SWC AND TS                                #
    ############################################################
    
    #Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
    construct_label=str(variable_to_fill+"_Con")
    CABLEandTOWER_DF[construct_label]=np.nan
    #add a column for the constructed QC flag
    #This will be 1 if valid data from the tower else 99
    construct_flag_label=str(variable_to_fill+"_Con_QCFlag")
    CABLEandTOWER_DF[construct_flag_label]=np.nan
    
    CABLEandTOWER_DF[construct_flag_label][CABLEandTOWER_DF[variable_to_fill].notnull()]=1
    CABLEandTOWER_DF[construct_label][CABLEandTOWER_DF[construct_flag_label]==1]=CABLEandTOWER_DF[variable_to_fill]
    #CABLEandTOWER_DF[construct_label][CABLEandTOWER_DF[variable_to_fill].notnull()&CABLEandTOWER_DF[xlabel].notnull()]=CABLEandTOWER_DF[variable_to_fill]

    
    #Apply the linear regression from the stats
    if CABLEData==True:
	CABLEandTOWER_DF[construct_flag_label][CABLEandTOWER_DF[construct_flag_label].isnull()]=99
	CABLEandTOWER_DF[construct_label][CABLEandTOWER_DF[construct_flag_label]==99]=CABLEandTOWER_DF[xlabel]*Var_stats[0] + Var_stats[1]
	#Create a series that is just the correlated output of CABLE adjusted to Tower
	corr_label=str(variable_to_fill+"_Corr")
	CABLEandTOWER_DF[corr_label]=CABLEandTOWER_DF[xlabel]*Var_stats[0] + Var_stats[1]            

    if CABLEData==False:	
	#Create a dummy Corr series so that the rest of the code doesnt fail.
	corr_label=str(variable_to_fill+"_Corr")
	CABLEandTOWER_DF[corr_label]=np.nan	
	
    ###########################################
    # Call function to do climatology gap fill#
    ###########################################

    CABLEandTOWER_DF=gap_fill_climatology.climatology_monthly_diurnal(CABLEandTOWER_DF,variable_to_fill)
    
    #####################
    # Finish up
    ######################
    
    print "FINISHED Gapfilling using CABLE data for variable ",variable_to_fill
    return CABLEandTOWER_DF




