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
import numpy as np
import os
import datetime as dt
import time
import pylab as pl
import meteorologicalfunctions as metfuncs
from scipy import stats

from pylab import figure, ioff, clf, contourf, ion, draw, show
from ffnet._tests import runtest
from ffnet import ffnet, mlgraph, readdata, tmlgraph, imlgraph
from numpy import array

#Import custom code modules required for adv processing
import Timezone_v3a as TimeZone
import Solar_calcs_v3 as Solar_Calcs



def Doplots_diurnal_monthly(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,item,index_str):
    ANN_label=str(item+"_NN")     #Do Monthly Plot
    print "Doing monthly diurnal plot for month and index ",index_str
    #Do Diurnal Plots for all 12 months
    #create an X axis series for all 24 hours
    t = np.arange(1, 25, 1)
    NN_label='Fc'
    Plottemp = PlottingDF[[NN_label,item]]#[PlottingDF['day_night']!=3]
    #Plottemp = PlottingDF[[NN_label,item]].dropna(how='any')
	    
    figure(1)
    pl.subplot(321)
    pl.title('Diurnal '+item+' month = 1')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==1)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==1)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')    
 
    pl.subplot(322)
    pl.title('Diurnal '+item+' month = 3')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==3)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==3)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')      

    pl.subplot(323)
    pl.title('Diurnal '+item+' month = 5')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==5)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==5)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')      
    
    pl.subplot(324)
    pl.title('Diurnal '+item+' month = 7')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==7)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==7)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.subplot(325)
    pl.title('Diurnal '+item+' month = 9')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==9)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==9)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.subplot(326)
    pl.title('Diurnal '+item+' month = 11')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==11)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==11)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.suptitle('Monthly ANN ensemble diurnal  for  '+item+' at '+Site_ID+ ' index '+index_str)
    pl.subplots_adjust(top=0.85)
    pl.tight_layout()  
    pl.savefig(mypathforResults+'/Monthly ANN ensemble diurnal  for  '+item+' at '+Site_ID+ ' index '+index_str)
    #pl.show() 
    pl.close()
    pl.close(1)
    time.sleep(1)
    
def Doplots_diurnal(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,item,index_str):
    ANN_label=str(item+"_NN")     #Do Monthly Plot
    print "Doing diurnal plot for month and index ",index_str
    #Do Diurnal Plots for all 12 months
    #create an X axis series for all 24 hours
    t = np.arange(1, 25, 1)
    NN_label='Fc'
    Plottemp = PlottingDF[[NN_label,item]]#[PlottingDF['day_night']!=3]
    #Plottemp = PlottingDF[[NN_label,item]].dropna(how='any')
	    
    figure(2)
    pl.title('Diurnal '+item )
    try:
	xdata1a=Plottemp[item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')    
    
    pl.suptitle('ANN ensemble diurnal  for  '+item+' at '+Site_ID+ ' index '+index_str)
    pl.subplots_adjust(top=0.85)
    pl.tight_layout()  
    pl.savefig(mypathforResults+'/ANN ensemble diurnal  for  '+item+' at '+Site_ID+ ' index '+index_str)
    #pl.show() 
    pl.close()
    pl.close(2)
    time.sleep(1)
	
def Doplots_monthly(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,item,index_str):   
    ANN_label=str(item+"_NN")     #Do Monthly Plots
    print "Doing Monthly  plot index", index_str
    #t = arange(1, 54, 1)
    NN_label='Fc'
    Plottemp = PlottingDF[[NN_label,item,'ustar','ustar_used','day_night']].dropna(how='any')
    Plottemp=Plottemp[Plottemp['day_night']!=1][Plottemp['ustar']>Plottemp['ustar_used']]
    
    figure(3)
    pl.title('Nightime ANN v Tower ustar filtered by year-month for '+item+' at '+Site_ID+ ' index '+index_str)

    try:
	xdata1a=Plottemp[item].groupby([lambda x: x.year,lambda x: x.month]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[NN_label].groupby([lambda x: x.year,lambda x: x.month]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[item][Plottemp['ustar']>Plottemp['ustar_used']].groupby([lambda x: x.year,lambda x: x.month]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False     

    if plotxdata1a==True:
	pl.plot(xdata1a.values,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(xdata1b.values,'b',label=NN_label)
    if plotxdata1c==True:
	pl.plot(xdata1c.values,'g',label='Fc ustar filt')
    pl.xlabel('Year - Month')       
    pl.legend()
    pl.savefig(mypathforResults+'/ANN and Tower plots by year and month for variable '+item+' at '+Site_ID+ ' index '+index_str)
    #pl.show()
    pl.close(3)
    time.sleep(1)
	
def regressionANN(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,item,index_str):    
    ANN_label=str(item+"_NN") 
    figure(5)
    graphtext1=str('slope      ' + str("{0:.2f}".format(float(regress[0][0]))) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(float(regress[0][1]))) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(float(regress[0][2]))) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(float(regress[0][3]))) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(float(regress[0][4])))       )  
    pl.figtext(0.7,0.6,graphtext1, bbox=dict())
    pl.plot(observed, predicted, 'o', label='targets vs. outputs')    
    slope = regress[0][0]; intercept = regress[0][1]

    x = np.linspace(min(observed),max(observed))
    y = slope * x + intercept
    pl.plot(x, y, linewidth = 2, label = 'regression line')
    pl.legend()
    pl.title('Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str)
    pl.xlabel('Tower ' + '('+units+')')
    pl.ylabel('ANN ' + '('+units+')')
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.savefig(mypathforResults+'/'+'Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str) 
    #pl.show()
    pl.close(5)    
    time.sleep(1)

def regressionANN2(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,item,index_str):    
    ANN_label=str(item+"_NN") 
    figure(5)
    predicted=np.reshape(predicted,(len(predicted)))
    observed=np.reshape(observed,(len(observed)))
    
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(predicted,observed)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.7,0.6,graphtext1, bbox=dict())
    pl.plot(observed, predicted, 'o', label='targets vs. outputs')    
    

    x = np.linspace(min(observed),max(observed))
    y = slope * x + intercept
    pl.plot(x, y, linewidth = 2, label = 'regression line')
    pl.legend()
    pl.title('Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str)
    pl.xlabel('Tower ' + '('+units+')')
    pl.ylabel('ANN ' + '('+units+')')
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.savefig(mypathforResults+'/'+'Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str) 
    #pl.show()
    pl.close(5)    
    time.sleep(1)

def mintimeseries_plot(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,targets,output,item,index_str):    
    ANN_label=str(item+"_NN")    
    figure(4)
    pl.plot( targets, 'b--' )
    pl.plot( output, 'k-' )
    pl.legend(('targets', 'output'))
    pl.xlabel('Time'); 
    pl.title('Outputs vs. target of trained network for '+item)
    pl.grid(True)
    pl.legend()
    pl.title('Tower vs ANN 30 min timeseries for '+item+' at ' +Site_ID+ ' index '+index_str)
    pl.ylabel(item + '('+units+')')
    pl.savefig(mypathforResults+'/'+'Tower vs ANN 30 min timeseries for '+item+' at ' +Site_ID+ ' index '+index_str) 
    #pl.show()
    pl.close(4)
    time.sleep(1)
    
def Fre_ANN_gapfill_func(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,latitude,longitude,index_str,ANN_label,frequency,evening,min_threshold,max_threshold,Ustar_filter_type):     
    
    if 'Fc' in list_out:
        units="umol.m-2.s-1"
    elif ('Fe' or 'Fh' or 'Fg') in list_out:
        units="W.m-2"
    else:
        units=" "
    
    ###### User-set IO file locations ######
    
    print "Starting ANN gap filling"
    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/ANN"):
        os.mkdir(myBaseforResults+"/ANN")
    mypathforResults=myBaseforResults+"/ANN"  
    if not os.path.isdir(myBaseforResults+"/ANN/Fre"):
        os.mkdir(myBaseforResults+"/ANN/Fre")
    mypathforResults=myBaseforResults+"/ANN/Fre" 
    
    #We need to update VPD for input here so also need e and es
    # Calculate vapour pressure from absolute humidity and temperature
    #  Ah - absolute humidity, g/m3
    #  Ta - air temperature, C
  
    number_of_inputs=len(list_in)
    number_of_outputs=len(list_out)
    #startdate=dt.date(2008,7,1)
    #enddate=dt.date(2008,8,1)
    alllist=list_in + list_out
    #Here now for Re we can further restrict the traing data to be noctural and ustar filtered
    #So first create a series with day_night based on solar geometry
    
    #Get timezone info, call routines from external code
    global timezone
    currentdate="2013-06-01"
    timezone,InDstNow=TimeZone.get_timezone_info(latitude,longitude,currentdate)
    print "AskGEO TimZone offset (hrs):  ", timezone
 
    #Start with blank series
    New_combined['day_night']=np.nan
     
    def day_night(x):
	#Get date from data group
	currentdate= x.index[0]
	currentyear= x.index[0].year
	currentmonth= x.index[0].month
	currentday= x.index[0].day
	basedate=dt.datetime(1900, 1, 1,0,0,0,0)    
	delta=(currentdate-basedate).days  
	#Calculate Approximate Solar noon, call routines from external code
	#Call Solar_Calcs.solar_calculations(Date_input_excel,latitude,longitude,timezone)
	solar_sunrise,solar_noon_for_date,solar_sunset  =Solar_Calcs.solar_calculations(delta,latitude,longitude,timezone)
	#return a fraction.  Convert to decimal hours
	solar_sunrise=solar_sunrise*24
	solar_sunset=solar_sunset*24
	daystart_hour= int(solar_sunrise)
	daystart_minute=int((solar_sunrise-daystart_hour)*60)
	daystart_dt=dt.datetime(currentyear,currentmonth,currentday,daystart_hour,daystart_minute,0)
	dayend_hour= int(solar_sunset)
	dayend_minute=int((solar_sunset-dayend_hour)*60)
	dayend_dt=dt.datetime(currentyear,currentmonth,currentday,dayend_hour,dayend_minute,0)	

	x['day_night'][daystart_dt:dayend_dt]=1
	#Define evening as 3 hours after sunset.  Needed for van Gorsel approach
	d=dt.timedelta(hours=3)
	eveningend_dt=dayend_dt+d
	x['day_night'][dayend_dt:eveningend_dt]=2
	#Else fill remainder with night
	x['day_night'][x['day_night'].isnull()]=3	
	#print x
	#print x[dayend_dt:eveningend_dt]
	#print solar_sunset
	#print daystart_dt
	#print dayend_dt
	#print eveningend_dt	
	return x
    
    #For each day of each year run the function to calculate day or night
    New_combined=New_combined.groupby([lambda x: x.year,lambda x: x.month,lambda x: x.day]).apply(day_night)
    #New_combined.to_csv(myBaseforResults+'/'+'Final_combined_'+Site_ID+'.csv', sep=',')
    
    ###############################
    #Define the data to select
    ###############################
    #Can select period of 'night' (2) and or 'evening' (3)
    if evening==True:
	xnow=New_combined[(New_combined['day_night']==3)]
    else:
	xnow=New_combined[(New_combined['day_night']==2)]
    # Use the actual ustar_used column which is defined from the type of ustar thershold approach chosen in the config
    xnow=xnow[xnow['ustar']>xnow['ustar_used']]

	
    #Remove -ve and +ve spikes which are uptake and inconsistent with Fre.
    xnow=xnow[xnow['Fc']>min_threshold][xnow['Fc']<max_threshold]
    xnow=xnow[alllist]
    #Drop nans and missing values so that Good data only is used in the training
    xnow=xnow.dropna(how='any')	

    xarray=np.array(xnow.dropna().reset_index(drop=True))
    #Define inputs and targets for NN from DF
    inputs =  xarray[:, :number_of_inputs] #first 2 columns
    lastcolums=(-1*number_of_outputs)
    targets = xarray[:, lastcolums:] #last column
    
    # Generate standard layered network architecture and create network
    #different network architectures avaiable
    #conec = mlgraph((number_of_inputs,24,16,number_of_outputs))  # Creates standard multilayer network architecture
    conec = tmlgraph((number_of_inputs,6,4,number_of_outputs))  # Creates multilayer network full connectivity list
    #conec = imlgraph((number_of_inputs,24,16,number_of_outputs))  # Creates multilayer architecture with independent outputs
    
    net = ffnet(conec)
    
    print "TRAINING NETWORK..."
    #net.train_tnc(inputs, targets, maxfun = iterations, messages=1)
    try:
	net.train_rprop(inputs, targets, maxiter=iterations)
    except:
	net.train_tnc(inputs, targets, maxfun = iterations, messages=1)
    #net.train_momentum(inputs, targets, maxfun = iterations, messages=1)
    #net.train_genetic(inputs, targets, maxfun = iterations, messages=1)
    #net.train_cg(inputs, targets, maxfun = iterations, messages=1)
    #net.train_bfgs(inputs, targets, maxfun = iterations, messages=1)
    
    
    # Test network
    print "TESTING NETWORK..."
    output, regression = net.test(inputs, targets, iprint = 0)
    
    print "R-squared:           %s  " %str(regression[0][2])
    #print "max. absolute error: %s  " %str(abs( array(output).reshape( len(output) ) - array(targets) ).max())
    output, regress = net.test(inputs, targets)
    
    #Create array for results. Then loop through elements on the original data to predict the ANN value
    predicted=np.empty((len(xarray),number_of_outputs))
    observed=np.empty((len(xarray),number_of_outputs))
    for index,rowdata in enumerate(xarray):
        predicted[index]=net([rowdata[0:number_of_inputs]])
        observed[index]=np.array(rowdata[(-1.0*number_of_outputs)])
    
    ############################################    
    # Generate output and return new variables
    ############################################
 
    for index, item in enumerate(list_out):
        New_combined[ANN_label]=net.call(New_combined[list_in])[:,index]
    
    #TEST
    #New_combined.to_csv("E:/My Dropbox/Dropbox/Data_flux_data/Site data processing/HowardSprings/Advanced/test_assertion.csv") 
    

    #####################################################
    #  Plots 
    #####################################################
    #Plot time series of all 30 minute data
    mintimeseries_plot(mypathforResults,predicted,observed,regress,item, Site_ID,units,targets,output,ANN_label,index_str)
    #Plot regression of Tower versus ANN
    regressionANN2(mypathforResults,predicted,observed,regress,item, Site_ID,units,ANN_label,index_str)
    #Plot diurnals for every second month 6 graphs
    Doplots_diurnal(mypathforResults,New_combined,item, Site_ID,units,ANN_label,index_str)
    #Plot diurnals for every second month 6 graphs
    if frequency=="all" or frequency=="annual":
	Doplots_diurnal_monthly(mypathforResults,New_combined,item, Site_ID,units,ANN_label,index_str)    
    #Plot timeseries of monthly over all periods
    Doplots_monthly(mypathforResults,New_combined,item, Site_ID,units,ANN_label,index_str)
		    
    ###################################################
    #  File stuff
    ###################################################
    
 
    return (New_combined)

###########################################################################################################
##                 START MAIN CODE
###########################################################################################################
#Not efficient but have written a wrapper around to get the grouping.
#Start by dividing the data up into slices.  Here can be all, annual, monthly or a column with a categorical variable (i.e. stage of LULCC or burning)

def Fre_ANN_gapfill(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,Tower_Lat,Tower_Long,frequency,evening,min_threshold,max_threshold,Ustar_filter_type):     
    
    if Site_ID=="RDMF":
	#define dates of phases
	startdate=dt.date(2011,9,11)
	phase1_end=dt.date(2012,3,2)
	phase2_end=dt.date(2012,3,6)
	phase3_end=dt.date(2012,8,6)
	phase4_end=dt.date(2012,8,28)
	phase5_end=dt.date(2013,1,22)
	phase6_end=dt.date(2013,5,23)
	phase7_end=dt.date(2013,7,23)
	
	New_combined['RDMF_Phase']=''
	New_combined['RDMF_Phase'][startdate:phase1_end]="Phase1"
	New_combined['RDMF_Phase'][phase1_end:phase2_end]="Phase2"
	New_combined['RDMF_Phase'][phase2_end:phase3_end]="Phase3"
	New_combined['RDMF_Phase'][phase3_end:phase4_end]="Phase4"
	New_combined['RDMF_Phase'][phase4_end:phase5_end]="Phase5"
	New_combined['RDMF_Phase'][phase5_end:phase6_end]="Phase6"
	New_combined['RDMF_Phase'][phase6_end:phase7_end]="Phase7"

    #Define the variable to be created
    item='Fre'
    #create and reset the variables
    ANN_label=str(item+"_NN")
    New_combined[ANN_label]=np.nan
    ANN_label_all=str(item+"_NN_all")
    New_combined[ANN_label_all]=np.nan
    
    if frequency == 'annual':
	New_combined_grouped=New_combined.groupby([lambda x: x.year])
    
    elif frequency == 'monthly':
	New_combined_grouped=New_combined.groupby([lambda x: x.year,lambda x: x.month])
    elif frequency=="all":
	pass
    else:
	New_combined_grouped=New_combined.groupby([frequency])
    
    #Do the ANN for the grouped analysis first. Then do 'ALL' data.
    #Need to fill in blocks with the annual ANN for blocks where there may not be enough data points
    #start with new column.  #Create a new variable called '_NN'
    #Fill with nans first
    
   
    
    ANN_label=str(item+"_NN")
    New_combined[ANN_label]=np.nan
    #ANN_flag_label=str(item+"_NN_flag")
    #New_combined[ANN_flag_label]=np.nan
    
    print "MEAN ", New_combined[ANN_label].mean()
    #Only do for GROUPED data otherwise skip
    if frequency!="all":
	i=0
	for index,group in New_combined_grouped:
	    #Define a string for the index.  if it is monthly then need to make a string includes year and month 2 dimensions
	    #Otherwise just use the DF index 1 dimension.
	    if frequency=='monthly':
		DFyear=int(index[0])
		DFmonth=int(index[1])  
		index_str=str(DFyear*1000+DFmonth)
	    else:
		index_str=str(index)
	    print "Performing ANN on index ",index_str
	    
	    #Here try ANN on the block.  If the block doesnt have enough values it will return an error and that 
	    #block is effectively skipped and will be blank
	    try:
		temp = Fre_ANN_gapfill_func(myBaseforResults,group,Site_ID,list_in,list_out,iterations,Tower_Lat,Tower_Long,index_str,ANN_label,frequency,evening,min_threshold,max_threshold,Ustar_filter_type)
	
	    except:
		print "Something wrong pass back unchanged for index", index_str
		temp = group
		
		pass
	    #after the first block then start to concatenate parts together
	    if i==0:
		temp_concat=temp
	    else:
		temp_concat=pd.concat([temp_concat,temp],axis=0)
	    i=i+1
    
    
    #NOW do for  'ALL' data.
    #Need to fill in blocks with the annual ANN for blocks where there may not be enough data points
    #start with new column.  #Create a new variable called '_NN'  #Fill with nans first
    item='Fre'
    
    #Define a string for the index.
    index_str='all'
    print "Performing ANN on index ",index_str
       
    #If using the grouped data then use the concatenated file otherwise if using ALL then use the original DF
    if frequency=="all":
	temp_concat = Fre_ANN_gapfill_func(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,Tower_Lat,Tower_Long,index_str,ANN_label_all,frequency,evening,min_threshold,max_threshold,Ustar_filter_type)
    else:
	temp_concat = Fre_ANN_gapfill_func(myBaseforResults,temp_concat,Site_ID,list_in,list_out,iterations,Tower_Lat,Tower_Long,index_str,ANN_label_all,frequency,evening,min_threshold,max_threshold,Ustar_filter_type)
    
    #Now put together the results from the groups and the all ANN
    #If user choose all data then simply copy 'all' column to the NN column
    #Otherwise keep column from groups and check for nans and fill them with the 'all' column
    if frequency == 'all':
	temp_concat[ANN_label]=temp_concat[ANN_label_all]
    else:
	temp_concat[ANN_label][temp_concat[ANN_label].isnull()]=temp_concat[ANN_label_all]
	
    #Clean up columns
    #del temp_concat[ANN_label_all]
    
   
    print"Completely Finished Fre ANN outputs at "+ Site_ID  
    
    return temp_concat