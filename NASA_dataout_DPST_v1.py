############################################################################
#   
# Programmed by Jason (Dec 1, 2012)
############################################################################

import pandas as pd
import datetime as dt
import calendar
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
import meteorologicalfunctions as metfuncs



####################
##START MAIN CODE
####################

def Output_files(New_combined,myBaseforResults,Site_ID,versionID,Ws_variable_name):
    #Do any calculations on the whole datasewt before grouping
    #Calculate RH_con
    New_combined['RH_Con']=metfuncs.RHfromabsolutehumidity(New_combined['Ah_Con'],New_combined['Ta_Con'])
    #Convert VPD in kPa to hPa.
    #We need to update VPD for input here so also need e and es
    # Calculate vapour pressure from absolute humidity and temperature
    #  Ah - absolute humidity, g/m3
    #  Ta - air temperature, C
    New_combined['VPD_kPa_Con']=(metfuncs.es(New_combined['Ta_Con']))-(metfuncs.vapourpressure(New_combined['Ah_Con'],New_combined['Ta_Con']))   

    #Do timestamp operations
    #Make a copy of timestamp to the df
    #Take mean first (equal to mid day) then convert to DOY, day, month and year
    New_combined['DTcopy']=New_combined.index
    New_combined['Year']=New_combined['DTcopy'].apply(lambda x: int(x.strftime('%Y')))
    
 

    #Group DF by year
    New_combined_grouped=New_combined.groupby([lambda x: x.year])
    
    for year_index in New_combined_grouped:
        print year_index[0]
        
        print "Starting output for NASA"
        #Check for place to put results - does it exist? If not create
        if not os.path.isdir(myBaseforResults):
            os.mkdir(myBaseforResults)
        #Then subdirectories
        if not os.path.isdir(myBaseforResults+"/NASA_out"):
            os.mkdir(myBaseforResults+"/NASA_out")
        mypathforResults=myBaseforResults+"/NASA_out/"
        
        #Subset the DF to make it easier
	#WD removed here as its not required for NASA and for some sites SD variable names are not standard
	if Ws_variable_name=="Ws_CSAT":
	    REddyProc_DF=New_combined[['DTcopy','Year','Ah_Con','Cc','eta','Fa','Fc_ustar','GPP_Con','Fre_Con','Fe_Con','Fg_Con','Fh_Con','Fld_Con','Flu_Con','Fm','Fn_Con','Fsd_Con','Fsu_Con','ps_Con','Precip_Con','Sws_Con','Ta_Con','Ts_Con','ustar','Ws_CSAT_Con','RH_Con','VPD_kPa_Con','Ah_Con_QCFlag','Fc_Con_QCFlag','Fe_Con_QCFlag','Fg_Con_QCFlag','Fh_Con_QCFlag','Fld_Con_QCFlag','Flu_Con_QCFlag','Fn_Con_QCFlag','Fsd_Con_QCFlag','Fsu_Con_QCFlag','ps_Con_QCFlag','Precip_Con_QCFlag','Sws_Con_QCFlag','Ta_Con_QCFlag','Ts_Con_QCFlag']]
	else:
	    REddyProc_DF=New_combined[['DTcopy','Year','Ah_Con','Cc','eta','Fa','Fc_ustar','GPP_Con','Fre_Con','Fe_Con','Fg_Con','Fh_Con','Fld_Con','Flu_Con','Fm','Fn_Con','Fsd_Con','Fsu_Con','ps_Con','Precip_Con','Sws_Con','Ta_Con','Ts_Con','ustar','Ws_Con','RH_Con','VPD_kPa_Con','Ah_Con_QCFlag','Fc_Con_QCFlag','Fe_Con_QCFlag','Fg_Con_QCFlag','Fh_Con_QCFlag','Fld_Con_QCFlag','Flu_Con_QCFlag','Fn_Con_QCFlag','Fsd_Con_QCFlag','Fsu_Con_QCFlag','ps_Con_QCFlag','Precip_Con_QCFlag','Sws_Con_QCFlag','Ta_Con_QCFlag','Ts_Con_QCFlag']]
	
        #Select current year of yaer only
        REddyProc_DF=REddyProc_DF[REddyProc_DF['Year']==year_index[0]]
        
        #Calculate some things for plots
        n_datapoints=len(REddyProc_DF)
        startdate= REddyProc_DF.index[0]
        enddate= REddyProc_DF.index[n_datapoints-1]
        print n_datapoints,startdate,enddate
        
        #Calculate the DAILY means/sums from the half hourly data

	tempDF_mean=REddyProc_DF.groupby(lambda x : x.dayofyear).mean().add_suffix('_mean')
	tempDF_sum=REddyProc_DF.groupby(lambda x : x.dayofyear).sum().add_suffix('_sum')
	
	tempDF=tempDF_mean.join(tempDF_sum,how='left') 
			
	#Add QC counts to the means DF
	#Good QC value not gap filled is 1.  Get sall values ==1 then do a count.  Divide by  48 for 48 half hour periods in the day
	tempDF['Rn_qc']=REddyProc_DF['Fn_Con_QCFlag'][REddyProc_DF['Fn_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['Rn_qc'].fillna(value=0,inplace=True)
	tempDF['Rs_qc']=REddyProc_DF['Fsd_Con_QCFlag'][REddyProc_DF['Fsd_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['Rs_qc'].fillna(value=0,inplace=True)
	tempDF['Ta_qc']=REddyProc_DF['Ta_Con_QCFlag'][REddyProc_DF['Ta_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['Ta_qc'].fillna(value=0,inplace=True)
	tempDF['VPD_qc']=REddyProc_DF['Ah_Con_QCFlag'][REddyProc_DF['Ah_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['VPD_qc'].fillna(value=0,inplace=True)
	tempDF['Ts_qc']=REddyProc_DF['Ts_Con_QCFlag'][REddyProc_DF['Ts_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['Ts_qc'].fillna(value=0,inplace=True)
	tempDF['NEE_qc']=REddyProc_DF['Fc_Con_QCFlag'][REddyProc_DF['Fc_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['NEE_qc'].fillna(value=0,inplace=True)
	tempDF['GPP_qc']=REddyProc_DF['Fc_Con_QCFlag'][REddyProc_DF['Fc_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['GPP_qc'].fillna(value=0,inplace=True)
	tempDF['Reco_qc']=REddyProc_DF['Fc_Con_QCFlag'][REddyProc_DF['Fc_Con_QCFlag']==1.].groupby(lambda x : x.dayofyear).count()/48
	tempDF['Reco_qc'].fillna(value=0,inplace=True)	
	
	#add a site lable to columns
	tempDF['Site_ID']=Site_ID
	
	tempDF['DTmean']=REddyProc_DF['DTcopy'].groupby(lambda x : x.dayofyear).min()
	tempDF['Day']=tempDF['DTmean'].apply(lambda x: int(x.strftime('%d')))
	tempDF['Month']=tempDF['DTmean'].apply(lambda x: int(x.strftime('%m')))        
	tempDF['Year']=tempDF['DTmean'].apply(lambda x: int(x.strftime('%Y')))
	# Jan the 1st is day 1           
	#tempDF['DOY'] = (tempDF['DTmean'] - dt.datetime(year_index[0], 1, 1))
	
	tempDF['DOY'] = tempDF['DTmean'].apply(lambda x: int(x.strftime('%j')))

	#Do conversions for Carbon variables (convert from umol to g C for NASA)
	tempDF['Fc_ustar_mean']=tempDF['Fc_ustar_mean']*60*60*24/1000000*12
	tempDF['GPP_Con_mean']=tempDF['GPP_Con_mean']*60*60*24/1000000*12
	tempDF['Fre_Con_mean']=tempDF['Fre_Con_mean']*60*60*24/1000000*12
		
	#Do conversions for Radiation variables (convert Wm-2 to MJ m-2 day-1)
	tempDF['Fsd_Con_mean']=tempDF['Fsd_Con_mean']*60*60*24/1000000
	tempDF['Fn_Con_mean']=tempDF['Fn_Con_mean']  *60*60*24/1000000

	newline2="ID, Year, Mo, Day, DOY, Rn_f, Rn_qc, Rs_f, Rs_qc, Ta, Ta_qc, VPD, VPD_qc, Ts_f, Ts_qc, PREC, SWC, NEE, NEE_qc, GPP, GPP_qc, Reco, Reco_qc, PRESS, SNOWD"
	newline3="-, -, -, -, -, MJ m-2 day-1, -, MJ m-2 day-1, -, oC, -, kPa, -, oC, -, mm day-1, m3/m3, gC m-2 day-1, -, gC m-2 day-1, -, gC m-2 day-1, -, MPa day-1, mm"
	columns_out  = ['Site_ID','Year','Month','Day','DOY', 'Fn_Con_mean','Rn_qc','Fsd_Con_mean','Rs_qc', 'Ta_Con_mean', 'Ta_qc', 'VPD_kPa_Con_mean', 'VPD_qc', 'Ts_Con_mean', 'Ts_qc', 'Precip_Con_sum', 'Sws_Con_mean', 'Fc_ustar_mean', 'NEE_qc', 'GPP_Con_mean', 'GPP_qc','Fre_Con_mean', 'Reco_qc', 'ps_Con_mean']
	output_temp_filename=mypathforResults+'/NASA_SMAP_temp_'+Site_ID+'_'+str(year_index[0])+'_'+versionID +'.csv'
	output_filename=mypathforResults+'/NASA_SMAP_'+Site_ID+'_'+str(year_index[0])+'_'+versionID +'.csv'
	tempDF[columns_out].to_csv(output_temp_filename, na_rep='-9999', float_format='%.3f', header=False, index=False, index_label=None, mode='w')

        
        #Now add another line with units
        #Open txt file
        with open(output_temp_filename) as infile:
	    with open(output_filename,"w") as outfile:
		#outfile.write(newline1+"\n")
		outfile.write(newline2+"\n")
		outfile.write(newline3+"\n")                   
                for i,line in enumerate(infile):
                        outfile.write(line)
        
        os.remove(output_temp_filename)
        
    #####################
    # Finish up
    ######################
    
    print "FINISHED writing out files for use for NASA "

    
    

