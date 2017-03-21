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

def Output_files(New_combined,myBaseforResults,Site_ID,versionID):
    New_combined_grouped=New_combined.groupby([lambda x: x.year])
    
    for year_index in New_combined_grouped:
        print year_index[0]
        
        print "Starting output for ARM"
        #Check for place to put results - does it exist? If not create
        if not os.path.isdir(myBaseforResults):
            os.mkdir(myBaseforResults)
        #Then subdirectories
        if not os.path.isdir(myBaseforResults+"/ARM"):
            os.mkdir(myBaseforResults+"/ARM")
        mypathforResults=myBaseforResults+"/ARM/"
        
        #Calculate RH_con
        New_combined['RH_Con']=metfuncs.RHfromabsolutehumidity(New_combined['Ah_Con'],New_combined['Ta_Con'])
        
        #Convert VPD in kPa to hPa.
        #We need to update VPD for input here so also need e and es
        # Calculate vapour pressure from absolute humidity and temperature
        #  Ah - absolute humidity, g/m3
        #  Ta - air temperature, C
        New_combined['VPD_hPa_Con']=(metfuncs.es(New_combined['Ta_Con']))-(metfuncs.vapourpressure(New_combined['Ah_Con'],New_combined['Ta_Con']))/10
        REddyProc_DF=New_combined[['Ah','Cc_7500_Av','Cc','eta','Fa','Fc','Fe','Fg','Fh','Fld','Flu','Fm','Fn','Fsd','Fsu','ps','Sws','Sws40cm','Sws100cm','Ta','theta','Ts','ustar','Ws_CSAT','Wd_CSAT']]
              
        #The date/time components are separated into columns. E.g. first column: julian day, second column: decimal hour. 
        #Possible date formats are indicated in the input form. Never use an hour of 24 with the time 
        #format 'year', 'month', 'day', 'hour', 'minute' (use 0 instead). Hour '0' is interpreted as first hour of the day, 
        #i.e. when you have transition from one to another it must be like (day, 23 --> day+1, 0) not like (day, 23 --> day, 0),
        #because then the data set is not chronological (this misunderstanding happened before).
        
        #REddyProc_DF['DT1','DT2','DT3','DT4','DT5','DT5']=REddyProc_DF.index.timetuple()
        REddyProc_DF['DTcopy']=REddyProc_DF.index
        
        REddyProc_DF['Day']=REddyProc_DF['DTcopy'].apply(lambda x: int(x.strftime('%j')))
        REddyProc_DF['Year']=REddyProc_DF['DTcopy'].apply(lambda x: int(x.strftime('%Y')))
        REddyProc_DF['Hour']=REddyProc_DF['DTcopy'].apply(lambda x: float(x.strftime('%H'))+(float(x.strftime('%M'))/60))
        
        #Select current year of yaer only
        REddyProc_DF=REddyProc_DF[REddyProc_DF['Year']==year_index[0]]
        
        #Calculate some things for plots
        n_datapoints=len(REddyProc_DF)
        startdate= REddyProc_DF.index[0]
        enddate= REddyProc_DF.index[n_datapoints-1]
        print n_datapoints,startdate,enddate
        
        newline1="TIMESTAMP,Merged from Ah_HMP_23m Ah_7500_Av Ah_HMP_2m,CO2 concentration average,Merged from Cc_7500_Av converted to umol/mol,Horizontal rotation angle,Available energy using Fn Fg,CO2 flux rotated to natural wind coordinates WPL corrected Fc converted to umol/m2/s,Latent heat flux rotated to natural wind coordinates WPL corrected Fe,Element-wise average of series Fg_8cma Fg_8cmb Fg_8cmc Fg_8cmd Soil heat flux corrected for storage,Sensible heat flux rotated to natural wind coordinates Fh rotated and converted from virtual heat flux,Down-welling long wave,Up-welling long wave,Momentum flux rotated to natural wind coordinates,Merged from Fn_KZ Fn_NR	Down-welling short wave	Up-welling short wave,Air pressure standard deviation,Element-wise average of series Sws_10cma Sws_10cmb,Soil water fraction sensor 2a,Soil water fraction sensor 3a,Merged from Ta_HMP_23m Ta_CSAT Ta_HMP_2m,Vertical rotation angle,Element-wise average of series Ts_8cma,Friction velocity rotated to natural wind coordinates,ustar filtered for low turbulence conditions (<0.25),Wind speed,Wind direction"
        newline2="TIMESTAMP,g/m3,mg/m3,umol/mol,deg,W/m2,umol/m2/s,W/m2,W/m2,W/m2,W/m2,W/m2,kg/m/s2,W/m2,W/m2,W/m2,kPa,frac,frac,frac,C,deg,C,m/s,m/s,m/s,deg"
        newline3= "TIMESTAMP,Ah,Cc_7500_Av,Cc,eta,Fa,Fc,Fe,Fg,Fh,Fld,Flu,Fm,Fn,Fsd,Fsu,ps,Sws,Sws_40cma,Sws_100cma,Ta,theta,Ts,ustar,ustar_filtered,Ws_CSAT,Wd_CSAT"
         
        columns_out  = ['DTcopy','Ah','Cc_7500_Av','Cc','eta','Fa','Fc','Fe','Fg','Fh','Fld','Flu','Fm','Fn','Fsd','Fsu','ps','Sws','Sws40cm','Sws100cm','Ta','theta','Ts','ustar','ustar','Ws_CSAT','Wd_CSAT']             
        
        output_temp_filename=mypathforResults+'/ARM_temp_'+Site_ID+'_'+str(year_index[0])+'_'+versionID + '.csv'
        output_filename=mypathforResults+'/ARM_'+Site_ID+'_'+str(year_index[0])+'_'+versionID +'.csv'
        
        
        
        REddyProc_DF.to_csv(output_temp_filename, na_rep='-9999', float_format='%.3f', cols=columns_out, header=False, index=False, index_label=None, mode='w')
        
        #Now add another line with units
        #Open txt file
        with open(output_temp_filename) as infile:
            with open(output_filename,"w") as outfile:
		outfile.write(newline1+"\n")
		outfile.write(newline2+"\n")
		outfile.write(newline3+"\n")
		for i,line in enumerate(infile):
		    outfile.write(line)
        
        os.remove(output_temp_filename)
   

    #####################
    # Finish up
    ######################
    
    print "FINISHED writing out files for use for ARM "
    
    
    
    

