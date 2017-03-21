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
        
        print "Starting output for EddyProc MPI online tool"
        #Check for place to put results - does it exist? If not create
        if not os.path.isdir(myBaseforResults):
            os.mkdir(myBaseforResults)
        #Then subdirectories
        if not os.path.isdir(myBaseforResults+"/REddyProc"):
            os.mkdir(myBaseforResults+"/REddyProc")
        mypathforResults=myBaseforResults+"/REddyProc/"
        
        #Calculate RH_con
        New_combined['RH_Con']=metfuncs.RHfromabsolutehumidity(New_combined['Ah_Con'],New_combined['Ta_Con'])
        
        #Convert VPD in kPa to hPa.
        #We need to update VPD for input here so also need e and es
        # Calculate vapour pressure from absolute humidity and temperature
        #  Ah - absolute humidity, g/m3
        #  Ta - air temperature, C
        New_combined['VPD_hPa_Con']=(metfuncs.es(New_combined['Ta_Con']))-(metfuncs.vapourpressure(New_combined['Ah_Con'],New_combined['Ta_Con']))/10
        REddyProc_DF=New_combined[['Fc','Fe','Fh','Fg','Ta_Con','Ts_Con','RH_Con','VPD_hPa_Con','ustar']]
              
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
        
        
        header_names = ['Year','Day','Hour','NEE','LE','H' ,'Rg','Tair'  ,'Tsoil' ,'rH',    'VPD',        'Ustar']
        columns_out  = ['Year','Day','Hour','Fc', 'Fe','Fh','Fg','Ta_Con','Ts_Con','RH_Con','VPD_hPa_Con','ustar']
        
        newline1='Year \t DoY \t Hour \t NEE \t LE \t H \t Rg \t Tair \t Tsoil \t rH \t VPD \t Ustar'
        newline2=" -- \t -- \t -- \t umolm-2s-1 \t Wm-2 \t Wm-2 \t Wm-2 \t degC \t degC \t % \t hPa \t ms-1"
        
        #newline1='Year,Day,Hour,NEE,LE,H,Rg,Tair,Tsoil,rH,VPD,Ustar'
        #newline2="--,--,--,umolm-2s-1,Wm-2,Wm-2,Wm-2,degC,degC,%,hPa,ms-1"
        
        output_temp_filename=mypathforResults+'/REddyProc_temp_'+Site_ID+'_'+str(year_index[0])+'_'+versionID +'.txt'
        output_filename=mypathforResults+'/REddyProc_'+Site_ID+'_'+str(year_index[0])+'_'+versionID +'.txt'
        
        
        
        REddyProc_DF.to_csv(output_temp_filename, sep='\t', na_rep='-9999', float_format='%.3f', cols=columns_out, header=False, index=False, index_label=None, mode='w')
        
        #Now add another line with units
        #Open txt file
        with open(output_temp_filename) as infile:
            with open(output_filename,"w") as outfile:
                outfile.write(newline1+"\n")
                outfile.write(newline2+"\n")  
                for i,line in enumerate(infile):
                    outfile.write(line)
        
        os.remove(output_temp_filename)
        
    #####################
    # Finish up
    ######################
    
    print "FINISHED writing out files for use in EddyProc MPI online tool "
    
    
    
    

