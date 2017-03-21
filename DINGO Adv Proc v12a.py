############################################################################
# Dynamic INtegrated Gap-filling and partitioning for OzFlux (DINGO)
############################################################################
#The overall approach used in DINGO is to take the L3 OzFluxQC data, which has gaps from data processing 
#(data excluded due to values out of range, spike detection or manual exclusion of date and time ranges) 
#or from site issues (instrument or power failure, herbivores, fire, eagles nests, cows, lightning, PI on 
#sabbatical, etc.) and gap-fill and partition the data using a variety of data sources (Fig. 1).  DINGO is 
#programmed in python 2.7 and is currently at version 12a and publically available on GitHub (https://github.com/jberinge/DINGO12). 
#It should be noted that v13 is scheduled to be released in July 2016 and will include uncertainty 
#and footprint analysis .   It is designed
#to work with OzFlux data produced in NetCDF format by the OzFluxQC (Isaac et al., 2016) and draws on 
#Australian AWS data but could be adapted for other data sources across other flux sites.  The primary
#interface for the user is through a text based control file that has information on site characteristics
#(name, latitude and longitude, the frequency of the flux measurements (30 or 60 minutes) and elevation),
#file paths (to the OzFluxQC NetCDF files and other ancillary data inputs), data processing options and
#data plotting and output formats.  In general, prior to the processing steps below, any gaps in fluxes
#or meteorological quantities of less than two hours are filled by DINGO using linear interpolation.  
#
#DINGO is documented in  BGD paper 2016.
#
# Programmed by Jason Beringer updated 5/5/2016
############################################################################

import pandas as pd
import os
import datetime as dt
import pickle
from configobj import ConfigObj
import Tkinter, tkFileDialog
import numpy as np

#Import custom code modules required for adv processing
import AWS_fetch_v3a as Fetch_AWS
import Global_Solar_fetch_v1 as Fetch_Global_Solar
import FLUX_nc_fetch_v5 as Fetch_Flux
#import Met_construct_v4 as Met_construct
import Met_construct_v4_linux as Met_construct
 
import BIOS_nc_fetch_v5 as CABLE_nc_fetch
import CABLE_gapfil_var_v4a as CABLE_gapfil_var
import MODIS_ingest_and_interpolate_with_flux_tower_NC_file_v5b as MODIS
import Rain_gapfill_v6 as Rainstuff
import IncomingSW_gapfill as Shortwave_incoming
import OutgoingSW_gapfill as Shortwave_outgoing
import IncomingLW_gapfill as Longwave_incoming
import OutgoingLW_gapfill as Longwave_outgoing
import Net_radiation_v1 as Net_radiation
import Light_response_curve_v25c as Fre_Lasslop

import FFNET_v5a as ANN
import ustar_filtering_Reichsteinv2a as ustar
import FFNET_Fre_v5b as Fre_ANN
import GPP_calc_and_ustar_sensitivity_v2 as GPP
import ISAACfingerprintv2 as  Fingerprints
import Diagnostics_Results_v5e as Diagnostics
import QCCPD_DINGO_v1b as CPD

#Modules for data formatting and output
import REddyProc_dataout_v2a as REddyProc
import ARM_dataout_v1 as ARMProc
import WAVES_dataout_DPST_v1 as WAVESProc
import NASA_dataout_DPST_v1 as NASAProc

def dataframe_check(Dataframe, FluxFreq):
    #Check dataframe for duplicates, pad as necessary and sort
    Dataframe.sort(inplace=True)
    Dataframe["index"] = Dataframe.index
    Dataframe.drop_duplicates('index', take_last=True, inplace=True)
    del Dataframe["index"]	
    Dataframe=Dataframe.asfreq(FluxFreq, method=None)      
    return Dataframe

def main():
    #Set version ID
    versionID="v12a"
    #Define string of Rainfall name.
    Rain_label_variable_to_fill = 'Precip'
    
    #Open Config file to input USER configuration data
    print "Opening Configuration file"
    #Prepare to open window to prompt for filename of Config file
    root = Tkinter.Tk(); root.withdraw()
    CFname = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    #Returns CF object that we can get info from later
    cf = ConfigObj(CFname) 
    #Input site details
    Site_ID=cf['Site']['Site_ID']
    Tower_Lat = float(cf['Site']['Tower_Lat'])
    Tower_Long = float(cf['Site']['Tower_Long'])
    #Read in frequency of fluxdata
    FluxFreq = cf['Site']['FluxFreq']
    altitude = cf['Site']['Altitude']
    Ws_variable_name = cf['Options']['Ws_variable_name']
    Ustar_filter_type = cf['Options']['Ustar_filter_type']
    Ustar_Barr_calculate= cf['Options']['Ustar_Barr_calculate']
    Ustar_Barr_bootstraps= int(cf['Options']['Ustar_Barr_bootstraps'])
        
    print "Start advanced processing for "+Site_ID
    print "Site coordinates latitude "+str(Tower_Lat)+" and longitude "+str(Tower_Long)
    
    #File stuff here - Get data from dropbox folder etc
    mypathforFluxdata=cf['Files']['mypathforFluxdata']+Site_ID
    myBaseforResults=cf['Files']['mypathforFluxdata']+Site_ID+ "/Advanced_"+versionID
    mypathforAWSdata=cf['Files']['mypathforAWSdata']
    mypathforGlobalSolardata=cf['Files']['mypathforGlobalSolardata']
    mypathforAWAPdata=cf['Files']['mypathforAWAPdata']
    FileName_AWAPData=cf['Files']['FileName_AWAPData']
    BIOS2filename=cf['Files']['BIOS2filename']
    inputMODIS_base = cf['Files']['inputMODIS_base'] 
    BoMfilename="BoMstations.txt"
    
    #Get option data from config files
    #either 'all','annual','monthly'
    corr_freq=cf['Options']['corr_freq'] 
    #The Bom AWS has time stamp starting at the hour and Loggers time stamp at end of period. So we need to offset by 30 minutes
    if cf['Options']['OffsetAWS30min']=='Yes':
	OffsetAWS30min=True
    
    ##================================================================
    ##                Main code started here
    ##================================================================
    #Check for place to put resuolts - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)

    # Start by getting the AWS flux varibles first.
    # Will return DF with 3 closest AWS station that are current and have 30 minute data.
    # Will also return a list of top three AWS BoM station ID
    
    DFfilename=myBaseforResults+'/Advanced_processed_data_'+Site_ID+'_'+versionID+'.df'
    CSVfilename=myBaseforResults+'/Advanced_processed_data_'+Site_ID+'_'+versionID+'.csv'
    
    #next Get the Flux tower nc file
    if cf['Options']['GetL3FluxTowerData']=='Yes':   
	print "Getting the Flux tower nc file"
	#Get all variables
	Fluxfilepath=mypathforFluxdata+'/'+cf['Files']['FLUXfilename']
	#Call the module to retrieve the L3 nc file and return a DF and a dictionary "Flux_Variable_Units" containing units from metadata
	FLUXDataframe, Flux_Variable_Units = Fetch_Flux.fetch_flux(Fluxfilepath,FluxFreq,Site_ID)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	FLUXDataframe = dataframe_check(FLUXDataframe, FluxFreq)
	# save dataframe. 
	FLUXDataframe.to_csv(myBaseforResults+'/FLUXDataframe_'+Site_ID+'.csv')  
	FLUXDataframe.to_pickle(myBaseforResults+'/FLUXDataframe_'+Site_ID+'.df') 
  
	
    if cf['Options']['GetBoMAWSfiles']=='Yes':
	print "Starting processing of BoM AWS files"
	FLUXDataframe= pd.read_pickle(myBaseforResults+'/FLUXDataframe_'+Site_ID+'.df')
	#Call routines to injest the gridded daily global solar data from BoM subscription
	#save it to the flux data frame
	Flux_plus_solar = Fetch_Global_Solar.get_Global_Solar_data(BoMfilename,Site_ID,Tower_Lat,Tower_Long,mypathforGlobalSolardata,myBaseforResults,OffsetAWS30min,FLUXDataframe,FluxFreq)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	Flux_plus_solar = dataframe_check(Flux_plus_solar, FluxFreq)	
	#Save DF
	Flux_plus_solar.to_pickle(myBaseforResults+'/FLUXDataframe_'+Site_ID+'.df')

	AWS_combined, bestAWS_ID = Fetch_AWS.get_AWS_data(BoMfilename,Site_ID,Tower_Lat,Tower_Long,mypathforAWSdata,myBaseforResults, OffsetAWS30min,Flux_plus_solar,FluxFreq)
	
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	AWS_combined = dataframe_check(AWS_combined, FluxFreq)
	
	# save dataframe and CSV files
	AWS_combined.to_pickle(myBaseforResults+'/AWS_combined_'+Site_ID+'.df') 
	AWS_combined.to_csv(myBaseforResults+'/'+'AWS_combined_'+Site_ID+'.csv', sep=',')
	pickle.dump(bestAWS_ID, open((myBaseforResults+'/AWS/bestAWS_ID_'+str(Site_ID)), 'wb'))
	
	#Combine into single file, UUse 'inner' here so that it effectively cuts off any data that is not matching
	print "Combining AWS and FLUX files"
	ALL_combined2=Flux_plus_solar.join(AWS_combined,how='left')    
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	ALL_combined2 = dataframe_check(ALL_combined2, FluxFreq)  
	# Temp statement to save dataframe. 
	ALL_combined2.to_pickle(DFfilename) 
	ALL_combined2.to_csv(CSVfilename)
	
    #Do correlation and plots for AWS variables
    #-------------------------------------------
    #This will return a dataframe with the CONSTRCUCTED variable  'XXX_Con' and the FLAG 'XXX_Con' and  
    #variable that is correlated variable 'XXX__Corr' and 
    #At the moment the variables lists in the Bom AWS and Flux files are different. Often just a capital (or not)
    #We need to specifiy them here and pass to the constructor
    if cf['Options']['ConstructMetfromBoMAWS']=='Yes':   
	# Load in previous data step data
	ALL_combined= pd.read_pickle(DFfilename)   
	FLUXDataframe= pd.read_pickle(myBaseforResults+'/FLUXDataframe_'+Site_ID+'.df') 
	bestAWS_ID = pickle.load(open((myBaseforResults+'/AWS/bestAWS_ID_'+str(Site_ID)), 'rb'))
	
	#Define the variable names to be processed.  These are static but may change with Level of input data or previously processed data
	#Change if necessary
	AWS_variables=['Ta','WS','P','Ah']
	Flux_variables=['Ta',Ws_variable_name,'ps','Ah']
	construct=list(xrange(len(AWS_variables)))
	for index, items in enumerate(AWS_variables):
	    print "Calling construct for variable "+items
	    construct[index]=Met_construct.construct_data(ALL_combined, Flux_variables[index] ,AWS_variables[index], bestAWS_ID, Site_ID,corr_freq,myBaseforResults)
	#then join the returned data  
	#Concatenate all the series together
	ConstructFrame=pd.concat((construct[:]),axis=1)		
	New_combined=FLUXDataframe.join(ConstructFrame,how='left')
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)  
    
    #Call routines to ingest CABLE LSM data.  Returns dataframe with CABLE variables
    #CABLE data processed using previously defined levels for soil layers
    #Also hard coded the start date of the CABLE file as 1/1/1990.  This can be changed in the routines
    #Need to specify the CABLE file name at the start of this code whcih is passed to the routine
    #Variables 'Sws_CABLE' and 'Ts_CABLE' are returned which are the relevant layers of CABLE averaged to compare with Tower data
    #Version 11 of main code now uses BIOS2 input
    if cf['Options']['Fetch_CABLE_data_NC']=='Yes':   
    	# Load in previous data step data
	New_combined=pd.read_pickle(DFfilename) 
	#Call fetch CABLE data routines
	New_combined=CABLE_nc_fetch.fetch_CABLE(New_combined,myBaseforResults,BIOS2filename,Site_ID,Tower_Lat,Tower_Long)  
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq) 
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)    
	
    #Call routine to gapfill Ts.  Will return variable 'XXX_con' combining Tower and Cable.
    #Flag is 1 for good Tower data and 99 for Gap filled
    if cf['Options']['Construct_gapfill_Ts_Swc']=='Yes':   
	# Load in previous data step data
	New_combined=pd.read_pickle(DFfilename)   	
	#Define variables to fetch and gap fill from CABLE data.  Can be modified
	variables_to_fill = ['Ts', 'Sws']
	for variable_to_fill in variables_to_fill:
	    New_combined=CABLE_gapfil_var.CABLE_gapfill(variable_to_fill,New_combined,myBaseforResults,Site_ID)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq) 
	#Save dataframe
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename) 	  
	
    #Call routine to gapfill Rainfall.  Will return the DF with variables added, including 'XXX_con' combining Tower and Cable
    #Flag is 1 for good Tower data and 99 for Gap filled using AWS and 100 for gap filled using CABLE/AWAP gridded met.
    #In some L3 files the label is Precip and others Rain.  Choose accordingly.  Also if you want to set to include 
    #AWAP/CABLE rainfall in the comparison then select True (default)
    if cf['Options']['Gapfill_Rainfall']=='Yes':   
	# Load in previous data step data
	New_combined=pd.read_pickle(DFfilename)  
	bestAWS_ID = pickle.load(open((myBaseforResults+'/AWS/bestAWS_ID_'+str(Site_ID)), 'rb'))
	#Do you want to use CABLE data for rainfall or just AWS data.  Using CAble is good esp for remote sites
	#Currently cable rainfall only up to 2010. So set to False if required
	if cf['Options']['CABLEdata_for_rainfallgap_too']=='Yes': 
	    cable_rain=True
	else:
	    cable_rain=False
	New_combined=Rainstuff.Rain_gapfill(Rain_label_variable_to_fill,New_combined,myBaseforResults,Site_ID,bestAWS_ID,cable_rain,FluxFreq,Tower_Lat,Tower_Long,versionID)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  
	#to_pickle new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)    
    #elif cf['Options']['Gapfill_Rainfall']=='No':   
	## Just copy observed to Con and Corr columns
	#New_combined=pd.read_pickle(DFfilename)  
	#New_combined['Precip_Con']=np.nan
	#New_combined['Precip_Corr']=np.nan
	#New_combined['Precip_Con_QCFlag']=np.nan
	#New_combined['Rainf_CABLE_mm']=np.nan 
	##to_pickle new data
	#New_combined.to_pickle(DFfilename) 
	#New_combined.to_csv(CSVfilename)   
	
    #Call MODIS ingestion routines.  Will return all MODIS products into the DF that are listed under the MODIS directory    
    if cf['Options']['Ingest_MODIS']=='Yes':   
	# Load in previous data step data
	New_combined=pd.read_pickle(DFfilename) 
	#Call MODIS routines
	MODIS_key=cf['Options']['MODIS_key']
	New_combined=MODIS.MODIS_ingest_interp(New_combined,myBaseforResults,inputMODIS_base,Site_ID,MODIS_key,FluxFreq)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  	
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename) 
    
    if cf['Options']['Gapfill_radiation']=='Yes':      
	#Run Ians code to gapfill Fsd, Downscaling routine for estimating half-hourly radiation from AWAP daily estimates
	#Equation of time was excluded for reasons of simplicity (and doesn't have large effect)
	#Outputs cloudiness index to be used as input for long wave incoming estimation routine
	# Load in previous data step data
	New_combined=pd.read_pickle(DFfilename) 
	#Define the exact string for the variable to Fill
	variable_to_fill = 'Fsd'    
	New_combined=Shortwave_incoming.Fsd_gapfill(variable_to_fill,myBaseforResults,mypathforAWAPdata,FileName_AWAPData,New_combined,Site_ID,Tower_Lat,Tower_Long,altitude,FluxFreq)
	variable_to_fill = 'Fsu'  
	New_combined= Shortwave_outgoing.Fsu_gapfill(New_combined, variable_to_fill,FluxFreq)	
    	#Run Ians code to gapfill Fld, Downscaling routine for estimating half-hourly radiation from MODIS  estimates
	#To Run the Flu gap filling code we need to have ingested the MODIS data first (above)
	variable_to_fill = 'Flu'    
	New_combined=Longwave_outgoing.Flu_gapfill(New_combined,variable_to_fill,FluxFreq)
	#Run Ians code to gapfill Fld
	variable_to_fill = 'Fld'    
	New_combined=Longwave_incoming.Fld_gapfill(variable_to_fill,myBaseforResults,mypathforAWAPdata,FileName_AWAPData,New_combined,Site_ID,FluxFreq)
	#Run small script to gapfill Fn
	variable_to_fill = 'Fn'    	
	New_combined=Net_radiation.Fn_gapfill(New_combined,variable_to_fill)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)    
    
    #Now Perform the ANN gap filling for Fc, Fh, Fe, Fg.  Call the routine. It is coded so that it can accept multiple inputs and outputs
    if cf['Options']['Gapfill_fluxes_ANN']=='Yes':      
    	#The basic model is an tmlgraph((number_of_inputs,24,16,number_of_outputs))  # Creates multilayer network full connectivity list
	# Training is done using TNC learning.  If that doesnt work well you can change the actual ANN model in the routine
	#Pass list of inputs and outputs as well as number of iterations.  About 500-1000 is OK.  The more the better but takes longer
	# Returns dataframe and columns with ANN output 'XXX_NN' constructed series 'XXX_Con and 'XXX_Con_QCFlag' = 1 if valid data from the tower else 99   
	# Load in previous data step data
	print "Starting ANN gapfilling" 
	frequency = cf['Options']['ANN_gapfill_freq']
	Use_Fc_Storage=cf['Options']['Use_Fc_Storage']
	New_combined= pd.read_pickle(DFfilename) 
	bestAWS_ID = pickle.load(open((myBaseforResults+'/AWS/bestAWS_ID_'+str(Site_ID)), 'rb'))
	#Set number of iterations for ANN
    	iterations=1000
	#Call as many times as needed for each variable or set of variables
	
	list_in=['Fsd_Con','VPD_Con','Sws_Con','Ts_Con',Ws_variable_name+'_Con','250m_16_days_EVI_new_interp']
	target=['Fc']
	#Call ANN routines
	New_combined=ANN.ANN_gapfill(myBaseforResults,New_combined,Site_ID,list_in,target,iterations,frequency,Use_Fc_Storage)
    
	list_in=['Fsd_Con','VPD_Con','Sws_Con','Ts_Con',Ws_variable_name+'_Con','250m_16_days_EVI_new_interp']
	target=['Fe','Fh','Fg']
	New_combined=ANN.ANN_gapfill(myBaseforResults,New_combined,Site_ID,list_in,target,iterations,frequency,Use_Fc_Storage)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)  
    
    #Calculate Barr et al. ustar threshold and uncertainty. It saves output to 'annual_statistics.csv'
    if cf['Options']['Ustar_Barr_calculate']=='Yes':   	
	print "Starting Barr et al.  Get some coffee!"
	Fluxfilepath=mypathforFluxdata+'/'+cf['Files']['FLUXfilename']
	CPD.CPD_main(Fluxfilepath, myBaseforResults, Ustar_Barr_bootstraps)
    
    #If user chooses 'ustar_Barr' or 'auto' as ustar methods
    #This is done seperately here.  If you have already done the calculations and produced the annual stats then you can 
    #selected NO to calculate but still choose ustar_Barr and it will use the stats already calculated.
    
    #OPen main DF and create a new column for Barr with Nans	
    New_combined= pd.read_pickle(DFfilename) 
    New_combined['ustar_Barr'] = np.nan
    
    #Try to open stats file.  In any vase write nans to the ustar Barr column
    try:
	annual_stats= pd.read_csv(myBaseforResults+"/CPD/Results/annual_statistics.csv",index_col=0)
	# Now take the results and apply that to the dataframe and add a new column of Barr ustar for each year.
	for i in annual_stats.index:
	    New_combined['ustar_Barr'].ix[str(i)]=annual_stats['ustar_mean'].ix[i]
	#Sometimes there is a years missing from the annual stats sheet.  Average ustar and gap fill any missing values 
	New_combined['ustar_Barr'].fillna(annual_stats['ustar_mean'].mean(), inplace=True)	
    except:
	pass
       
    #Save new data
    New_combined.to_pickle(DFfilename) 
    New_combined.to_csv(CSVfilename)

    if cf['Options']['Ustar_Reichstein_calculate']=='Yes':   
	#For the u*-filtering the data set (storage corrected flux is assumed) is split into 6 temperature classes of the same sample size 
	#(according to quantiles) and for each temperature class the set is split into 20 u*-classes. The threshold is defined as the u*-class 
	#where the night-time flux reaches more than 95% of the average flux at the higher u*-classes. The threshold is only accepted if for the 
	#temperature class if temperature and u* are not or only weakly correlated (|r| < 0.3). The final threshold is defined as the median of the 
	#thresholds of the (up-to) six temperature classed. This procedure is applied to the subsets of four three-months periods to account for 
	#seasonal variation of vegetation structure. For each period the u*-threshold is reported, but the whole data set is filtered according to 
	#the highest threshold found (conservative appraoch). In cases where no u*-threshold could be found it is set to 0.4 . A minimum threshold 
	#is set to 0.1. Each half-hourly value of NEE with the corresponding u*-below the threshold and each succeeding NEE value are removed.
	# Routine will return 3 month moving window ustar value 'ustar_threshold' and max ustar for period 'ustar_max'
	# Load in previous data step data
	New_combined=pd.read_pickle(DFfilename) 
	#Call Ustar routines
	New_combined=ustar.ustar_filter(New_combined,myBaseforResults,Site_ID)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  	
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)
	
    # Apply selected ustar threshold type
    # Here options for Ustar_filter_type are "auto" which uses default ustar threshold scheme which is currently Barr  et al.
    # Other options are a numerical threshold manually set (i.e. 0.17), 
    # Or set for Reichstein et al approach
    # Or use Reichstein with maximum value over entire timeseries which can be multiple years "ustar_Reichstein_Max", 
    # Or 1 month window maximum allows ustar threshold to vary "ustar_Reichstein_window"    
    # Define a column that states the actual ustar threshold used  for ANN and later GPP calculations, 
    # Defined here first based on Ustar_filter_type
    Ustar_filter_type=cf['Options']['Ustar_filter_type']    
    # Load in previous data step data
    New_combined=pd.read_pickle(DFfilename)    
    if Ustar_filter_type == "ustar_Reichstein_Max": 
	New_combined['ustar_used'] = New_combined['ustar_Reich_max']
    elif Ustar_filter_type == "ustar_Reichstein_window": 
	New_combined['ustar_used'] = New_combined['ustar_Reich_var']
    elif Ustar_filter_type == "ustar_Barr": 
	New_combined['ustar_used'] = New_combined['ustar_Barr']   
    elif Ustar_filter_type == "auto": 
	#Current definition for auto is Barr approach
	New_combined['ustar_used'] = New_combined['ustar_Barr']
    else:
	#Create a new column with the manual ustar value if this is set from Ustar_filter_type
	New_combined['ustar_manual'] = float(Ustar_filter_type)
	New_combined['ustar_used'] = float(Ustar_filter_type)    
    #Save new data
    New_combined.to_pickle(DFfilename) 	
	
    if cf['Options']['Calculate_Fre_ANN']=='Yes':      
	#Here calculate Fre (Re) using the ANN networks as above.  Code is the same except use nightime values
	# and ustar filtered data to traing to get Re=Fc
	print "Starting Fre calculation using ANN"
	print "Using the following ustar threshold " + str(Ustar_filter_type)
	frequency = cf['Options']['ANN_gapfill_freq']
	New_combined= pd.read_pickle(DFfilename) 
	#Set number of iterations for ANN
	iterations=300
	#set other options.  Use only evening data.  Can help if getting a bulge in early hours due to morning flush of CO2
	#Here chose whether to use use flux data from the first 3 hours after sunset where the canopy is still coupled with the atmosphere, as shown in Van Gorsel et al. (2007)
	evening = True
	#Set min and max threshold Fc in UMOL .m-2 s-1! 
	min_threshold = -1.0                    #UMOL .m-2 s-1
	max_threshold = 15.0                    #UMOL .m-2 s-1
	#Call as many times as needed for each variable or set of variables
	list_in=['Sws_Con','Ts_Con','Ta_Con','250m_16_days_EVI_new_interp']
	list_out=['Fc']
	#Call ANN routines # This time need Lat and long for solar calcs to get sunrise and set times
	New_combined=Fre_ANN.Fre_ANN_gapfill(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,Tower_Lat,Tower_Long,frequency,evening,min_threshold,max_threshold,Ustar_filter_type)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)  
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename) 
	
    if cf['Options']['Calculate_Fre_Lasslop']=='Yes':      
	#Here calculate Fre (Re) using the LASSLOP approach.  
	print "Starting Fre using Lasslop"
	New_combined= pd.read_pickle(DFfilename) 
	#Call GPP routine
	New_combined = Fre_Lasslop.Lasslop(New_combined) 
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)  	
	
	
    if cf['Options']['Calculate_GPP']=='Yes':      
	#Here calculate GPP
	print "Starting GPP calculation using Ustar and Fre from ANN"
	New_combined= pd.read_pickle(DFfilename) 
	#Call GPP routine
	#Calculate the range of GPP, Re and NEP from different range of ustar methods and sensitivityies (i.e. 0.5 ustar). Read in value from CF
	ustarrange=cf['Options']['Calculate_GPP_ustar_range']
	New_combined=GPP.GPP_calculate(myBaseforResults, New_combined, Site_ID, versionID, ustarrange)
	#Check dataframe for duplicates, set freq as required and pad as necessary and sort.  Call function
	New_combined = dataframe_check(New_combined, FluxFreq)	
	#Save new data
	New_combined.to_pickle(DFfilename) 
	New_combined.to_csv(CSVfilename)  
	
    if cf['Options']['Isaac_Fingerprint_plots']=='Yes':      
	#Here Plot Fingerprint plots as per the control file
	New_combined= pd.read_pickle(DFfilename) 
	CFname=cf['Options']['Fingerprint_file']
	Fingerprints.fingerprint_plots(myBaseforResults,New_combined,Site_ID,CFname,versionID)

    if cf['Options']['Diagnostics']=='Yes':      
	#Here do diagnostics from gap filling
	print "Starting Diagnostics"
	New_combined= pd.read_pickle(DFfilename) 
	# set of variable list to do diagnostics
	list_in=['Sws_Con','Ts_Con']
	do_results=True
	#Call ANN routines # This time need Lat and long for solar calcs to get sunrise and set times
	Diagnostics.basic_diags(myBaseforResults,New_combined,Site_ID,list_in, Ws_variable_name,do_results,Rain_label_variable_to_fill,versionID)

    if cf['Options']['Output_yearly_EddyProc_files']=='Yes':      
	#Here do diagnostics from gap filling
	New_combined= pd.read_pickle(DFfilename) 
	#Will produce yearly files for inpout to MPI EddyProc online tools in correct format, calc VPD, tab del, time stamp etc
	REddyProc.Output_files(New_combined,myBaseforResults,Site_ID,versionID)  

    if cf['Options']['Output_yearly_ARM_files']=='Yes':      
	#Here do diagnostics from gap filling
	New_combined= pd.read_pickle(DFfilename) 
	#Will produce yearly files for ARM
	ARMProc.Output_files(New_combined,myBaseforResults,Site_ID,versionID)  

    if cf['Options']['Output_yearly_WAVES_files']=='Yes':      
	#Here do diagnostics from gap filling
	New_combined= pd.read_pickle(DFfilename) 
	#Will produce yearly files for WAVES model
	WAVESProc.Output_files(New_combined,myBaseforResults,Site_ID,versionID)  
	
    if cf['Options']['Output_yearly_NASA_files']=='Yes':      
	#Here do diagnostics from gap filling
	New_combined= pd.read_pickle(DFfilename) 
	#Will produce yearly files for NASA SMAP model
	Ws_variable_name =  cf['Options']['Ws_variable_name']
	NASAProc.Output_files(New_combined,myBaseforResults,Site_ID,versionID,Ws_variable_name)      
	
    print "***************** Finished "+Site_ID+" Hooray ************************"
    
    
main()