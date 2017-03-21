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
import struct
import math
import os

from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
import reportlab.pdfgen 
from reportlab.lib.units import inch

def dataframe_check(Dataframe, FluxFreq):
    #Check dataframe for duplicates, pad as necessary and sort
    Dataframe.sort(inplace=True)
    Dataframe["index"] = Dataframe.index
    Dataframe.drop_duplicates('index', take_last=True, inplace=True)
    del Dataframe["index"]	
    Dataframe=Dataframe.asfreq(FluxFreq, method=None)      
    return Dataframe

def excel_to_pydate(exceldate):
    datemode=0           # datemode: 0 for 1900-based, 1 for 1904-based
    pyear, pmonth, pday, phour, pminute, psecond = xlrd.xldate_as_tuple(exceldate, datemode)
    py_date = dt.datetime(pyear, pmonth, pday, phour, pminute, psecond)
    return(py_date)

def distance_on_unit_sphere(lat1, long1, lat2, long2):

    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
        
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
        
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
        
    # Compute spherical distance from spherical coordinates.
        
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
    
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )*6373

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc


def get_BoM_stations(BoMfilename):
    #File is a list of AWS stations from the Bureau of Meteorology
    #ftp://ftp.bom.gov.au/anon2/home/ncc/metadata/sitelists/stations.zip
    #the file has then had the header and footer removed for ease of import
    #Open file and read lines from file in one statement
    with open(BoMfilename, 'r') as stations_file:
        lines = stations_file.readlines()    
    #The Vraibale headings for this file are "SiteID  Dist  Sitename Start End Lat Lon Source  STA Height (m)   Bar_ht    WMO"
    #Define the field widths of each variable in the file
    fieldwidths = (8, 6, 44,8,5, 10,9,15,4,11,9,6)
    fmtstring = ''.join('%ds' % f for f in fieldwidths)
    parse = struct.Struct(fmtstring).unpack_from

    
    #Go through each line in the file and parse it to get the variables then append to new list
    stations = [] 
    for line in lines:
        p = parse(line)
        stations.append(p)
    return stations

def dofilestuff(myBaseforResults):
    #Check for place to put resuolts - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/AWS"):
        os.mkdir(myBaseforResults+"/AWS")
    mypathforResults=myBaseforResults+"/AWS"
    return mypathforResults

def get_AWS_data(BoMfilename,Site_ID,Tower_Lat,Tower_Long,mypathforAWSdata,myBaseforResults, OffsetAWS30min,FLUXDataframe,FluxFreq):
    #Set the start year for getting AWS data otherwise it goes back a LONG way and takes too long.  Injest data required only
    Flux_start_year=FLUXDataframe.index[0].year
    print "Flux start year is ",Flux_start_year
    #Set paths and directories
    mypathforResults=dofilestuff(myBaseforResults)
    print "Getting AWS files and get station ID that we have from Dropbox data folder"
    templist=os.listdir(mypathforAWSdata)
    AWS_files=[]
    #Go through the temp list and take only the first 6 characters which is the station ID number
    #This is the 30 minute AWS data we have on Dropbox
    for filen in templist:
        stationID=filen[:6]
        AWS_files.append(stationID)   
    
    #Call function to return the parsed BoM text file
    print "Call function to return the parsed BoM text file"
    BoMstations = get_BoM_stations(BoMfilename)

    #Define new list
    #Convert parsed file into tuple to be passed to table maker
    print "Convert parsed file into tuple to be passed to table maker"
    distances=[]
    for station in BoMstations:
        # For each line look for values that are actually ".." in the file.  These can not be converted to float
        #Convert values to floats or integers along the way        
        BoM_Lat=float(station[5])
        BoM_Long=float(station[6])
        BoM_ID=str.strip(str(station[0]))
        BoM_name=(station[2])
        #Now some variables have ".." in the file which can not be processed. Check first then assign to blank if ".."

        if str.strip(station[3])=='..':
            BoM_start=''
        else:
            BoM_start=int(station[3])
        #If the end date is ".." thenit is still current and assign string "present" which is used later 
        #This is because we only want used AWS data when we HAVE it and when it is still present
        
        if str.strip(station[4])=='..':
            BoM_end='present' 
        else:
            BoM_end=int(station[4])
        if str.strip(station[9])=='..':
            BoM_elev='' 
        else:
            BoM_elev=float(station[9])
        if str.strip(station[0]) in AWS_files:
            if BoM_start<=Flux_start_year:
                haveAWS=True
            else:
                haveAWS=False
        else:
            haveAWS=False
        
        #Here some AWS stations dont have all the data (like P) and this causes program to crash
        #So for now have an excludelist for station ID.  Should program this later to check for valid
        #data from each column of each station.
        excludeAWS=['081049','080091','083085','088164','088051','072161','086371','066212','067119','067117','066059','066194','066062',
                    '070349','088109','086351','004021','097024','094191','094220','086104','040284','014046','014072','014277','014031',
                    '010917','009240','067113','086395','086396','086397','086398','088001','009214','009265','082139','074037','095003','094087','097083','083024','085296']

        if BoM_ID in excludeAWS:
                haveAWS=False
                
        #Call function to calculate distance      
        dist=distance_on_unit_sphere(Tower_Lat, Tower_Long, BoM_Lat, BoM_Long)
        BoMtuple=(BoM_ID,dist,haveAWS,BoM_Lat,BoM_Long,BoM_name,BoM_start,BoM_end,BoM_elev)
        distances.append(BoMtuple)    
    
    #Sort the list by distance
    print "Sort list of AWS files from dropbox"   
    sorteddist=sorted(distances, key=lambda ID: ID[1])
    #Filter the list by data when we HAVE it and when it is still present
    sorteddistfile = open('sorteddistfile.txt', 'w')    

    print "Filter list of AWS files from dropbox" 
    #Get the top sites that we have data for AND are current
    #And where the Bom Data coveres the start year of the flux data
    AWSsorted = [t for t in sorteddist if t[2]==True and t[7]=='present']

    #Output a list with top 3 ID and names
    bestAWS_ID=[]    

    for item in AWSsorted[:3]:
        bestAWS_ID.append(item[0])
    print " These are the top 3 matches in closest distance", bestAWS_ID


    #Outputs (Tables and PDF)
    #===============================
       
    #sorted list by distance of all BoM stations regardless of whether we have the data or not
    #Define the number of rows to output
    print "Create AWS location Table in PDF output" 
    rowsto_output1=33
    doc1name=mypathforResults+'/'+"Sorted list AWS stations by dist to flux tower_"+Site_ID+".pdf"
    doc1 = SimpleDocTemplate(doc1name, pagesize=A4)
    # container for the 'Flowable' objects
    elements1 = []
    t1=Table(sorteddist[0:rowsto_output1])
    t1.setStyle(TableStyle([('BACKGROUND',(0,0),(0,rowsto_output1),colors.green),('BACKGROUND',(1,0),(1,rowsto_output1),colors.yellow)]))
    #set_column_titles(array('ID','Dist', 'Have Data','Lat','Long','Name','Start', 'End','Elev'))    
    elements1.append(t1)
    # write the document to disk
    doc1.build(elements1)    
    
    #Define the number of rows to output
    rowsto_output2=10    
    # container for the 'Flowable' objects
    elements2 = []
    #Here output top ten in sorted and filtered list
    t2=Table(AWSsorted[0:rowsto_output2])
    t2.setStyle(TableStyle([('BACKGROUND',(0,0),(0,rowsto_output2),colors.green),('BACKGROUND',(1,0),(1,rowsto_output2),colors.yellow)]))
    elements2.append(t2)
    
    #Output PDF Canvas
    doc3name=mypathforResults+'/'+"Top 10 AWS sorted for flux tower_"+Site_ID+".pdf"    
    doc3 = SimpleDocTemplate(doc3name,pagesize=A4,
                        rightMargin=40,leftMargin=40,
                        topMargin=20,bottomMargin=20)
    width, height = A4
    container=[]
    ptext = ' <font size=12>The data shows sites that are still current and that we have 30 minute data for\
        The table shows the station ID, Distance (km), Do we have data, Lat, Long, Name, Start, End, Elev</font>'
    styles=getSampleStyleSheet()    
    container.append(Paragraph('<font size=14><b>TABLE comparison Flux Tower Site with BoM AWSstations</b></font>', styles["Normal"]))
    container.append(Spacer(1, 10))
    container.append(Paragraph(ptext, styles["Normal"]))
    container.append(Spacer(1, 16))
    container.append(t2)
    doc3.build(container)      


    #Now get the AWS files and import to Pandas
    #==========================================
    aws_df1_filename=mypathforAWSdata+ chr(92)+str(bestAWS_ID[0])+'.txt'
    aws_df2_filename=mypathforAWSdata+ chr(92)+str(bestAWS_ID[1])+'.txt'
    aws_df3_filename=mypathforAWSdata+ chr(92)+str(bestAWS_ID[2])+'.txt'
    
    #Here create a little function that gets passed the columns from the df and then combines into one datetime
    parse = lambda x: dt.datetime.strptime(x, '%Y %m %d %H %M')
    #Define the heading names
    headings=['hm','ID','Name','Rain','Rain_QC','Ta','Ta_QC','Tw','Tw_QC','DP','DP_QC','RH','RH_QC','WS','WS_QC','WD','WD_QC','GUST','GUST_QC','P','P_QC','AWSflag','Misc']
    #Then make custom headings that include ID for each of the three AWS files
    headings1=[]
    headings2=[]
    headings3=[]
    for item in headings:
        tempitem =item+'_'+str(bestAWS_ID[0])
        headings1.append(tempitem)
        tempitem =item+'_'+str(bestAWS_ID[1])
        headings2.append(tempitem)
        tempitem =item+'_'+str(bestAWS_ID[2])
        headings3.append(tempitem)      
        
    print "Parsing AWS files for dates"      
    aws_df1a = pd.read_csv(aws_df1_filename, parse_dates = [[3,4,5,6,7]], date_parser=parse, index_col=0,na_values="")
    aws_df2a = pd.read_csv(aws_df2_filename, parse_dates = [[3,4,5,6,7]], date_parser=parse, index_col=0,na_values="")
    aws_df3a = pd.read_csv(aws_df3_filename, parse_dates = [[3,4,5,6,7]], date_parser=parse, index_col=0,na_values="")
    
    #Apply the headings to the columns
    aws_df1a.columns=headings1
    aws_df2a.columns=headings2
    aws_df3a.columns=headings3
    
    ##Resample and set freq to 30 M.    
    #print "Resampling AWS files to Flux data frequency"    
    #aws_df1b=aws_df1a.asfreq(FluxFreq)
    #aws_df2b=aws_df2a.asfreq(FluxFreq)
    #aws_df3b=aws_df3a.asfreq(FluxFreq)
    
    ##Check for duplicates
    #print "Checking for duplicates in AWS files"
    #aws_df1c=aws_df1b.groupby(aws_df1b.index).first()
    #aws_df2c=aws_df2b.groupby(aws_df2b.index).first()
    #aws_df3c=aws_df3b.groupby(aws_df3b.index).first()   
    
    #Check dataframe for duplicates, pad as necessary and sort
    aws_df1c = dataframe_check(aws_df1a, FluxFreq)
    aws_df2c = dataframe_check(aws_df2a, FluxFreq)
    aws_df3c = dataframe_check(aws_df3a, FluxFreq)
        
    print "Station ID     Number items changed in DF check"
    print "=================================================================================="
    print str(bestAWS_ID[0]) + "                        "+ str((len(aws_df1c)-len(aws_df1a)))
    print str(bestAWS_ID[1]) + "                        "+ str((len(aws_df2c)-len(aws_df2a)))
    print str(bestAWS_ID[2]) + "                        "+ str((len(aws_df3c)-len(aws_df3a)))
            
    #Combine into single file, Use union of keys from both frames (outer)
    print "Combining AWS files"
    AWS_combined1=aws_df1c.join(aws_df2c,how='outer')
    AWS_combined=AWS_combined1.join(aws_df3c,how='outer')

    #Check dataframe for duplicates, pad as necessary and sort
    AWS_combined = dataframe_check(AWS_combined, FluxFreq)

    #Ossfet by 30 minutes if variablke is True
    if OffsetAWS30min==True:
        AWS_combined.shift(-1)
        
    print "Finished." 
    return(AWS_combined, bestAWS_ID)
   
    
