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
#   The index name is expected to be found on the first row.
#
# Programmed by Jason (Dec 1, 2012)
############################################################################

import datetime as dt
import numpy as np
import time
import urllib2
import string
import xml.etree.ElementTree as ET
import math



def get_timezone_info(latitude,longitude,currentdate):
    #get time zone info for site using http://api.askgeo.com/ Web API
    #Jason has account.  Need to enter those details in the http query string

    #Query the timezone and Astronamy databases.  Will return somethhing like
    # <TimeZone InDstNow="true" WindowsStandardName="Eastern Standard Time" CurrentOffsetMs="-14400000" 
    # ShortName="EDT" TimeZoneId="America/New_York" MinDistanceKm="0.0" AskGeoId="12571" IsInside="true"/>
    
    # Each database query returns a large number of variables.  We just print a few and return ONLY timezone
    # But you could return other variables if that was useful. http://askgeo.com/
    
	#Note in v3a we send the query the date and time during winter (2013/07/01) so that the timezone returned is not daylight saving
    url="http://api.askgeo.com/v1/841/9d6d6265d142bb57d292df6b2c1c7a7c1a4b78f2dd3734f475f6abf698bd0c52/query.xml?points="+str(latitude)+"%2C"+str(longitude)+"&databases=Point%2CTimeZone%2CAstronomy&dateTime="+str(currentdate)
    #Attempt to get web data.  Give three goes before giving an error.
    try:
        result = urllib2.urlopen(url)
        #Parse the returned XML
        tree = ET.parse(result)
        root = tree.getroot()
        #Look for the following elements
        a= tree.find('.//Astronomy')
        b= tree.find(".//TimeZone")
    except urllib2.URLError, e:
        try:
	    print " Failed 1st time trying web API trying again in 30 seconds"
	    time.sleep(100)
	    result = urllib2.urlopen(url)
	    #Parse the returned XML
	    tree = ET.parse(result)
	    root = tree.getroot()
	    #Look for the following elements
	    a= tree.find('.//Astronomy')
	    b= tree.find(".//TimeZone")       
	except urllib2.URLError, e:
	    try:
		print " Failed 2nd time trying web API trying again in 30 seconds"
		time.sleep(500)
		result = urllib2.urlopen(url)
		#Parse the returned XML
		tree = ET.parse(result)
		root = tree.getroot()
		#Look for the following elements
		a= tree.find('.//Astronomy')
		b= tree.find(".//TimeZone")
	    except urllib2.URLError, e:
		try:
		    print " Failed 3rd time trying web API trying again in 30 seconds"
		    time.sleep(1000)
		    result = urllib2.urlopen(url)
		    #Parse the returned XML
		    tree = ET.parse(result)
		    root = tree.getroot()
		    #Look for the following elements
		    a= tree.find('.//Astronomy')
		    b= tree.find(".//TimeZone")   		
		except urllib2.URLError, e:
		    print "Couldnt get the data from the web API. Tried 4 times"
		    handleError(e)
    
    #Start formatting the returned element to put into dictionary  
    tempTimeZonelist = ET.tostring(b)
    TimeZonelist = tempTimeZonelist.replace('\"','')
    TimeZonelist = tempTimeZonelist[10:-3]
    TimeZonelist=TimeZonelist.split()
    
    #Create a dictionary
    TimeZonedic = {}
    for entry in TimeZonelist[1:8]:
        key, val = entry.split('=')
        TimeZonedic[key] = val
    
    #Start formatting the returned element to put into dictioonary  
    tempAstronomylist = ET.tostring(a)
    Astronomylist = tempAstronomylist.replace('\"','')
    Astronomylist = tempAstronomylist[10:-3]
    Astronomylist = Astronomylist.split()
                      
    Astronomydic = {}
    for entry in Astronomylist:
        key, val = entry.split('=')
        Astronomydic[key] = val
	
    #Convert timezone to hours
    timezone_str= TimeZonedic.get('CurrentOffsetMs')
    timezone_str = timezone_str.replace('\"','')
    timezone=(float(timezone_str))/1000/60/60
    
    #get the variable we need from the list
    #PRINT THE VARIABLES WE WANT.  At the moment just a selection
    print "------------Data from Timezone Database"
    #print "CurrentOffsetMs of Time Zone",TimeZonedic.get('CurrentOffsetMs')
    #print "TimeZoneId",TimeZonedic.get('TimeZoneId')
    print "WindowsStandardName ",TimeZonedic.get('WindowsStandardName') + ",  " + str(timezone) + "   (hrs)"
    #print "ShortName",TimeZonedic.get('ShortName')
    #print "Is currently daylight savings ",TimeZonedic.get('InDstNow')
    #print "------------Data from Astronomy Database"
    #print "CurrentDateTimeIso8601", Astronomydic.get('CurrentDateTimeIso8601')
    print "TodaySolarNoonUtcMsecs", Astronomydic.get('TodaySolarNoonUtcMsecs')    
    print ""
    
    #return only time zone but can pass anything back
    return (timezone,TimeZonedic.get('InDstNow'))
