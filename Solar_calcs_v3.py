############################################################################
# Coded by Jason Beringer
# Based on NOAA solar calculations spreadsheet
# Inputs = Date_input_excel = date input i.e. excel 1/1/2010 = 40179. SO year 1900 is =1
# Inputs = lat/long
# Run the Timezone_v3.py to get timezone offset (in hours) for your lat/long first.  Then pass to this function
# For calculating solar noon we only need the date (not time)
# However capacity to pas time as well and add calculations for diurnal variables like sun angle etc.
# In that case the time during the day would vary below.  At the moment it is invariant
#
# Programmed by Jason (July, 2013)
############################################################################

import datetime as dt
import time
import string
import math

def solar_calculations(Date_input_excel,latitude,longitude,timezone):
      
    Time_past_midnight_local=0.5
    
    #calculations
    Julian_day=Date_input_excel+2415018.5+Time_past_midnight_local-timezone/24
    Julian_century=(Julian_day-2451545)/36525
    Mean_obliqu_ecliptic=23+(26+((21.448-Julian_century*(46.815+Julian_century*(0.00059-Julian_century*0.001813))))/60)/60
    Obliq_corr=Mean_obliqu_ecliptic+0.00256*math.cos(math.radians(125.04-1934.136*Julian_century))
    Geom_mean_long_sun=(280.46646+Julian_century*(36000.76983 + Julian_century*0.0003032))%(360)
    Geom_mean_anomaly_sun=357.52911+Julian_century*(35999.05029 - 0.0001537*Julian_century)
    var_y=math.tan(math.radians(Obliq_corr/2))*math.tan(math.radians(Obliq_corr/2))
    Eccent_earth_orbit=0.016708634-Julian_century*(0.000042037+0.0000001267*Julian_century)
    eq_of_time=4*math.degrees(var_y*math.sin(2*math.radians(Geom_mean_long_sun))-2*Eccent_earth_orbit*math.sin(math.radians(Geom_mean_anomaly_sun))+4*Eccent_earth_orbit*var_y*math.sin(math.radians(Geom_mean_anomaly_sun))*math.cos(2*math.radians(Geom_mean_long_sun))-0.5*var_y*var_y*math.sin(4*math.radians(Geom_mean_long_sun))-1.25*Eccent_earth_orbit*Eccent_earth_orbit*math.sin(2*math.radians(Geom_mean_anomaly_sun)))
    Sun_eq_of_cntr=math.sin(math.radians(Geom_mean_anomaly_sun))*(1.914602-Julian_century*(0.004817+0.000014*Julian_century))+math.sin(math.radians(2*Geom_mean_anomaly_sun))*(0.019993-0.000101*Julian_century)+math.sin(math.radians(3*Geom_mean_anomaly_sun))*0.000289
    Sun_true_long=Geom_mean_long_sun+Sun_eq_of_cntr
    Sun_app_long=Sun_true_long-0.00569-0.00478*math.sin(math.radians(125.04-1934.136*Julian_century))
    sun_declination=math.degrees(math.asin(math.sin(math.radians(Obliq_corr))*math.sin(math.radians(Sun_app_long))))
    HA_sunrise=math.degrees(math.acos(math.cos(math.radians(90.833))/(math.cos(math.radians(latitude))*math.cos(math.radians(sun_declination)))-math.tan(math.radians(latitude))*math.tan(math.radians(sun_declination))))    
    
    #outputs
    solar_noon=(720-4*longitude-eq_of_time+timezone*60)/1440
    solar_sunrise=(solar_noon*1440-HA_sunrise*4)/1440
    solar_sunset=(solar_noon*1440+HA_sunrise*4)/1440
    true_solar_time=(Time_past_midnight_local*1440+eq_of_time+4*longitude-60*timezone)%(1440)
    
    #Return solar noon but could calculate sunset or setset and pass that back too.  There is capacity to add further calculations for zenith angle etc.
    #Returned time is id fraction of day (0-1)
    
    return (solar_sunrise,solar_noon,solar_sunset)

