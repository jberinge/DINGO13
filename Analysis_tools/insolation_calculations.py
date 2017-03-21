# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 16:32:06 2014

@author: imchugh
"""
import numpy as np
import pdb

def Insol(DOY,k,d):
    
    # For each day calculate equation of time correction, solar noon, declination and TOA radiation
    array_EqofTime=0.17*np.sin(4*np.pi*(DOY-80)/373)-0.129*np.sin(2*np.pi*(DOY-8)/355) # DiLaura (1984)
    array_solar_noon=12+(d['GMT_zone']*15.0-d['long_decdeg'])/360*24-array_EqofTime # Me
    array_decl=np.radians(23.4)*np.sin((DOY+284)/365.0*2*np.pi) # Oke (1987)
    array_TOArad=(1+0.034*np.cos(DOY/365.25*2*np.pi))*1367.0 # Duffie and Beckman (1980)
    
    # Create an hour angle array for each minute of day and each day of year
    array_h=np.tile(np.linspace(0,1439.0/1440*24,num=1440),(len(DOY),1))
    array_h=abs(np.radians((array_solar_noon.reshape(len(DOY),1)-array_h)*15))
    
    # Duplicate declination array for each time of day
    array_decl=np.tile(array_decl,(1440,1)).T

    # Calculate zenith angles
    array_z=np.arccos(np.sin(np.radians(d['lat_decdeg']))*np.sin(array_decl)+
            np.cos(np.radians(d['lat_decdeg']))*np.cos(array_decl)*np.cos(array_h))
    array_z_msk=np.ma.masked_greater_equal(array_z,np.pi/2) # Mask night values    
  
    # Calculate optical air mass term for all valid Z 
    array_m=(np.exp(-1*d['ALT_m']/8343.5)/(np.cos(array_z_msk)+0.15*
            (np.degrees(90-array_z_msk)+3.855)**-1.253)) # Wunderlich (1972)           
    
    # Instantaneous clear sky surface radiation in Wm-2 for each minute of the day
    array_Kdown_clr_mins=array_TOArad.reshape(len(array_TOArad),1)*np.exp(-k*array_m)*np.cos(array_z_msk)
    
    # Aggregate one-minute instantaneous clear sky rad to period average
    array_Kdown_clr_hr=np.empty([len(DOY),1440/d['rec_length']])
    for i in xrange(len(DOY)):
        array_temp=array_Kdown_clr_mins[i][:].reshape(1440/d['rec_length'],d['rec_length']) # Temporary bins
        array_Kdown_clr_hr[i][:]=np.ma.mean(array_temp,axis=1) # Average bin content  
    
    return array_Kdown_clr_hr