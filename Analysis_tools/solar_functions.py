# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 15:59:38 2016

@author: imchugh
"""

import ephem
import datetime as dt
import numpy as np
import pdb

def get_ephem_solar(data_dict, lat, lon, alt, GMT_zone, return_var = 'zenith'):
    """
    Pass as positional args:
        1) 'data_dict': dictionary containing key / value pairs of:
                - 'date_time': numpy array of datetimes
                - 'T': numpy array of air temperature (float); optional but 
                       affects calculation of sunrise and sunset due to 
                       refraction effects
                - 'P': numpy array of air pressure (float); optional but 
                       affects calculation of sunrise and sunset due to 
                       refraction effects
        2) 'lat': latitude either in float (decimal radians) or str (decimal
                  degrees); this is just ephem convention!
        3) 'lon': longitude either in float (decimal radians) or str (decimal
                  degrees); this is just ephem convention!
        4) 'alt': altitude in m (int or float)
        5) 'GMT_zone': time zone relative to Greenwich (int or float)
    pass as kwarg:
        1) 'return_var': type of observation required (str); choices are:
                - 'zenith' (solar zenith for all datetimes in series)
                - 'altitude' (solar elevation for all datetimes in series)
                - 'rise' (sunrise time - only one per date)
                - 'set' (sunset time - only one per date)
    Returns:
        dictionary containing either ...
    """
    
    # Map function variable names to ephem sun properties
    var_dict = {'zenith': 'alt',
                'altitude': 'alt',
                'rise': 'rise_time',
                'set': 'set_time'}

    # Check the requested variable is available for output
    if not return_var in var_dict:
        raise Exception('Invalid input variable! Valid variables are zenith, '\
                        'altitude, rise and set... exiting')

    # Check if user wants sunrise or sunset times
    time_bool = True if return_var == 'rise' or return_var == 'set' else False

    # set local variable datetime
    date_time = data_dict['date_time']

    # Create and populate local observer (lat and long must be strings)
    obs = ephem.Observer()
    obs.elev = alt
    obs.lat = str(lat)
    obs.long = str(lon)

    # Specify body
    sun = ephem.Sun(obs)

    # Convert to array if single value
    try: 
        iter(date_time)
    except:
        date_time = np.array([date_time, ])

    # Convert to UTC
    UTC_datetime = date_time - dt.timedelta(hours = GMT_zone)

    # Convert to array if single value
    try: 
        iter(UTC_datetime)
    except:
        UTC_datetime = np.array([UTC_datetime, ])

#    pdb.set_trace()

    # Convert to dates if user wants sunrise or sunset times
    if time_bool:
        UTC_datetime = np.unique(np.array([this_date.date() for 
                                           this_date in UTC_datetime]))
        my_datetime = np.unique(np.array([this_date.date() for 
                                          this_date in date_time]))

    # Get data for each date_time
    var_list = []
    eval_string = 'sun.' + var_dict[return_var]
    for i, this_dt in enumerate(UTC_datetime):
        if 'T' in data_dict.keys():
            obs.temp = data_dict['T'][i]
        if 'P' in data_dict.keys():
            obs.pressure = data_dict['P'][i]
        obs.date = this_dt
        sun.compute(obs)
        this_val = eval(eval_string)
        var_list.append(this_val)    

    # Convert to python datetime.time if user wants sunrise or sunset times
    if time_bool: 
        for i, this_val in enumerate(var_list):
            (y, m, d, h, mins, s) = this_val.tuple()
            this_time = dt.time(h, mins, int(s))
            var_list[i] = (dt.datetime.combine(dt.date.today(), this_time) + 
                           dt.timedelta(hours = GMT_zone)).time()
        return {'date': my_datetime,
                return_var: np.array(var_list)}
    else:
        out_array = np.array(var_list)
        if return_var == 'zenith':
            out_array = np.pi / 2 - out_array
        return {return_var: out_array}
    
# Estimate clear sky radiation
def Insol_calc(data_dict, GMT_zone, latit, longit, ALT_m, k, use_ephem = False):
    """
    Pass args: 
        1) data_dict: must contain equal-length arrays, with a minimum of 
                      key / value pair of:
                          - 'date_time': array of datetimes (python datetime);
                      additionally, if using pyephem for zenith calculation and
                      accurate refraction correction is required, include two 
                      additional key / value pairs of: 
                          - 'T': array of temperatures in C (float)
                          - 'P': array of pressures in hPa (float)
        2) GMT_zone: time zone in decimal hours (int or float),
        3) latit: latitude in decimal degrees (int or float),
        4) longit: longitude in decimal degrees (int or float),
        5) alt: altitude in mASL (int or float); 
        6) k: extinction coefficient - unitless with range 0-1 (int or float)

    Optional kwargs:
        1) use_ephem - use pyephem instead of algorithms here (boolean); if 
                       valid entries for T and P are found in data_dict, the 
                       atmospheric refraction correction will be more accurate
                       (otherwise defaults to )
    
    The algorithms used are from the references below:    
   
    DiLaura, D. L. (1984), IES Calculation Procedures Committee Recommended
    practice for the calculation of daylight availability, J. Illuminating
    Engineering Soc. of North America, 13(4), 381-392.
    
    Duffie, J. and W. Beckman (1980). Solar Engineering of Thermal Processes. 
    New York, John Wiley and Sons.

    Kasten, F. and Young, A.T. (1989) Revised optical air mass tables and 
    approximation formula. Applied Optics 28: 4735-4738.
    
    Note that solar disk diameter or refraction corrections are not yet 
    included; if u want these, set use_ephem to True
    """

    # Get date and time components
    date_time = data_dict['date_time']    
    try: 
        iter(date_time)
    except:
        date_time = [date_time]
        
    DOY = np.array([i.timetuple().tm_yday for i in date_time])
    hour = np.array([i.hour for i in date_time])
    minute = np.array([i.minute + 15 for i in date_time])

    # Calculate equation of time correction, solar noon, declination and TOA radiation
    EqofTime = (0.17 * np.sin(4 * np.pi * (DOY-80) / 373) 
                - 0.129 * np.sin(2 * np.pi *(DOY-8) / 355)) # DiLaura (1984)
    solar_noon = 12 + (GMT_zone * 15.0 - longit) / 360 * 24 - EqofTime # Me
    decl = np.radians(23.4) * np.sin((DOY + 284) / 365.0 * 2 * np.pi) # Oke (1987)
    TOArad = (1 + 0.034 * np.cos(DOY / 365.25 * 2 * np.pi)) * 1367.0 # Duffie and Beckman (1980)

    # Calculate hour angle    
    hr_angle = abs(np.radians((solar_noon - (minute/60.0 + hour)) * 15))

    # Calculate solar zenith angle
    if use_ephem:
        date_time = np.array([this_dt - dt.timedelta(minutes = 15) 
                              for this_dt in date_time])
        zenith = get_ephem_solar(data_dict, latit, longit, ALT_m, GMT_zone)
    else:
        zenith = np.arccos(np.sin(np.radians(latit)) * np.sin(decl) + 
                 np.cos(np.radians(latit)) * np.cos(decl) * np.cos(hr_angle))
    zenith_msk = np.ma.masked_greater_equal(zenith, np.pi / 2) # Mask night values

    # Calculate optical air mass term (correct for optical effects of reduced 
    # air density [by calculating ratio of pressure at altitude to pressure at
    # sea level])
    # Kasten and Young (1989)
    m_sl = 1 / (np.cos(zenith_msk) + 0.50572 * 
                (96.07995 - np.degrees(zenith_msk)) ** -1.6364)
    p_p0 = np.exp(-0.0001184 * ALT_m) # Generic pressure/height approximation
    m = p_p0 * m_sl

    # Instantaneous clear sky surface radiation in Wm-2 (Beer-Lambert variant)
    Kdown = (TOArad * np.exp(-k * m) * np.cos (zenith_msk)).filled(0)
    m = m.filled(np.nan)

    # Make result dict    
    d = {}    
    d['solar_noon'] = solar_noon
    d['declination'] = decl
    d['TOA_radiation'] = TOArad
    d['zenith'] = zenith
    d['optical_air_mass'] = m
    d['Kdown'] = Kdown
    
    return d