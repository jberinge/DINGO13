# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 11:16:19 2014

@author: imchugh
"""
import pdb
import numpy as np

# Humidity functions and conversions
def es(t=20):
    """Pass temperature (t) in celsius; returns saturation vapour pressure in kPa"""
    return 0.611*10**(7.5*t/(237.3+t))
    
def VPD(e,t):
    """Pass temperature (t) in celsius and vapour pressure (e) in kPa;
       returns vapour pressure deficit in kPa"""    
    return 0.611*10**(7.5*t/(237.3+t)) - e
    
def q_to_e(q, p=101.3):
    """Pass specific humidity (q) in g.kg-1, barometric pressure (p) in kpa;
       returns vapour pressure in kPa"""
    return (q*28.97)/28.97*p
    
def Ah_to_e(Ah, Ta):
    """Pass absolute humidity (Ah) in g.m-3 and air temperature in celsius;
       returns vapour pressure in kPa"""
    return Ah / 18 * (Ta + 273.15) * 8.3143 / 10**3

def Lv(t):
    """Pass temperature in C; returns Lv in J g-1;
    note: valid for range -25 to +40C"""
   
    return 2500.8 -2.36 * t + 0.0016 * t**2 - 0.00006 * t**3  