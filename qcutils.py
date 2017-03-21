# the following line needed for unicode character in convert_anglestring
# -*- coding: latin-1 -*-
import ast
import constants as c
import datetime
import dateutil
import logging
import math
import meteorologicalfunctions as mf
import netCDF4
import numpy
import os
import platform
import pytz
import sys
import time
import Tkinter,tkSimpleDialog
import xlrd
import xlwt

log = logging.getLogger('qc.utils')

def bp(fx,tao):
    """
    Function to calculate the b and p coeficients of the Massman frequency correction.
    """
    bp = 2 * c.Pi * fx * tao
    return bp

def cfkeycheck(cf,Base='Variables',ThisOne=[],key=[]):
    if len(ThisOne) == 0:
        return
    if len(key) == 0:
        if Base in cf.keys() and ThisOne in cf[Base].keys():
            return ThisOne in cf[Base].keys()
        else:
            return
    else:
        if Base in cf.keys() and ThisOne in cf[Base].keys():
            return key in cf[Base][ThisOne].keys()
        else:
            return

def cfoptionskeylogical(cf,Key='',default=False):
    if 'Options' in cf:
        if Key in cf['Options']:
            returnValue = cf.get('Options').as_bool(Key)
            #if str(cf['Options'][Key]).lower()=="true" or str(cf['Options'][Key]).lower()=="yes":
                #returnValue = True
            #else:
                #returnValue = False
        else:
            returnValue = default
    else:
        returnValue = default
    return returnValue

#def CheckQCFlags(ds):
    #"""
    #Purpose:
     #Make sure that all values of -9999 in a data series have a non-zero QC flag value.
    #Usage:
     #qcutils.CheckQCFlags(ds)
    #Author: PRI
    #Date: August 2014
    #"""
    #for ThisOne in ds.series.keys():
        #data = numpy.ma.masked_values(ds.series[ThisOne]["Data"],-9999)
        #flag = numpy.ma.masked_equal(ds.series[ThisOne]["Flag"],0)
        #mask = data.mask&flag.mask
        #index = numpy.ma.where(mask==True)[0]
        #ds.series[ThisOne]["Flag"][index] = numpy.int32(8)

def CheckQCFlags(ds):
    """
    Purpose:
     Make sure that all values of -9999 in a data series have a non-zero QC flag value.
    Usage:
     qcutils.CheckQCFlags(ds)
    Author: PRI
    Date: August 2014
    """
    # force any values of -9999 with QC flags of 0 to have a QC flag of 8
    for ThisOne in ds.series.keys():
        data = numpy.ma.masked_values(ds.series[ThisOne]["Data"],-9999)
        flag = numpy.ma.masked_equal(numpy.mod(ds.series[ThisOne]["Flag"],10),0)
        mask = data.mask&flag.mask
        index = numpy.ma.where(mask==True)[0]
        ds.series[ThisOne]["Flag"][index] = numpy.int32(8)
    # force all values != -9999 to have QC flag = 0, 10, 20 etc
    for ThisOne in ds.series.keys():
        index = numpy.where((abs(ds.series[ThisOne]['Data']-numpy.float64(c.missing_value))>c.eps)&
                            (numpy.mod(ds.series[ThisOne]["Flag"],10)!=0))
        ds.series[ThisOne]["Flag"][index] = numpy.int32(0)

def CheckTimeStep(ds):
    """
    Purpose:
     Checks the datetime series in the data structure ds to see if there are
     any missing time stamps.
     This function returns a logical variable that is true if any gaps exist
     in the time stamp.
    Useage:
     has_gaps = CheckTimeSTep(ds)
     if has_gaps:
         <do something about missing time stamps>
    Author: PRI
    Date: April 2013
    """
    # set the has_gaps logical
    has_gaps = False
    # get the number of records
    nRecs = int(ds.globalattributes["nc_nrecs"])
    # get the time step
    ts = int(ds.globalattributes["time_step"])
    # time step between records in seconds
    dt = get_timestep(ds)
    # indices of elements where time step not equal to default
    index = numpy.where(dt!=ts*60)[0]
    # check to see if ww have any time step problems
    if len(index)!=0:
        has_gaps = True
        log.warning(" CheckTimeStep: "+str(len(index))+" problems found with the time stamp")
    return has_gaps

def contiguous_regions(condition):
    """
    Purpose:
     Finds contiguous True regions of the boolean array "condition". Returns
     a 2D array where the first column is the start index of the region and the
     second column is the end index.
    Author: Joe Kington (via StackOverflow)
    Date: September 2014
    """
    # Find the indicies of changes in "condition"
    d = numpy.diff(condition)
    idx, = d.nonzero() 
    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = numpy.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = numpy.r_[idx, condition.size] # Edit
    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

def ConvertCO2Units(cf,ds,Cc='Cc'):
    Cc_units_out = "mg/m3"            # default value
    Cc_units_in = ds.series[Cc]['Attr']['units']
    if 'Options' in cf:
        if 'CO2Units' in cf['Options']:
            Cc_units_out = str(cf['Options']['CO2Units'])
    if Cc_units_out!=Cc_units_in:
        log.info(' Converting CO2 concentration from '+Cc_units_in+' to '+Cc_units_out)
        if Cc_units_out=="umol/mol" and Cc_units_in=="mg/m3":
            c_mgpm3,flag,attr = GetSeriesasMA(ds,Cc)
            T,f,a = GetSeriesasMA(ds,'Ta')
            p,f,a = GetSeriesasMA(ds,'ps')
            c_ppm = mf.co2_ppmfrommgpm3(c_mgpm3,T,p)
            attr["long_name"] = attr["long_name"]+", converted to umol/mol"
            attr["units"] = Cc_units_out
            attr["standard_name"] = "mole_concentration_of_carbon_dioxide_in_air"
            CreateSeries(ds,Cc,c_ppm,Flag=flag,Attr=attr)
        elif Cc_units_out=="mg/m3" and Cc_units_in=="umol/mol":
            c_ppm,flag,attr = GetSeriesasMA(ds,Cc)
            T,f,a = GetSeriesasMA(ds,'Ta')
            p,f,a = GetSeriesasMA(ds,'ps')
            c_mgpm3 = mf.co2_mgpm3fromppm(c_ppm,T,p)
            attr["long_name"] = attr["long_name"]+", converted to mg/m3"
            attr["units"] = Cc_units_out
            attr["standard_name"] = "mass_concentration_of_carbon_dioxide_in_air"
            CreateSeries(ds,Cc,c_mgpm3,Flag=flag,Attr=attr)
        else:
            log.info('  ConvertCO2Units: input or output units for CO2 concentration not recognised')

def ConvertFcUnits(cf,ds,Fc='Fc',Fc_storage='Fc_storage'):
    if 'Options' not in cf: return
    if 'FcUnits' not in cf['Options']: return
    # the user may want to change the units of Fc and Fc_storage
    Fc_units_out = str(cf['Options']['FcUnits'])
    # convert units of Fc if required
    if Fc in ds.series.keys():
        Fc_units_in = ds.series[Fc]['Attr']['units']
        if Fc_units_out!=Fc_units_in:
            log.info(' Converting CO2 flux from '+Fc_units_in+' to '+Fc_units_out)
            if Fc_units_out=="umol/m2/s" and Fc_units_in=="mg/m2/s":
                Fc_mgpm2ps,flag,attr = GetSeriesasMA(ds,Fc)
                Fc_umolpm2ps = mf.Fc_umolpm2psfrommgpm2ps(Fc_mgpm2ps)
                attr["long_name"] = attr["long_name"]+", converted to umol/m2/s"
                attr["units"] = Fc_units_out
                attr["standard_name"] = "surface_upward_mole_flux_of_carbon_dioxide"
                CreateSeries(ds,Fc,Fc_umolpm2ps,Flag=flag,Attr=attr)
            elif Fc_units_out=="mg/m2/s" and Fc_units_in=="umol/m2/s":
                Fc_umolpm2ps,f,a = GetSeriesasMA(ds,Fc)
                Fc_mgpm2ps = mf.Fc_mgpm2psfromumolpm2ps(Fc_umolpm2ps)
                attr["long_name"] = attr["long_name"]+', converted to mg/m2/s'
                attr["units"] = Fc_units_out
                attr["standard_name"] = "not defined"
                CreateSeries(ds,Fc,Fc_mgpm2ps,Flag=flag,Attr=attr)
            else:
                log.info('  ConvertFcUnits: input or output units for Fc unrecognised')
    # convert units of Fc_storage if required, just go with boiler plate for now
    if Fc_storage in ds.series.keys():
        Fc_storage_units_in = ds.series[Fc_storage]['Attr']['units']
        if Fc_units_out!=Fc_storage_units_in:
            log.info(' Converting CO2 storage flux from '+Fc_storage_units_in+' to '+Fc_units_out)
            if Fc_units_out=="umol/m2/s" and Fc_storage_units_in=="mg/m2/s":
                Fc_storage_mgpm2ps,flag,attr = GetSeriesasMA(ds,Fc_storage)
                Fc_storage_umolpm2ps = mf.Fc_umolpm2psfrommgpm2ps(Fc_storage_mgpm2ps)
                attr["long_name"] = attr["long_name"]+", converted to umol/m2/s"
                attr["units"] = Fc_units_out
                CreateSeries(ds,Fc_storage,Fc_storage_umolpm2ps,Flag=flag,Attr=attr)
            elif Fc_units_out=="mg/m2/s" and Fc_storage_units_in=="umol/m2/s":
                Fc_storage_umolpm2ps,f,a = GetSeriesasMA(ds,Fc_storage)
                Fc_storage_mgpm2ps = mf.Fc_mgpm2psfromumolpm2ps(Fc_storage_umolpm2ps)
                attr["long_name"] = attr["long_name"]+", converted to mg/m2/s"
                attr["units"] = Fc_units_out
                CreateSeries(ds,Fc_storage,Fc_storage_mgpm2ps,Flag=flag,Attr=attr)
            else:
                log.info('  ConvertFcUnits: input or output units for Fc_storage unrecognised')

def convertunits(old_data,old_units,new_units,ts,mode="quiet"):
    """
    Purpose:
     Generic routine for changing units.
     Nothing is done if the original units are the same as the requested units.
    Usage:
     new_data = qcutils.convertunits(old_data,old_units,new_units)
     where old_data is a 1D array of data in the original units
           old_units are the units of the original data
           new_units are the units of the new data
           ts is the time step
    Author: PRI
    Date: July 2015
    """
    if old_units==new_units: return
    # check the units are something we understand
    # add more lists here to cope with water etc
    co2_list = ["umol/m2/s","gC/m2"]
    ok_list = co2_list
    # parse the original units
    if old_units=="umol/m^2/s": old_units="umol/m2/s"
    if old_units.replace(" ","")=="umolm-2s-1": old_units="umol/m2/s"
    if old_units not in ok_list:
        msg = " Unrecognised units in quantity provided ("+old_units+")"
        log.error(msg)
        new_data = old_data
    elif new_units not in ok_list:
        msg = " Unrecognised units requested ("+new_units+")"
        log.error(msg)
        new_data = old_data
    elif old_units=="umol/m2/s" and new_units=="gC/m2":
        # do the units change
        new_data = old_data*12.01*ts*60/1E6
    elif old_units=="gC/m2" and new_units=="umol/m2/s":
        # do the units change
        new_data = old_data*1E6/(12.01*ts*60)
    else:
        msg = "Unrecognised units combination "+old_units+" and "+new_units
        log.error(msg)
    return new_data

def convert_anglestring(anglestring):
    """
    Purpose:
     Attempt to convert an angle string to a float.
    Usage:
     a = qcutils.convert_anglestring(astr)
     Acceptable input formats:
      astr = '''34 12' 24" S'''
      astr = '''34 12 24S'''
      astr = '''34 12'24.123"S''
      astr = '''34.123 S'''
      astr = '''-34.123'''
    """
    quadlist=["N","E","S","W"]
    direction = {'N':1, 'S':-1, 'E': 1, 'W':-1}
    try:
        # simple casting may work, who knows?
        return float(anglestring)
    except ValueError:
        # replace the degrees, minutes and seconds symbols with spaces
        new = anglestring.replace(u'\B0',' ').replace('\'',' ').replace('"',' ')
        # check there is a space between the quadrant letter (assumed to be one of N, E, W or S)
        # and the next character to the left
        # find out which of N, E, S, or W is in the string
        for item in quadlist:
            if item in new: quadletter=item
        # now get the index of this character in the string
        i=new.index(quadletter)
        # check that the next character to the left is a space character
        if new[i-1] != " ": new = new[0:i]+" "+new[i:]
        # now split the string on space characters
        new = new.split()
        # get the quadrant letter
        new_dir = new.pop()
        # make sure we have 3 parts
        new.extend([0,0,0])
        # return with the string converted to a float
        return (float(new[0])+float(new[1])/60.0+float(new[2])/3600.0) * direction[new_dir]    

def convert_WsWdtoUV(Ws,Wd):
    """
    Purpose:
     Convert wind speed and direction to U and V conponents.
     This routine follows the meteorological convention:
      - wind direction is positive going clockwise from north
      - U is positive towards east
      - V is positive towards north
    Usage:
     u,v = qcutils.convert_WsWdtoUV(Ws,Wd)
    Author: PRI
    Date: February 2015
    """
    u = -Ws*numpy.sin(numpy.radians(Wd))
    v = -Ws*numpy.cos(numpy.radians(Wd))
    return u,v

def convert_UVtoWsWd(u,v):
    """
    Purpose:
     Convert U and V conponents to wind speed and direction
     This routine follows the meteorological convention:
      - wind direction is positive going clockwise from north
      - U is positive towards east
      - V is positive towards north
    Usage:
     Ws,Wd = qcutils.convert_UVtoWsWd(U,V)
    Author: PRI
    Date: February 2015
    """
    Wd = float(270) - (numpy.degrees(numpy.arctan2(v,u)))
    Wd = numpy.mod(Wd,360)
    Ws = numpy.sqrt(u*u + v*v)
    return Ws,Wd

def CreateSeries(ds,Label,Data,FList=None,Flag=None,Attr=None):
    """
    Purpose:
     Create a series (1d array) of data in the data structure.
     If the series already exists in the data structure, data values and QC flags will be
     overwritten but attributes will be preserved.  However, the long_name and units attributes
     are treated differently.  The existing long_name will have long_name appended to it.  The
     existing units will be overwritten with units.
     This utility is the prefered method for creating or updating a data series because
     it implements a consistent method for creating series in the data structure.  Direct
     writes to the contents of the data structure are discouraged (unless PRI wrote the code:=P).
    Usage:
     Fsd,flag,attr = qcutils.GetSeriesasMA(ds,"Fsd")
      ... do something to Fsd here ...
     qcutils.CreateSeries(ds,"Fsd",Fsd,Flag=flag,Attr=attr)
    Author: PRI
    Date: Back in the day
    """
    ds.series['_tmp_'] = {}                       # create a temporary series to avoid premature overwrites
    # put the data into the temporary series
    if numpy.ma.isMA(Data):
        ds.series['_tmp_']['Data'] = numpy.ma.filled(Data,float(c.missing_value))
    else:
        ds.series['_tmp_']['Data'] = numpy.array(Data)
    # copy or make the QC flag
    if Flag is None:
        ds.series['_tmp_']['Flag'] = MakeQCFlag(ds,FList)
    else:
        ds.series['_tmp_']['Flag'] = Flag.astype(numpy.int32)
    # do the attributes
    ds.series['_tmp_']['Attr'] = {}
    if Label in ds.series.keys():                 # check to see if the series already exists
        for attr in ds.series[Label]['Attr']:     # if it does, copy the existing attributes
            if attr in Attr and ds.series[Label]['Attr'][attr]!=Attr[attr]:
                ds.series['_tmp_']['Attr'][attr] = Attr[attr]
            else:
                ds.series['_tmp_']['Attr'][attr] = ds.series[Label]['Attr'][attr]
    else:
        for item in Attr:
            ds.series['_tmp_']['Attr'][item] = Attr[item]
    ds.series[unicode(Label)] = ds.series['_tmp_']     # copy temporary series to new series
    del ds.series['_tmp_']                        # delete the temporary series

def CreateDatetimeRange(start,stop,step=datetime.timedelta(minutes=30)):
    '''
    Purpose:
     Create a series of datetimes between the "start" and "stop" datetimes
     and with a time step of "step".
    Useage:
     dt = ds.series['DateTime']['Data']
     ts = ds.globaleattributes['time_step']
     dt_evenlyspaced = CreateDatetimeRange(dt[0],dt[-1],step=datetime.timedelta(minutes=ts))]
    Author: PRI
    Date: December 2013
    '''
    result = []
    while start<stop:
        result.append(start)
        start = start + step
    return result

def file_exists(filename,mode="verbose"):
    if not os.path.exists(filename):
        if mode=="verbose":
            log.error(' File '+filename+' not found')
        return False
    else:
        return True

def FindIndicesOfBInA(a,b):
    """
    Purpose:
     Find the indices of elements in b that also occur in a.
     The routine is intended for use only with lists of Python datetime
     values.  This ensures the input series are monotonically increasing
     (though this is not a requirement) and contain no duplicates (which
     is required, or at least not handled).
    Limitations:
     Argument a is converted to a set to greatly speed the comparison
     of b elements with a.  This means that duplicates in a will be
     dropped and hence only 1 index will be returned for each value
     in b.
    Usage:
     indices = qcutils.FindIndicesOfBInA(a,b)
     where a is a list of Python datetime objects
           b is a list of Python datetime objects
           indices is a list of indices in b where the elements of b
                also occur in a
    Author: PRI
    Date: July 2015
    Comments: Replaces find_indices used up to V2.9.3.
    """
    if len(set(a))!=len(a):
        msg = " FindIndicesOfBInA: first argument contains duplicate values"
        log.warning(msg)
    tmpset = set(a)
    indices = [i for i,item in enumerate(b) if item in tmpset]
    return indices
    
def RemoveDuplicateRecords(ds):
    """ Remove duplicate records."""
    # the ds.series["DateTime"]["Data"] series is actually a list
    for item in ["DateTime","DateTime_UTC"]:
        if item in ds.series.keys():
            ldt,ldt_flag,ldt_attr = GetSeries(ds,item)
            # ldt_nodups is returned as an ndarray
            ldt_nodups,idx_nodups = numpy.unique(numpy.array(ldt),return_index=True)
            # now get ldt_nodups as a list
            ldt_nodups = ldt_nodups.tolist()
            # and put it back into the data structure
            ds.series[item]["Data"] = ldt_nodups
            ds.series[item]["Flag"] = ldt_flag[idx_nodups]
    # get a list of the series in the data structure
    series_list = [item for item in ds.series.keys() if '_QCFlag' not in item]
    # remove the DateTime
    for item in ["DateTime","DateTime_UTC"]:
        if item in series_list: series_list.remove(item)
    # loop over the series in the data structure
    for ThisOne in series_list:
        data_dups,flag_dups,attr = GetSeriesasMA(ds,ThisOne)
        data_nodups = data_dups[idx_nodups]
        flag_nodups = flag_dups[idx_nodups]
        CreateSeries(ds,ThisOne,data_nodups,Flag=flag_nodups,Attr=attr)
    ds.globalattributes['nc_nrecs'] = len(ds.series["DateTime"]["Data"])

def FixNonIntegralTimeSteps(ds,fixtimestepmethod=""):
    """
    Purpose:
     Fix time steps that are not an integral number of the default time step.
     The default time step is read from the "time_step" global attribute which is read from
     the L1 control file and written to the L1 netCDF file.
     The most common cause of non-integral time steps is drift in logger time stamp or
     rounding errors in Excel's treatment of datetimes.
    Usage:
     FixNonIntegralTimeSteps(ds)
    Called By: CheckTimeStep
    Author: PRI
    Date: February 2015
    To do:
     Implement [I]nterpolate
    """
    ts = int(ds.globalattributes["time_step"])
    ldt = ds.series["DateTime"]["Data"]
    dt_diffs = numpy.array([(ldt[i]-rounddttots(ldt[i],ts=ts)).total_seconds() for i in range(1,len(ldt))])
    log.info(" Maximum drift is "+str(numpy.max(dt_diffs))+" seconds, minimum drift is "+str(numpy.min(dt_diffs))+" seconds")
    ans = fixtimestepmethod
    if ans=="": ans = raw_input("Do you want to [Q]uit, [I]nterploate or [R]ound? ")
    if ans.lower()[0]=="q":
        print "Quiting ..."
        sys.exit()
    if ans.lower()[0]=="i":
        print "Interpolation to regular time step not implemented yet ..."
        sys.exit()
    if ans.lower()[0]=="r":
        log.info(" Rounding to the nearest time step")
        ldt_rounded = [rounddttots(dt,ts=ts) for dt in ldt]
        rdt = numpy.array([(ldt_rounded[i]-ldt_rounded[i-1]).total_seconds() for i in range(1,len(ldt))])
        log.info(" Maximum time step is now "+str(numpy.max(rdt))+" seconds, minimum time step is now "+str(numpy.min(rdt)))
        # replace the existing datetime series with the datetime series rounded to the nearest time step
        ds.series["DateTime"]["Data"] = ldt_rounded
    ds.globalattributes['nc_nrecs'] = len(ds.series["DateTime"]["Data"])
    
def FixTimeGaps(ds):
    """
    Purpose:
     Fix gaps in datetime series found by CheckTimeStep.
    Useage:
     has_gaps = CheckTimeStep(ds)
     if has_gaps:
         FixTimeGaps(ds)
    Author: PRI
    Date: April 2013
    Modified:
     September 2014 - rewrite for clarity and efficiency
     February 2015 - and again ...
    """
    ts = int(ds.globalattributes["time_step"])
    #ldt_gaps,ldt_flag,ldt_attr = GetSeries(ds,"DateTime")
    ldt_gaps = ds.series["DateTime"]["Data"]
    # generate a datetime list from the start datetime to the end datetime
    ldt_start = ldt_gaps[0]
    ldt_end = ldt_gaps[-1]
    ldt_nogaps = [result for result in perdelta(ldt_start,ldt_end,datetime.timedelta(minutes=ts))]
    # update the global attribute containing the number of records
    nRecs = len(ldt_nogaps)
    ds.globalattributes['nc_nrecs'] = nRecs
    # find the indices of the no-gap data in the original data
    idx_gaps = FindIndicesOfBInA(ldt_gaps,ldt_nogaps)
    # update the series of Python datetimes
    ds.series['DateTime']['Data'] = ldt_nogaps
    org_flag = ds.series['DateTime']['Flag'].astype(numpy.int32)
    ds.series['DateTime']['Flag'] = numpy.ones(nRecs,dtype=numpy.int32)
    ds.series['DateTime']['Flag'][idx_gaps] = org_flag
    # get a list of series in the data structure
    series_list = [item for item in ds.series.keys() if '_QCFlag' not in item]
    # remove the datetime-related series from data structure
    datetime_list = ["DateTime","DateTime_UTC"]
    for item in datetime_list:
        if item in series_list: series_list.remove(item)
    # now loop over the rest of the series in the data structure
    for ThisOne in series_list:
        data_nogaps = numpy.ones(nRecs,dtype=numpy.float64)*float(-9999)
        flag_nogaps = numpy.ones(nRecs,dtype=numpy.int32)
        data_gaps,flag_gaps,attr = GetSeriesasMA(ds,ThisOne)
        data_nogaps[idx_gaps] = data_gaps
        flag_nogaps[idx_gaps] = flag_gaps
        CreateSeries(ds,ThisOne,data_nogaps,Flag=flag_nogaps,Attr=attr)

def FixTimeStep(ds,fixtimestepmethod="round"):
    """
    Purpose:
     Fix problems with the time stamp.
    Useage:
     qcutils.FixTimeStep(ds,fixtimestepmethod=fixtimestepmethod)
    Author: PRI
    Date: April 2013
    Modified:
     February 2015 - split check and fix functions into different routines
    """
    # get the number of records
    nRecs = int(ds.globalattributes["nc_nrecs"])
    # get the time step
    ts = int(ds.globalattributes["time_step"])
    # time step between records in seconds
    dt = get_timestep(ds)
    dtmin = numpy.min(dt)
    dtmax = numpy.max(dt)
    if dtmin < ts*60:
        # duplicate or overlapping times found
        log.info(' FixTimeStep: duplicate or overlapping times found, removing ...')
        RemoveDuplicateRecords(ds)
        dt = get_timestep(ds)
        dtmin = numpy.min(dt)
        dtmax = numpy.max(dt)
        #log.info("After RemoveDuplicateRecords:"+str(dtmin)+" "+str(dtmax))
    if numpy.min(numpy.mod(dt,ts*60))!=0 or numpy.max(numpy.mod(dt,ts*60))!=0:
        # non-integral time steps found
        # indices of elements where time step not equal to default
        index = numpy.where(numpy.min(numpy.mod(dt,ts*60))!=0 or numpy.max(numpy.mod(dt,ts*60))!=0)[0]
        log.info(" FixTimeStep: Non-integral time steps found "+str(len(index))+" times out of "+str(nRecs))
        log.info(" FixTimeStep: Maximum time step was "+str(numpy.max(dt))+" seconds, minimum time step was "+str(numpy.min(dt)))
        FixNonIntegralTimeSteps(ds,fixtimestepmethod=fixtimestepmethod)
        dt = get_timestep(ds)
        dtmin = numpy.min(dt)
        dtmax = numpy.max(dt)
        #log.info("After FixNonIntegralTimeSteps:"+str(dtmin)+" "+str(dtmax))
    if dtmax > ts*60:
        # time gaps found
        log.info(' FixTimeStep: one or more time gaps found, inserting times ...')
        FixTimeGaps(ds)
        dt = get_timestep(ds)
        dtmin = numpy.min(dt)
        dtmax = numpy.max(dt)
        #log.info("After FixTimeGaps: "+str(dtmin)+" "+str(dtmax))

def GetAverageSeriesKeys(cf,ThisOne):
    if incf(cf,ThisOne) and haskey(cf,ThisOne,'AverageSeries'):
        if 'Source' in cf['Variables'][ThisOne]['AverageSeries'].keys():
            alist = ast.literal_eval(cf['Variables'][ThisOne]['AverageSeries']['Source'])
        else:
            log.error('  GetAverageSeriesKeys: key "Source" not in control file AverageSeries section for '+ThisOne)
            alist = []
        if 'standard_name' in cf['Variables'][ThisOne]['AverageSeries'].keys():
            standardname = str(cf['Variables'][ThisOne]['AverageSeries']['standard_name'])
        else:
            standardname = "not defined"
    else:
        standardname = "not defined"
        log.info('  GetAverageSeriesKeys: '+ThisOne+ ' not in control file or it does not have the "AverageSeries" key')
        alist = []
    return alist, standardname

def GetAltName(cf,ds,ThisOne):
    '''
    Check to see if the specified variable name is in the data structure (ds).
    If it is, return the variable name unchanged.
    If it isn't, check the control file to see if an alternate name has been specified
     and return the alternate name if one exists.
    '''
    if ThisOne not in ds.series.keys():
        if ThisOne in cf['Variables'].keys():
            ThisOne = cf['Variables'][ThisOne]['AltVarName']
            if ThisOne not in ds.series.keys():
                log.error('GetAltName: alternate variable name not in ds')
        else:
            log.error('GetAltName: cant find ',ThisOne,' in ds or control file')
    return ThisOne

def GetAltNameFromCF(cf,ThisOne):
    '''
    Get an alternate variable name from the control file.
    '''
    if ThisOne in cf['Variables'].keys():
        if 'AltVarName' in cf['Variables'][ThisOne].keys():
            ThisOne = str(cf['Variables'][ThisOne]['AltVarName'])
        else:
            print 'GetAltNameFromCF: AltVarName key not in control file for '+str(ThisOne)
    else:
        print 'GetAltNameFromCF: '+str(ThisOne)+' not in control file'
    return ThisOne

def GetAttributeDictionary(ds,ThisOne):
    attr = {}
    # if series ThisOne is in the data structure
    if ThisOne in ds.series.keys():
        attr = ds.series[ThisOne]['Attr']
    else:
        attr = MakeAttributeDictionary()
    return attr

def GetcbTicksFromCF(cf,ThisOne):
    '''
    Get colour bar tick labels from the control file.
    '''
    if ThisOne in cf['Variables'].keys():
        if 'Ticks' in cf['Variables'][ThisOne].keys():
            Ticks = eval(cf['Variables'][ThisOne]['Ticks'])
        else:
            print 'GetcbTicksFromCF: Ticks key not in control file for '+str(ThisOne)
    else:
        print 'GetcbTicksFromCF: '+str(ThisOne)+' not in control file'
    return Ticks

def GetRangesFromCF(cf,ThisOne,mode="verbose"):
    '''
    Get lower and upper range limits from the control file.
    '''
    if ThisOne in cf['Variables'].keys():
        if 'Lower' in cf['Variables'][ThisOne].keys():
            lower = float(cf['Variables'][ThisOne]['Lower'])
        else:
            if mode.lower()!="quiet":
                msg = "GetRangesFromCF: Lower key not in control file for "+str(ThisOne)
                log.info(msg)
            lower = None
        if 'Upper' in cf['Variables'][ThisOne].keys():
            upper = float(cf['Variables'][ThisOne]['Upper'])
        else:
            if mode.lower()!="quiet":
                msg = "GetRangesFromCF: Upper key not in control file for "+str(ThisOne)
                log.info(msg)
            upper = None
    else:
        if mode.lower()!="quiet":
            msg = "GetRangesFromCF: "+str(ThisOne)+" not in control file"
            log.info(msg)
        lower, upper = None
    return lower, upper

def GetDateIndex(dts,date,ts=30,default=0,match='exact'):
    """
    Purpose:
     Return the index of a date/datetime string in an array of datetime objects
    Usage:
     si = qcutils.GetDateIndex(datetimeseries,date_str,ts=30,default=0,match='exact')
    where
     dts      - array of datetime objects
     date_str - a date or date/time string in a format dateutils can parse
     ts       - time step for the data, optional (integer)
     default  - default value, optional (integer)
     match    - type of match (string) options are:
                "exact"            - finds the specified datetime and returns
                                     the index
                "startnextday"     - returns the index of the first time period
                                     in the next day
                "endpreviousday"   - returns the index of the last time period
                                     in the previous day
                "startnexthour"    - returns the index of the first time period
                                     in the next hour
                "endprevioushour"  - returns the index of the last time period
                                     in the previous hour
                "startnextmonth"   - returns the index of the first time period
                                     in the next month
                "endpreviousmonth" - returns the index of the last time period
                                     in the previous month
                NOTE: "startnextday" and "endpreviousday" can be used to pick
                    out time periods with an integer number of days
    Author: PRI
    Date: Back in the day
    """
    try:
        if len(date)!=0:
            i = dts.index(dateutil.parser.parse(date))
        else:
            if default==-1:
                i = len(dts)-1
            else:
                i = default
    except ValueError:
        if default==-1:
            i = len(dts)-1
        else:
            i = default
    if match=="exact":
        # if an exact match is required, do nothing
        pass
    elif match=="startnextmonth":
        # get to the start of the next day
        while abs(dts[i].hour+float(dts[i].minute)/60-float(ts)/60)>c.eps:
            i = i + 1
        while dts[i].day!=1:
            i = i + int(float(24)/(float(ts)/60))
    elif match=='startnextday':
        while abs(dts[i].hour+float(dts[i].minute)/60-float(ts)/60)>c.eps:
            i = i + 1
    elif match=="startnexthour":
        # check the time step value
        if int(ts)!=60:
            # if the time step is 60 then it is always the start of the next hour
            # we assume here that the time period ends on the datetime stamp
            while dts[i].minute!=ts:
                # iterate until the minutes equal the time step
                i = i + 1
    elif match=='endpreviousmonth':
        while abs(dts[i].hour+float(dts[i].minute)/60)>c.eps:
            i = i - 1
        while dts[i].day!=1:
            i = i - int(float(24)/(float(ts)/60))
    elif match=='endpreviousday':
        while abs(dts[i].hour+float(dts[i].minute)/60)>c.eps:
            i = i - 1
    elif match=="endprevioushour":
        # check the time step value
        if int(ts)!=60:
            # if the time step is 60 then it is always the end of the previous hour
            # we assume here that the time period ends on the datetime stamp
            while dts[i].minute!=0:
                # iterate until the minutes equal 0
                i = i - 1
    else:
        log.error("GetDateIndex: Unrecognised match option")
    return i

def GetGlobalAttributeValue(cf,ds,ThisOne):
    if ThisOne not in ds.globalattributes.keys():
        if ThisOne in cf['General'].keys():
            ds.globalattributes[ThisOne] = cf['General'][ThisOne]
        else:
            log.error('  GetGlobalAttributeValue: global attribute '+ThisOne+' was not found in the netCDF file or in the control file')
            ds.globalattributes[ThisOne] = None
    return ds.globalattributes[ThisOne]

def GetMergeSeriesKeys(cf,ThisOne,section=''):
    if len(section)==0: section = 'Variables'
    if 'Source' in cf[section][ThisOne]['MergeSeries'].keys():
        mlist = ast.literal_eval(cf[section][ThisOne]['MergeSeries']['Source'])
    else:
        log.error('  GetMergeSeriesKeys: key "Source" not in control file MergeSeries section for '+ThisOne)
        mlist = []
    if 'standard_name' in cf[section][ThisOne]['MergeSeries'].keys():
        standardname = str(cf[section][ThisOne]['MergeSeries']['standard_name'])
    else:
        standardname = 'not defined'
    return mlist, standardname

def GetPlotTitleFromCF(cf, nFig):
    if 'Plots' in cf:
        if str(nFig) in cf['Plots']:
            if 'Title' in cf['Plots'][str(nFig)]:
                Title = str(cf['Plots'][str(nFig)]['Title'])
            else:
                print 'GetPlotTitleFromCF: Variables key not in control file for plot '+str(nFig)
        else:
            print 'GetPlotTitleFromCF: '+str(nFig)+' key not in Plots section of control file'
    else:
        print 'GetPlotTitleFromCF: Plots key not in control file'
    return Title

def GetPlotVariableNamesFromCF(cf, n):
    if 'Plots' in cf:
        if str(n) in cf['Plots']:
            if 'Variables' in cf['Plots'][str(n)]:
                SeriesList = eval(cf['Plots'][str(n)]['Variables'])
            else:
                print 'GetPlotVariableNamesFromCF: Variables key not in control file for plot '+str(n)
        else:
            print 'GetPlotVariableNamesFromCF: '+str(n)+' key not in Plots section of control file'
    else:
        print 'GetPlotVariableNamesFromCF: Plots key not in control file'
    return SeriesList

def GetSeries(ds,ThisOne,si=0,ei=-1,mode="truncate"):
    """ Returns the data, QC flag and attributes of a series from the data structure."""
    # number of records
    nRecs = int(ds.globalattributes["nc_nrecs"])
    # check the series requested is in the data structure
    if ThisOne in ds.series.keys():
        # series is in the data structure
        if isinstance(ds.series[ThisOne]['Data'],list):
            # return a list if the series is a list
            Series = list(ds.series[ThisOne]['Data'])
        elif isinstance(ds.series[ThisOne]['Data'],numpy.ndarray):
            # return a numpy array if series is an array
            Series = ds.series[ThisOne]['Data'].copy()
        # now get the QC flag
        if 'Flag' in ds.series[ThisOne].keys():
            # return the QC flag if it exists
            Flag = ds.series[ThisOne]['Flag'].copy()
        else:
            # create a QC flag if one does not exist
            Flag = numpy.zeros(nRecs,dtype=numpy.int32)
        # now get the attribute dictionary
        if "Attr" in ds.series[ThisOne].keys():
            Attr = GetAttributeDictionary(ds,ThisOne)
        else:
            Attr = MakeAttributeDictionary()
    else:
        # make an empty series if the requested series does not exist in the data structure
        Series,Flag,Attr = MakeEmptySeries(ds,ThisOne)
    # tidy up
    if ei==-1: ei = nRecs - 1
    if mode=="truncate":
        # truncate to the requested start and end indices
        si = max(0,si)                  # clip start index at 0
        ei = min(nRecs,ei)              # clip end index to nRecs
        Series = Series[si:ei+1]        # truncate the data
        Flag = Flag[si:ei+1]            # truncate the QC flag
    elif mode=="pad":
        # pad with missing data at the start and/or the end of the series
        if si<0 and ei>nRecs-1:
            # pad at the start
            Series = numpy.append(float(c.missing_value)*numpy.ones(abs(si),dtype=numpy.float64),Series)
            Flag = numpy.append(numpy.ones(abs(si),dtype=numpy.int32),Flag)
            # pad at the end
            Series = numpy.append(Series,float(c.missing_value)*numpy.ones((ei-(nRecs-1)),dtype=numpy.float64))
            Flag = numpy.append(Flag,numpy.ones((ei-(nRecs-1)),dtype=numpy.int32))
        elif si<0 and ei<=nRecs-1:
            # pad at the start, truncate the end
            Series = numpy.append(float(c.missing_value)*numpy.ones(abs(si),dtype=numpy.float64),Series[:ei+1])
            Flag = numpy.append(numpy.ones(abs(si),dtype=numpy.int32),Flag[:ei+1])
        elif si>=0 and ei>nRecs-1:
            # truncate at the start, pad at the end
            Series = numpy.append(Series[si:],float(c.missing_value)*numpy.ones((ei-(nRecs-1)),numpy.float64))
            Flag = numpy.append(Flag[si:],numpy.ones((ei-(nRecs-1)),dtype=numpy.int32))
        elif si>=0 and ei<=nRecs-1:
            # truncate at the start and end
            Series = Series[si:ei+1]
            Flag = Flag[si:ei+1]
        else:
            msg = 'GetSeries: unrecognised combination of si ('+str(si)+') and ei ('+str(ei)+')'
            raise ValueError(msg)
    elif mode=="mirror":
        # reflect data about end boundaries if si or ei are out of bounds
        if si<0 and ei>nRecs-1:
            # mirror at the start
            Series = numpy.append(numpy.fliplr([Series[1:abs(si)+1]])[0],Series)
            Flag = numpy.append(numpy.fliplr([Flag[1:abs(si)+1]])[0],Flag)
            # mirror at the end
            sim = 2*nRecs-1-ei
            eim = nRecs-1
            Series = numpy.append(Series,numpy.fliplr([Series[sim:eim]])[0])
            Flag = numpy.append(Flag,numpy.fliplr([Flag[sim:eim]])[0])
        elif si<0 and ei<=nRecs-1:
            # mirror at start, truncate at end
            Series = numpy.append(numpy.fliplr([Series[1:abs(si)+1]])[0],Series[:ei+1])
            Flag = numpy.append(numpy.fliplr([Flag[1:abs(si)+1]])[0],Flag[:ei+1])
        elif si>=0 and ei>nRecs-1:
            # truncate at start, mirror at end
            sim = 2*nRecs-1-ei
            eim = nRecs
            Series = numpy.append(Series[si:],numpy.fliplr([Series[sim:eim]])[0])
            Flag = numpy.append(Flag[si:],numpy.fliplr([Flag[sim:eim]])[0])
        elif si>=0 and ei<=nRecs-1:
            # truncate at the start and end
            Series = Series[si:ei+1]
            Flag = Flag[si:ei+1]
        else:
            msg = 'GetSeries: unrecognised combination of si ('+str(si)+') and ei ('+str(ei)+')'
            raise ValueError(msg)            
    else:
        raise ValueError("GetSeries: unrecognised mode option "+str(mode))
    return Series,Flag,Attr

def MakeEmptySeries(ds,ThisOne):
    nRecs = int(ds.globalattributes['nc_nrecs'])
    Series = float(c.missing_value)*numpy.ones(nRecs,dtype=numpy.float64)
    Flag = numpy.ones(nRecs,dtype=numpy.int32)
    Attr = MakeAttributeDictionary()
    return Series,Flag,Attr

def GetSeriesasMA(ds,ThisOne,si=0,ei=-1,mode="truncate"):
    """
    Purpose:
     Returns a data series and the QC flag series from the data structure.
    Usage:
     data,flag,attr = qcutils.GetSeriesasMA(ds,label,si=0,ei=-1)
    where the arguments are;
      ds    - the data structure (dict)
      label - label of the data series in ds (string)
      si    - start index (integer), default 0
      ei    - end index (integer), default -1
    and the returned values are;
      data - values for the requested series in ds
             (numpy masked array, float64)
      flag - QC flag for the requested series in ds
             (numpy masked array, int32)
      attr - attribute dictionary for series
    Example:
     The code snippet below will return the incoming shortwave data values
     (Fsd) and the associated QC flag (f) as numpy masked arrays;
      ds = qcio.nc_read_series("HowardSprings_2011_L3.nc")
      Fsd,f,a = qcutils.GetSeriesasMA(ds,"Fsd")
    Author: PRI
    """
    Series,Flag,Attr = GetSeries(ds,ThisOne,si=si,ei=ei,mode=mode)
    Series,WasND = SeriestoMA(Series)
    return Series,Flag,Attr

def GetUnitsFromds(ds, ThisOne):
    units = ds.series[ThisOne]['Attr']['units']
    return units

def get_cfsection(cf,series='',mode='quiet'):
    '''
    Find the section in the control file that contains an entry for the series "series".
    USEAGE:  section = qcutils.get_cfsection(cf,series=<series_name>)
    INPUT:   cf            - a control file object (from ConfigObj)
             <series_name> - the name of the series (string)
    RETURNS: section       - the name of the section containing an entry for <series_name> (string)
    Note that the returned section name is an empty string if there is no entry for <series_name> in
    the control file.
    '''
    section = ''
    sectionlist = ['Variables','Drivers','Fluxes','Respiration','Partition','ER','GPP','NEE']
    if len(series)==0:
        msgtxt = ' get_cfsection: no input series specified'
        if mode!='quiet': log.info(msgtxt)
        return section
    for ThisSection in sectionlist:
        if ThisSection in cf.keys():
            if series in cf[ThisSection]: section = ThisSection
    if len(section)==0:
        msgtxt = ' get_cfsection: series '+str(series)+' not found in control file'
        if mode!='quiet': log.info(msgtxt)
    return section

def get_coverage_groups(ds,rad=None,met=None,flux=None,soil=None):
    level = "L1"
    if "nc_level" in ds.globalattributes:
        level = str(ds.globalattributes["nc_level"])
    rad = ['Fsd','Fsu','Fld','Flu','Fn']
    met = ['Ah','Cc','Precip','ps','Ta','Ws','Wd']
    flux = ['Fm','ustar','Fh','Fe','Fc']
    soil = ['Fg','Ts','Sws']
    for ThisGroup, ThisLabel in zip([rad,met,flux,soil],['radiation','meteorology','flux','soil']):
        sum_coverage = float(0); count = float(0)
        for ThisOne in ThisGroup:
            if ThisOne in ds.series.keys():
                sum_coverage = sum_coverage + float(ds.series[ThisOne]['Attr']['coverage_'+level])
                count = count + 1
        if count!=0:
            coverage_group = sum_coverage/count
        else:
            coverage_group = 0
        ds.globalattributes['coverage_'+ThisLabel+'_'+level] = str('%d'%coverage_group)

def get_coverage_individual(ds):
    level = "L1"
    if "nc_level" in ds.globalattributes:
        level = str(ds.globalattributes["nc_level"])
    SeriesList = ds.series.keys()
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in SeriesList: SeriesList.remove(ThisOne)
    for ThisOne in SeriesList:
        num_good = len(numpy.where(abs(ds.series[ThisOne]['Data']-float(c.missing_value))>c.eps)[0])
        coverage = 100*float(num_good)/float(ds.globalattributes['nc_nrecs'])
        ds.series[ThisOne]['Attr']['coverage_'+level] = str('%d'%coverage)

def get_datetimefromnctime(ds,time,time_units):
    """
    Purpose:
     Create a series of datetime objects from the time read from a netCDF file.
    Usage:
     qcutils.get_datetimefromnctime(ds,time,time_units)
    Side effects:
     Creates a Python datetime series in the data structure
    Author: PRI
    Date: September 2014
    """
    ts = int(ds.globalattributes["time_step"])
    nRecs = int(ds.globalattributes["nc_nrecs"])
    dt = netCDF4.num2date(time,time_units)
    ds.series[unicode("DateTime")] = {}
    ds.series["DateTime"]["Data"] = list(dt)
    ds.series["DateTime"]["Flag"] = numpy.zeros(nRecs)
    ds.series["DateTime"]["Attr"] = {}
    ds.series["DateTime"]["Attr"]["long_name"] = "Datetime in local timezone"
    ds.series["DateTime"]["Attr"]["units"] = "None"

def get_datetimefromxldate(ds):
    ''' Creates a series of Python datetime objects from the Excel date read from the Excel file.
        Thanks to John Machin for the quick and dirty code
         see http://stackoverflow.com/questions/1108428/how-do-i-read-a-date-in-excel-format-in-python'''

    log.info(' Getting the Python datetime series from the Excel datetime')
    xldate = ds.series['xlDateTime']['Data']
    nRecs = len(ds.series['xlDateTime']['Data'])
    datemode = int(ds.globalattributes['xl_datemode'])
    ds.series[unicode('DateTime')] = {}
    ds.series['DateTime']['Data'] = [None]*nRecs
    basedate = datetime.datetime(1899, 12, 30)
    #ldt = [basedate + datetime.timedelta(days=xldate[i] + 1462 * datemode) for i in range(nRecs)]
    #ds.series['DateTime']['Data'][i] = ldt
    for i in range(nRecs):
        ds.series['DateTime']['Data'][i] = basedate + datetime.timedelta(days=xldate[i] + 1462 * datemode)
    ds.series['DateTime']['Flag'] = numpy.zeros(nRecs)
    ds.series['DateTime']['Attr'] = {}
    ds.series['DateTime']['Attr']['long_name'] = 'Datetime in local timezone'
    ds.series['DateTime']['Attr']['units'] = 'None'

def get_datetimefromymdhms(ds):
    ''' Creates a series of Python datetime objects from the year, month,
    day, hour, minute and second series stored in the netCDF file.'''
    SeriesList = ds.series.keys()
    if 'Year' not in SeriesList or 'Month' not in SeriesList or 'Day' not in SeriesList or 'Hour' not in SeriesList or 'Minute' not in SeriesList or 'Second' not in SeriesList:
        log.info(' get_datetimefromymdhms: unable to find all datetime fields required')
        return
    log.info(' Getting the date and time series')
    nRecs = get_nrecs(ds)
    ts = ds.globalattributes["time_step"]
    ds.series[unicode('DateTime')] = {}
    ds.series['DateTime']['Data'] = [None]*nRecs
    if "Microseconds" in ds.series.keys():
        microseconds = ds.series["Microseconds"]["Data"]
    else:
        microseconds = numpy.zeros(nRecs,dtype=numpy.float64)
    for i in range(nRecs):
        #print i,int(ds.series['Year']['Data'][i]),int(ds.series['Month']['Data'][i]),int(ds.series['Day']['Data'][i])
        #print i,int(ds.series['Hour']['Data'][i]),int(ds.series['Minute']['Data'][i]),int(ds.series['Second']['Data'][i])
        ds.series['DateTime']['Data'][i] = datetime.datetime(int(ds.series['Year']['Data'][i]),
                                                       int(ds.series['Month']['Data'][i]),
                                                       int(ds.series['Day']['Data'][i]),
                                                       int(ds.series['Hour']['Data'][i]),
                                                       int(ds.series['Minute']['Data'][i]),
                                                       int(ds.series['Second']['Data'][i]),
                                                       int(microseconds[i]))
    ds.series['DateTime']['Flag'] = numpy.zeros(nRecs)
    ds.series['DateTime']['Attr'] = {}
    ds.series['DateTime']['Attr']['long_name'] = 'Date-time object'
    ds.series['DateTime']['Attr']['units'] = 'None'

def get_keyvaluefromcf(cf,sections,key,default=None,mode="quiet"):
    """
    Purpose:
     General return a keyword value from a control file.
    Usage:
     keyval = qcutils.get_keyvaluefromcf(cf,sections,key,default=default)
     where
      cf is a control file object from ConfigObj
      sections is a list of sections and nested sub-sections to search
      key is the keyword
      default is a default value
    Example:
     ncOutFileName = qcutils.get_keyvaluefromcf(cf,["Files","Out"],"ncFileName",default="")
     The example above will return the value for ncFileName from the ["Files"]["Out"] sub-section
     in the control file.
    Author: PRI
    Date: February 2015
    """
    if len(sections)<1:
        msg = " get_keyvaluefromsections: no sections specified"
        if mode.lower()!="quiet": log.info(msg)
    if sections[0] in cf:
        section = cf[sections[0]]
        if len(sections)>1:
            for item in sections[1:]:
                if item in section:
                    section = section[item]
                else:
                    msg = " get_keyvaluefromcf: Sub section "+item+" not found in control file, used default ("+str(default)+")"
                    if mode.lower()!="quiet": log.info(msg)
                    value = default
        if key in section:
            value = section[key]
        else:
            msg = " get_keyvaluefromcf: Key "+key+" not found in section, used default ("+str(default)+")"
            if mode.lower()!="quiet": log.info(msg)
            value = default
    else:
        msg = " get_keyvaluefromcf: Section "+sections[0]+" not found in control file, used default ("+str(default)+")"
        if mode.lower()!="quiet": log.error(msg)
        value = default
    return value

def get_missingingapfilledseries(ds):
    """
    Purpose:
     Check series in data structure and print a message to the screen if missing points are found.
    Usage:
     gfalternate_checkformissing(ds,series_list=series_list)
      where ds is a data structure
            series_list is a list of series to check
    Author: PRI
    Date: March 2015
    """
    # get a local pointer to the datetime
    ldt = ds.series["DateTime"]["Data"]
    # create an empty list
    alt_list = []
    # check to see if there was any gap filling using data from alternate sources
    if "alternate" in dir(ds):
        # if so, get a list of the quantities gap filled from alternate sources
        alt_list = list(set([ds.alternate[item]["label_tower"] for item in ds.alternate.keys()]))
    # create an empty list
    cli_list = []
    # check to see if there was any gap filling from climatology
    if "climatology" in dir(ds):
        # if so, get a list of the quantities gap filled using climatology
        cli_list = list(set([ds.climatology[item]["label_tower"] for item in ds.climatology.keys()]))
    # one list to rule them, one list to bind them ...
    gf_list = list(set(alt_list+cli_list))
    # clear out if there was no gap filling
    if len(gf_list)==0: return
    # loop over the series to be checked
    for series in gf_list:
        if series not in ds.series.keys(): continue
        data,flag,attr = GetSeriesasMA(ds,series)
        idx = numpy.ma.where(data.mask==True)[0]
        if len(idx)!=0:
            msg = " Missing points ("+str(len(idx))+") found in "+series
            log.error(msg)
            #ldt_missing = [ldt[i] for i in idx]
            #msg = " The first 10 missing data is at datetimes "+str(ldt_missing[0:9])
            #log.error(msg)

def get_nrecs(ds):
    if 'nc_nrecs' in ds.globalattributes.keys():
        nRecs = int(ds.globalattributes['nc_nrecs'])
    elif 'NumRecs' in ds.globalattributes.keys():
        nRecs = int(ds.globalattributes['NumRecs'])
    else:
        series_list = ds.series.keys()
        nRecs = len(ds.series[series_list[0]]['Data'])
    return nRecs

def get_timestep(ds):
    """
    Purpose:
     Return an array of time steps in seconds between records
    Useage:
     dt = qcutils.get_timestep(ds)
    Author: PRI
    Date: February 2015
    """
    # local pointer to the Python datetime series
    ldt = ds.series["DateTime"]["Data"]
    # time step between records in seconds
    dt = numpy.array([(ldt[i]-ldt[i-1]).total_seconds() for i in range(1,len(ldt))])
    return dt

def get_timezone(site_name,prompt="no"):
    """ Return the time zone based on the site name."""
    time_zone = ""
    found = False
    # strip out spaces and commas from the site name
    site_name = site_name.replace(" ","").replace(",","")
    for item in c.tz_dict.keys():
        if item in site_name.lower():
            time_zone = c.tz_dict[item]
            found = True
        else:
            # cant find the site in the dictionary so ask the user
            if prompt.lower()=="yes":
                root = Tkinter.Tk(); root.withdraw()
                time_zone = tkSimpleDialog.askstring("Time zone","Enter time zone eg Australia/Melbourne")
                root.destroy()
                found = True
    return time_zone,found

def get_UTCfromlocaltime(ds):
    '''
    Purpose:
     Creates a UTC datetime series in the data structure from the
     local datetime series.
    Usage:
     ldt_UTC = qcutils.get_UTCfromlocaltime(ds)
    Assumptions:
     No daylight savings used in the local datetime
    Author: PRI
    '''
    # check the time_zone global attribute is set, we cant continue without it
    if "time_zone" not in ds.globalattributes.keys():
        log.warning("get_UTCfromlocaltime: time_zone not in global attributes, checking elsewhere ...")
        if "site_name" in ds.globalattributes.keys():
            site_name = ds.globalattributes["site_name"]
        else:
            log.warning("get_UTCfromlocaltime: site_name not in global attributes, skipping UTC calculation ...")
            return
        time_zone,found = get_timezone(site_name,prompt="no")
        if not found:
            log.warning("get_UTCfromlocaltime: site_name not in time zone dictionary")
            return
        else:
            log.info("get_UTCfromlocaltime: time_zone found in time zone dictionary")
            ds.globalattributes["time_zone"] = time_zone
    log.info(' Getting the UTC datetime from the local datetime')
    # get the number of records
    nRecs = len(ds.series['xlDateTime']['Data'])
    # get the time zone
    tz = ds.globalattributes["time_zone"]
    # create a timezone object
    loc_tz = pytz.timezone(tz)
    # local pointer to the datetime series in ds
    ldt = ds.series["DateTime"]["Data"]
    # localise the datetime
    ldt_loc = [loc_tz.localize(dt) for dt in ldt]
    # convert to UTC
    ldt_utc = [dt.astimezone(pytz.utc) for dt in ldt_loc]
    return ldt_utc

def get_xldatefromdatetime(ds):
    '''
    Purpose:
     Returns a list of xldatetime (floating point number represent decimal days
     since 00:00 1/1/1900) from a list of Python datetimes
    Usage:
     qcutils.get_xldatefromdatetime(ds)
    Assumptions:
     The Excel datetime series ("xlDateTime") exists in the data structure ds.
    Author: PRI
    '''
    # get the datemode of the original Excel spreadsheet
    if "xl_datemode" in ds.globalattributes.keys():
        datemode = int(ds.globalattributes["xl_datemode"])
    else:
        datemode = int(0)
    nRecs = int(ds.globalattributes["nc_nrecs"])
    # get the Excel datetime series, flag and attributes
    if "xlDateTime" in ds.series.keys():
        xldt_org,xldt_flag,xldt_attr = GetSeriesasMA(ds,"xlDateTime")
    else:
        xldt_flag = numpy.zeros(nRecs,dtype=numpy.int32)
        xldt_attr = MakeAttributeDictionary(long_name="Date/time in Excel format",units="days since 1899-12-31 00:00:00")
    # get a local pointer to the Python DateTime series in ds
    ldt = ds.series["DateTime"]["Data"]
    # get a list of Excel datetimes from the Python datetime objects
    xldate = [xlrd.xldate.xldate_from_datetime_tuple((ldt[i].year,
                                                      ldt[i].month,
                                                      ldt[i].day,
                                                      ldt[i].hour,
                                                      ldt[i].minute,
                                                      ldt[i].second),
                                                      datemode) for i in range(0,len(ldt))]
    xldt_new = numpy.ma.array(xldate, dtype=numpy.float64)
    # overwrite the existing Excel datetime series
    CreateSeries(ds,"xlDateTime",xldt_new,Flag=xldt_flag,Attr=xldt_attr)

def get_ymdhmsfromdatetime(ds):
    '''
    Purpose:
     Gets the year, month, day, hour, minute and second from a list of
     Python datetimes.  The Python datetime series is read from
     the input data structure and the results are written back to the
     data structure.
    Usage:
     qcutils.get_ymdhmsfromdatetime(ds)
    Assumptions:
     None
    Author: PRI
    '''
    nRecs = int(ds.globalattributes["nc_nrecs"])
    dt = ds.series["DateTime"]["Data"]
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    Year = numpy.array([dt[i].year for i in range(0,nRecs)]).astype(numpy.int32)
    Month = numpy.array([dt[i].month for i in range(0,nRecs)]).astype(numpy.int32)
    Day = numpy.array([dt[i].day for i in range(0,nRecs)]).astype(numpy.int32)
    Hour = numpy.array([dt[i].hour for i in range(0,nRecs)]).astype(numpy.int32)
    Minute = numpy.array([dt[i].minute for i in range(0,nRecs)]).astype(numpy.int32)
    Second = numpy.array([dt[i].second for i in range(0,nRecs)]).astype(numpy.int32)
    Hdh = numpy.array([float(Hour[i])+float(Minute[i])/60. for i in range(0,nRecs)]).astype(numpy.float64)
    Ddd = numpy.array([(dt[i] - datetime.datetime(Year[i],1,1)).days+1+Hdh[i]/24. for i in range(0,nRecs)]).astype(numpy.float64)
    CreateSeries(ds,'Year',Year,Flag=flag,Attr=MakeAttributeDictionary(long_name='Year',units='none'))
    CreateSeries(ds,'Month',Month,Flag=flag,Attr=MakeAttributeDictionary(long_name='Month',units='none'))
    CreateSeries(ds,'Day',Day,Flag=flag,Attr=MakeAttributeDictionary(long_name='Day',units='none'))
    CreateSeries(ds,'Hour',Hour,Flag=flag,Attr=MakeAttributeDictionary(long_name='Hour',units='none'))
    CreateSeries(ds,'Minute',Minute,Flag=flag,Attr=MakeAttributeDictionary(long_name='Minute',units='none'))
    CreateSeries(ds,'Second',Second,Flag=flag,Attr=MakeAttributeDictionary(long_name='Second',units='none'))
    CreateSeries(ds,'Hdh',Hdh,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal hour of the day',units='none'))
    CreateSeries(ds,'Ddd',Ddd,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal day of the year',units='none'))

def get_ymdhmsfromxldate(ds):
    """
        Gets year, month, day, hour, and if available seconds, from
        excel-formatted Timestamp
        
        Usage qcts.get_ymdhmsfromxldate(ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting date and time variables')
    # get the date mode of the original Excel datetime
    datemode = int(ds.globalattributes['xl_datemode'])
    nRecs = len(ds.series['xlDateTime']['Data'])
    Year = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Month = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Day = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Hour = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Minute = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Second = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Hdh = numpy.array([c.missing_value]*nRecs,numpy.float64)
    Ddd = numpy.array([c.missing_value]*nRecs,numpy.float64)
    flag = numpy.zeros(nRecs)
    for i in range(nRecs):
        DateTuple = xlrd.xldate_as_tuple(ds.series['xlDateTime']['Data'][i],datemode)
        Year[i] = int(DateTuple[0])
        Month[i] = int(DateTuple[1])
        Day[i] = int(DateTuple[2])
        Hour[i] = int(DateTuple[3])
        Minute[i] = int(DateTuple[4])
        Second[i] = int(DateTuple[5])
        Hdh[i] = float(DateTuple[3])+float(DateTuple[4])/60.
        Ddd[i] = ds.series['xlDateTime']['Data'][i] - xlrd.xldate.xldate_from_date_tuple((Year[i],1,1),datemode) + 1
    CreateSeries(ds,'Year',Year,Flag=flag,Attr=MakeAttributeDictionary(long_name='Year',units='none'))
    CreateSeries(ds,'Month',Month,Flag=flag,Attr=MakeAttributeDictionary(long_name='Month',units='none'))
    CreateSeries(ds,'Day',Day,Flag=flag,Attr=MakeAttributeDictionary(long_name='Day',units='none'))
    CreateSeries(ds,'Hour',Hour,Flag=flag,Attr=MakeAttributeDictionary(long_name='Hour',units='none'))
    CreateSeries(ds,'Minute',Minute,Flag=flag,Attr=MakeAttributeDictionary(long_name='Minute',units='none'))
    CreateSeries(ds,'Second',Second,Flag=flag,Attr=MakeAttributeDictionary(long_name='Second',units='none'))
    CreateSeries(ds,'Hdh',Hdh,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal hour of the day',units='none'))
    CreateSeries(ds,'Ddd',Ddd,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal day of the year',units='none'))

def haskey(cf,ThisOne,key):
    return key in cf['Variables'][ThisOne].keys()

def incf(cf,ThisOne):
    return ThisOne in cf['Variables'].keys()

def linear_function(B,x):
    """
    Purpose:
     Linear function for use with orthogonal distance regression.
    Usage:
     linear = scipy.odr.Model(qcutils.linear_function)
     where B is a list of slope and offset values
           x is an array of x values
    """
    return B[0]*x + B[1]

def MakeAttributeDictionary(**kwargs):
    """
    Purpose:
     Make an attribute dictionary.
    Usage:
     attr_new = qcutils.MakeAttributeDictionary(long_name = "some string",attr_exist)
     where long_name is an attribute to be written to the new attribute dictionary
           attr_exist is an existing attribute dictionary
    Author: PRI
    Date: Back in the day
    """
    default_list = ['ancillary_variables','height','instrument','serial_number','standard_name','long_name','units']
    attr = {}
    for item in kwargs:
        if isinstance(item, dict):
            for entry in item: attr[entry] = item[entry]
        else:
            attr[item] = kwargs.get(item,'not defined')
        if item in default_list: default_list.remove(item)
    if len(default_list)!=0:
        for item in default_list: attr[item] = 'not defined'
    attr["missing_value"] = c.missing_value
    return attr

def MakeQCFlag(ds,SeriesList):
    flag = []
    if len(SeriesList)<=0:
        #log.info('  MakeQCFlag: no series list specified')
        pass
    if len(SeriesList)==1:
        if SeriesList[0] in ds.series.keys():
            flag = ds.series[SeriesList[0]]['Flag'].copy()
        else:
            log.error('  MakeQCFlag: series '+str(SeriesList[0])+' not in ds.series')
    if len(SeriesList)>1:
        for ThisOne in SeriesList:
            if ThisOne in ds.series.keys():
                if len(flag)==0:
                    #flag = numpy.ones(numpy.size(ds.series[ThisOne]['Flag']))
                    flag = ds.series[ThisOne]['Flag'].copy()
                else:
                    tmp_flag = ds.series[ThisOne]['Flag'].copy()      # get a temporary copy of the flag
                    index = numpy.where(numpy.mod(tmp_flag,10)==0)    # find the elements with flag = 0, 10, 20 etc
                    tmp_flag[index] = 0                               # set them all to 0
                    flag = numpy.maximum(flag,tmp_flag)               # now take the maximum
            else:
                log.error('  MakeQCFlag: series '+ThisOne+' not in ds.series')
    return flag.astype(numpy.int32)

def MAtoSeries(Series):
    """
    Convert a masked array to a numpy ndarray with masked elements set to c.missing_value.
    Useage:
     Series, WasMA = MAtoSeries(Series)
     where:
      Series (input)    is the data series to be converted.
      WasMA  (returned) is a logical, True if the input series was a masked array.
      Series (output)   is the input series convered to an ndarray with c.missing_value values
                        for missing data.
    """
    WasMA = False
    if numpy.ma.isMA(Series):
        WasMA = True
        Series = numpy.ma.filled(Series,float(c.missing_value))
    return Series, WasMA

def MergeQCFlag(QCFlag_list):
    """ Merge a list of QC flags by taking the element-wise maximum."""
    if len(QCFlag_list)==0: return None
    if len(QCFlag_list)==1: return QCFlag_list[0]
    flag = QCFlag_list[0].copy()                            # get a copy of the first flag
    for item in QCFlag_list[1:]:                            # loop over the list of flags
        tmp_flag = item.copy()                              # get a copy of the next flag
        index = numpy.where(numpy.mod(tmp_flag,10)==0)      # find the elements with flag = 0, 10, 20 etc
        tmp_flag[index] = 0                                 # set them all to 0
        flag = numpy.maximum(flag,tmp_flag)                 # now take the maximum
    return flag

def nxMom_nxScalar_alpha(zoL):
    nRecs = numpy.size(zoL)
    nxMom = numpy.ma.ones(nRecs) * 0.079
    nxScalar = numpy.ma.ones(nRecs) * 0.085
    alpha = numpy.ma.ones(nRecs) * 0.925
    #  get the index of stable conditions
    stable = numpy.ma.where(zoL>0)[0]
    #  now set the series to their stable values
    nxMom[stable] = 0.079 * (1 + 7.9 * zoL[stable]) ** 0.75
    nxScalar[stable] = 2.0 - 1.915 / (1 + 0.5 * zoL[stable])
    alpha[stable] = 1
    return nxMom, nxScalar, alpha

def perdelta(start, end, delta):
    """
    Yields an iterator of datetime objects from start to end with time step delta.
    """
    curr = start
    while curr <= end:
        yield curr
        curr += delta

def polyval(p,x):
    """
    Replacement for the polyval routine in numpy.  This version doesnt check the
    input variables to make sure they are array_like.  This means that when
    masked arrays are treated correctly when they are passed to this routine.
    Parameters
    ----------
     p : a 1D array of coefficients, highest order first
     x : a 1D array of points at which to evaluate the polynomial described by
         the coefficents in p
    Example
    -------
    >>> x = numpy.array([1,2,3])
    >>> p = numpy.array([2,0])
    >>> qcutils.polyval(p,x)
        array([2,4,6])
    >>> y = numpy.array([1,c.missing_value,3])
    >>> y = numpy.ma.masked_where(y==c.missing_value,y)
    >>> qcutils.polyval(p,y)
    masked_array(data = [2 -- 6],
                 mask = [False True False],
                 fill_value = 999999)
    """
    y = 0
    for i in range(len(p)):
        y = x*y + p[i]
    return y

def rounddttots(dt,ts=30):
    dt += datetime.timedelta(minutes=int(ts/2))
    dt -= datetime.timedelta(minutes=dt.minute % int(ts),seconds=dt.second,microseconds=dt.microsecond)
    return dt    

def rounddttoseconds(dt):
    dt += datetime.timedelta(seconds=0.5)
    dt -= datetime.timedelta(seconds=dt.second % 1,microseconds=dt.microsecond)
    return dt

def round_datetime(ds,mode="nearest_timestep"):
    """
    Purpose:
     Round the series of Python datetimes to the nearest time based on mode
    Usage:
     qcutils.round_datetime(ds,mode=mode)
     where;
      mode = "nearest_second" rounds to the nearesy second
      mode = "nearest_timestep" rounds to the nearest time step
    Author: PRI
    Date: February 2015
    """
    # local pointer to the datetime series
    ldt = ds.series["DateTime"]["Data"]
    # check which rounding option has been chosen
    if mode.lower()=="nearest_timestep":
        # get the time step
        if "time_step" in ds.globalattributes:
            ts = int(ds.globalattributes["time_step"])
        else:
            ts = numpy.mean(get_timestep(ds)/60)
            ts = roundtobase(ts,base=30)
            ds.globalattributes["time_step"] = ts
        # round to the nearest time step
        rldt = [rounddttots(dt,ts=ts) for dt in ldt]
    elif mode.lower()=="nearest_second":
        # round to the nearest second
        rldt = [rounddttoseconds(dt) for dt in ldt]
    else:
        # unrecognised option for mode, return original datetime series
        log.error(" round_datetime: unrecognised mode ("+str(mode)+")"+" ,returning original time series")
        rldt = ds.series["DateTime"]["Data"]
    # replace the original datetime series with the rounded one
    ds.series["DateTime"]["Data"] = rldt

def roundtobase(x,base=5):
    return int(base*round(float(x)/base))

def round2sig(x,sig=2):
    '''
    Round a float to a specified number of significant digits (default is 2).
    '''
    return round(x, sig-int(math.floor(math.log10(abs(x))))-1)

def r(b, p, alpha):
    """
    Function to calculate the r coeficient of the Massman frequency correction.
    """
    r = ((b ** alpha) / (b ** alpha + 1)) * \
           ((b ** alpha) / (b ** alpha + p ** alpha)) * \
           (1 / (p ** alpha + 1))
    return r

def SeriestoMA(Series):
    """
    Convert a numpy ndarray to a masked array.
    Useage:
     Series, WasND = SeriestoMA(Series)
     where:
      Series (input)    is the data series to be converted.
      WasND  (returned) is a logical, True if the input series was an ndarray
      Series (output)   is the input series convered to a masked array.
    """
    WasND = False
    if not numpy.ma.isMA(Series):
        WasND = True
        Series = numpy.ma.masked_where(abs(Series-numpy.float64(c.missing_value))<c.eps,Series)
    return Series, WasND

def SetUnitsInds(ds, ThisOne, units):
    ds.series[ThisOne]['Attr']['units'] = units

def startlog(loggername,loggerfile):
    logger = logging.getLogger(loggername)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(loggerfile)
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', '%H:%M:%S')
    #formatter = logging.Formatter('%(asctime)s %(name)-8s %(levelname)-6s %(message)s', '%d-%m-%y %H:%M')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def UpdateGlobalAttributes(cf,ds,level):
    ds.globalattributes["nc_level"] = str(level)
    ds.globalattributes["EPDversion"] = sys.version
    # put the control file name into the global attributes
    ds.globalattributes["controlfile_name"] = cf["controlfile_name"]
    if "Global" in cf:
        for item in cf["Global"].keys():
            if item not in ds.globalattributes.keys():
                ds.globalattributes[item] = cf["Global"][item].replace("\n"," ").replace("\r","")
