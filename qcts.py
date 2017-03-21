"""
    QC Data Function Module
    Used to perform the tasks queued by qcls.py
    """

import sys
import ast
import constants as c
import datetime
import inspect
from matplotlib.dates import date2num
import meteorologicalfunctions as mf
import numpy
import os
import qcck
import qcfunc
import qcio
import qcutils
from scipy import interpolate, signal
import time
import xlrd
from matplotlib.mlab import griddata
import xlwt
import logging
import pysolar

log = logging.getLogger('qc.ts')

def albedo(cf,ds):
    """
        Filter albedo measurements to:
            high solar angle specified by periods between 10.00 and 14.00, inclusive
            and
            full sunlight in which Fsd > 290 W/m2
        
        Usage qcts.albedo(ds)
        ds: data structure
        """
    log.info(' Applying albedo constraints')
    if 'albedo' not in ds.series.keys():
        if 'Fsd' in ds.series.keys() and 'Fsu' in ds.series.keys():
            Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
            Fsu,f,a = qcutils.GetSeriesasMA(ds,'Fsu')
            albedo = Fsu / Fsd
            attr = qcutils.MakeAttributeDictionary(long_name='solar albedo',units='none',standard_name='solar_albedo')
            qcutils.CreateSeries(ds,'albedo',albedo,FList=['Fsd','Fsu'],Attr=attr)
        else:
            log.warning('  Fsd or Fsu not in ds, albedo not calculated')
            return
    else:
        albedo,f,a = qcutils.GetSeriesasMA(ds,'albedo')
        if 'Fsd' in ds.series.keys():
            Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
        else:
            Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fn')
    
    if qcutils.cfkeycheck(cf,ThisOne='albedo',key='Threshold'):
        Fsdbase = float(cf['Variables']['albedo']['Threshold']['Fsd'])
        ds.series['albedo']['Attr']['FsdCutoff'] = Fsdbase
    else:
        Fsdbase = 290.
    index = numpy.ma.where((Fsd < Fsdbase) | (ds.series['Hdh']['Data'] < 10) | (ds.series['Hdh']['Data'] > 14))[0]
    index1 = numpy.ma.where(Fsd < Fsdbase)[0]
    index2 = numpy.ma.where((ds.series['Hdh']['Data'] < 10) | (ds.series['Hdh']['Data'] > 14))[0]
    albedo[index] = numpy.float64(c.missing_value)
    ds.series['albedo']['Flag'][index1] = numpy.int32(51)     # bad Fsd flag only if bad time flag not set
    ds.series['albedo']['Flag'][index2] = numpy.int32(52)     # bad time flag
    ds.series['albedo']['Data']=numpy.ma.filled(albedo,float(c.missing_value))

def ApplyLinear(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. Time period
        to apply the correction, slope and offset are specified in the control
        file.
        
        Usage qcts.ApplyLinear(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    if ThisOne not in ds.series.keys(): return
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'Linear'):
        log.info('  Applying linear correction to '+ThisOne)
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(c.missing_value),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag'].copy()
        ldt = ds.series['DateTime']['Data']
        LinearList = cf['Variables'][ThisOne]['Linear'].keys()
        for i in range(len(LinearList)):
            LinearItemList = ast.literal_eval(cf['Variables'][ThisOne]['Linear'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(LinearItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(LinearItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = float(LinearItemList[2])
            Offset = float(LinearItemList[3])
            data[si:ei] = Slope * data[si:ei] + Offset
            index = numpy.where(flag[si:ei]==0)[0]
            flag[si:ei][index] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,float(c.missing_value)).astype(numpy.float64)
            ds.series[ThisOne]['Flag'] = flag

def ApplyLinearDrift(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. The slope is
        interpolated for each 30-min period between the starting value at time 0
        and the ending value at time 1.  Slope0, Slope1 and Offset are defined
        in the control file.  This function applies to a dataset in which the
        start and end times in the control file are matched by the time period
        in the dataset.
        
        Usage qcts.ApplyLinearDrift(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    if ThisOne not in ds.series.keys(): return
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'Drift'):
        log.info('  Applying linear drift correction to '+ThisOne)
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(c.missing_value),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag']
        ldt = ds.series['DateTime']['Data']
        DriftList = cf['Variables'][ThisOne]['Drift'].keys()
        for i in range(len(DriftList)):
            DriftItemList = ast.literal_eval(cf['Variables'][ThisOne]['Drift'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(DriftItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(DriftItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.zeros(len(data))
            Slope0 = float(DriftItemList[2])
            Slope1 = float(DriftItemList[3])
            Offset = float(DriftItemList[4])
            nRecs = len(Slope[si:ei])
            for i in range(nRecs):
                ssi = si + i
                Slope[ssi] = ((((Slope1 - Slope0) / nRecs) * i) + Slope0)
            data[si:ei] = Slope[si:ei] * data[si:ei] + Offset
            flag[si:ei] = 10
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,float(c.missing_value))
            ds.series[ThisOne]['Flag'] = flag

def ApplyLinearDriftLocal(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. The slope is
        interpolated since the starting value at time 0 using a known 30-min
        increment.  Slope0, SlopeIncrement and Offset are defined in the control
        file.  This function applies to a dataset in which the start time in the
        control file is matched by dataset start time, but in which the end time
        in the control file extends beyond the dataset end.
        
        Usage qcts.ApplyLinearDriftLocal(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    if ThisOne not in ds.series.keys(): return
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'LocalDrift'):
        log.info('  Applying linear drift correction to '+ThisOne)
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(c.missing_value),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag']
        ldt = ds.series['DateTime']['Data']
        DriftList = cf['Variables'][ThisOne]['LocalDrift'].keys()
        for i in range(len(DriftList)):
            DriftItemList = ast.literal_eval(cf['Variables'][ThisOne]['LocalDrift'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(DriftItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(DriftItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.zeros(len(data))
            Slope0 = float(DriftItemList[2])
            SlopeIncrement = float(DriftItemList[3])
            Offset = float(DriftItemList[4])
            nRecs = len(Slope[si:ei])
            for i in range(nRecs):
                ssi = si + i
                Slope[ssi] = (SlopeIncrement * i) + Slope0
            data[si:ei] = Slope[si:ei] * data[si:ei] + Offset
            flag[si:ei] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,float(c.missing_value))
            ds.series[ThisOne]['Flag'] = flag

def AverageSeriesByElements(cf,ds,Av_out):
    """
        Calculates the average of multiple time series.  Multiple time series
        are entered and a single time series representing the average at each
        observational period is returned.
        
        Usage qcts.AverageSeriesByElements(cf,ds,Av_out)
        cf: control file object (must contain an entry for Av_out)
        ds: data structure
        Av_out: output variable to ds.  Example: 'Fg'
        Series_in: input variable series in ds.  Example: ['Fg_8cma','Fg_8cmb']
        """
    if Av_out not in cf['Variables'].keys(): return
    if Av_out in ds.averageserieslist: return
    srclist, standardname = qcutils.GetAverageSeriesKeys(cf,Av_out)
#    log.info(' Averaging series in '+str(srclist)+' into '+Av_out)
    log.info(' Averaging '+str(srclist)+'==>'+Av_out)
    
    nSeries = len(srclist)
    if nSeries==0:
        log.error('  AverageSeriesByElements: no input series specified for'+str(Av_out))
        return
    if nSeries==1:
        tmp_data = ds.series[srclist[0]]['Data'].copy()
        tmp_flag = ds.series[srclist[0]]['Flag'].copy()
        tmp_attr = ds.series[srclist[0]]['Attr'].copy()
        Av_data = numpy.ma.masked_where(tmp_data==float(c.missing_value),tmp_data)
        Mn_flag = tmp_flag
        SeriesNameString = srclist[0]
    else:
        tmp_data = ds.series[srclist[0]]['Data'].copy()
        tmp_flag = ds.series[srclist[0]]['Flag'].copy()

        index = numpy.where(numpy.mod(tmp_flag,10)==0)    # find the elements with flag = 0, 10, 20 etc
        tmp_flag[index] = 0                               # set them all to 0
        
        tmp_attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            SeriesNameString = SeriesNameString+', '+ThisOne
            tmp_data = numpy.vstack((tmp_data,ds.series[ThisOne]['Data'].copy()))
            tmp_flag = numpy.vstack((tmp_flag,ds.series[ThisOne]['Flag'].copy()))
        tmp_data = numpy.ma.masked_where(tmp_data==float(c.missing_value),tmp_data)
        Av_data = numpy.ma.average(tmp_data,axis=0)
        Mn_flag = numpy.min(tmp_flag,axis=0)
    ds.averageserieslist.append(Av_out)
    #attr = qcutils.MakeAttributeDictionary(long_name='Element-wise average of series '+SeriesNameString,
                                       #standard_name=standardname,units=ds.series[srclist[0]]['Attr']['units'])
    # this is a temporary fix, better to have a routine update the attr dictionary
    tmp_attr["long_name"] = tmp_attr["long_name"]+", element-wise average of series " + SeriesNameString
    qcutils.CreateSeries(ds,Av_out,Av_data,Flag=Mn_flag,Attr=tmp_attr)

def CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg'):
    """
        Calculate the average energy as Fn - G.
        
        Usage qcts.CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
        ds: data structure
        Fa_out: output available energy variable to ds.  Example: 'Fa'
        Fn_in: input net radiation in ds.  Example: 'Fn'
        Fg_in: input ground heat flux in ds.  Example: 'Fg'
        """
    log.info(' Calculating available energy from Fn and Fg')
    Fn,f,a = qcutils.GetSeriesasMA(ds,Fn_in)
    Fg,f,a = qcutils.GetSeriesasMA(ds,Fg_in)
    Fa_calc = Fn - Fg
    if Fa_out not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Available energy using '+Fn_in+','+Fg_in,units='W/m2')
        qcutils.CreateSeries(ds,Fa_out,Fa_calc,FList=[Fn_in,Fg_in],Attr=attr)
    else:
        Fa_exist,flag,attr = qcutils.GetSeriesasMA(ds,Fa_out)
        idx = numpy.where((numpy.ma.getmaskarray(Fa_exist)==True)&(numpy.ma.getmaskarray(Fa_calc)==False))[0]
        if len(idx)!=0:
            Fa_exist[idx] = Fa_calc[idx]
            flag[idx] = numpy.int32(20)
        qcutils.CreateSeries(ds,Fa_out,Fa_exist,Flag=flag,Attr=attr)

def CalculateFluxes(cf,ds):
    """
        Calculate the fluxes from the rotated covariances.
        
        Usage qcts.CalculateFluxes(ds)
        ds: data structure
        
        Pre-requisite: CoordRotation2D
        
        Accepts meteorological constants or variables
        """
    Ta,f,a = qcutils.GetSeriesasMA(ds,"Ta")
    ps,f,a = qcutils.GetSeriesasMA(ds,"ps")
    Ah,f,a = qcutils.GetSeriesasMA(ds,"Ah")
    rhom,f,a = qcutils.GetSeriesasMA(ds,"rhom")
    RhoCp,f,a = qcutils.GetSeriesasMA(ds,"RhoCp")
    Lv,f,a = qcutils.GetSeriesasMA(ds,"Lv")

    long_name = ''
    if 'Massman' in ds.globalattributes['Functions']:
        long_name = ' and frequency response corrected'
    
    log.info(" Calculating fluxes from covariances")
    if "wT" in ds.series.keys():
        ok_units = ["mC/s","Cm/s"]
        wT,flag,attr = qcutils.GetSeriesasMA(ds,"wT")
        if attr["units"] in ok_units:
            Fhv = RhoCp*wT
            attr["long_name"] = "Virtual heat flux, rotated to natural wind coordinates"+long_name
            attr["units"] = "W/m2"
            qcutils.CreateSeries(ds,"Fhv",Fhv,Flag=flag,Attr=attr)
        else:
            log.error(" CalculateFluxes: Incorrect units for wA, Fe not calculated")
    else:
        log.error("  CalculateFluxes: wT not found, Fh not calculated")
    if "wA" in ds.series.keys():
        wA,flag,attr = qcutils.GetSeriesasMA(ds,"wA")
        if attr["units"]=="g/m2/s":
            Fe = Lv*wA/float(1000)
            attr["long_name"] = "Latent heat flux, rotated to natural wind coordinates"+long_name
            attr["standard_name"] = "surface_upward_latent_heat_flux"
            attr["units"] = "W/m2"
            qcutils.CreateSeries(ds,"Fe",Fe,Flag=flag,Attr=attr)
        else:
            log.error(" CalculateFluxes: Incorrect units for wA, Fe not calculated")
    else:
        log.error("  CalculateFluxes: wA not found, Fe not calculated")
    if "wC" in ds.series.keys():
        wC,flag,attr = qcutils.GetSeriesasMA(ds,"wC")
        if attr["units"]=="mg/m2/s":
            Fc = wC
            attr["long_name"] = "CO2 flux, rotated to natural wind coordinates"+long_name
            attr["units"] = "mg/m2/s"
            qcutils.CreateSeries(ds,"Fc",Fc,Flag=flag,Attr=attr)
        else:
            log.error(" CalculateFluxes: Incorrect units for wC, Fc not calculated")
    else:
        log.error("  CalculateFluxes: wC not found, Fc not calculated")
    if "uw" in ds.series.keys():
        if "vw" in ds.series.keys():
            uw,f,a = qcutils.GetSeriesasMA(ds,"uw")
            vw,f,a = qcutils.GetSeriesasMA(ds,"vw")
            vs = uw*uw + vw*vw
            Fm = rhom*numpy.ma.sqrt(vs)
            us = numpy.ma.sqrt(numpy.ma.sqrt(vs))
            attr["long_name"] = "Momentum flux, rotated to natural wind coordinates"+long_name
            attr["units"] = "kg/m/s2"
            qcutils.CreateSeries(ds,"Fm",Fm,FList=["uw","vw"],Attr=attr)
            attr["long_name"] = "Friction velocity, rotated to natural wind coordinates"+long_name
            attr["units"] = "m/s"
            qcutils.CreateSeries(ds,"ustar",us,FList=["uw","vw"],Attr=attr)
        else:
            log.error("  CalculateFluxes: vw not found, Fm and ustar not calculated")
    else:
        log.error("  CalculateFluxes: uw not found, Fm and ustar not calculated")
    if 'CalculateFluxes' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CalculateFluxes'

def CalculateLongwave(ds,Fl_out,Fl_in,Tbody_in):
    """
        Calculate the longwave radiation given the raw thermopile output and the
        sensor body temperature.
        
        Usage qcts.CalculateLongwave(ds,Fl_out,Fl_in,Tbody_in)
        ds: data structure
        Fl_out: output longwave variable to ds.  Example: 'Flu'
        Fl_in: input longwave in ds.  Example: 'Flu_raw'
        Tbody_in: input sensor body temperature in ds.  Example: 'Tbody'
        """
    log.info(' Calculating longwave radiation')
    Fl_raw,f,a = qcutils.GetSeriesasMA(ds,Fl_in)
    Tbody,f,a = qcutils.GetSeriesasMA(ds,Tbody_in)
    Fl = Fl_raw + c.sb*(Tbody + 273.15)**4
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated longwave radiation using '+Fl_in+','+Tbody_in,units='W/m2')
    qcutils.CreateSeries(ds,Fl_out,Fl,FList=[Fl_in,Tbody_in],Attr=attr)

def CalculateHumidities(ds):
    """
    Purpose:
     Calculate any missing humidities from whatever is available.
     If absolute humidity (Ah) is available then;
      - calculate specific humidity (q) if it is not present
      - calculate relative humidity (RH) if it is not present
     If specific humidity (q) is available then;
      - calculate absolute humidity (Ah) if it is not present
      - calculate relative humidity (RH) if it is not present
     If reative humidity (RH) is available then;
      - calculate specific humidity (q) if it is not present
      - calculate relative humidity (RH) if it is not present
    Usage:
     qcts.CalculateHumidities(ds)
    Date:
     March 2015
    Author: PRI
    """
    if "Ah" not in ds.series.keys():
        if "q" in ds.series.keys():
            AbsoluteHumidityFromq(ds)    # calculate Ah from q
        elif "RH" in ds.series.keys():
            AbsoluteHumidityFromRH(ds)   # calculate Ah from RH
    if "q" not in ds.series.keys():
        if "Ah" in ds.series.keys():
            SpecificHumidityFromAh(ds)
        elif "RH" in ds.series.keys():
            SpecificHumidityFromRH(ds)
    if "RH" not in ds.series.keys():
        if "Ah" in ds.series.keys():
            RelativeHumidityFromAh(ds)
        elif "q" in ds.series.keys():
            RelativeHumidityFromq(ds)

def CalculateHumiditiesAfterGapFill(ds):
    """
    Purpose:
     Check to see which humidity quantities (Ah, RH or q) have been gap filled
     and, if necessary, calculate the other humidity quantities from the gap
     filled one.
    Usage:
     qcts.CalculateHumiditiesAfterGapFill(ds)
     where ds is a data structure
    Author: PRI
    Date: April 2015
    """
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
    # check to see if absolute humidity (Ah) was gap filled ...
    if "Ah" in gf_list:
        if "q" not in gf_list: SpecificHumidityFromAh(ds)
        if "RH" not in gf_list: RelativeHumidityFromAh(ds)
    # ... or was relative humidity (RH) gap filled ...
    elif "RH" in gf_list:
        if "Ah" not in gf_list: AbsoluteHumidityFromRH(ds)
        if "q" not in gf_list: SpecificHumidityFromRH(ds)
    # ... or was specific humidity (q) gap filled ...
    elif "q" in gf_list:
        if "Ah" not in gf_list: AbsoluteHumidityFromq(ds)
        if "RH" not in gf_list: RelativeHumidityFromq(ds)
    else:
        msg = "No humidities were gap filled!"
        log.warning(msg)

def AbsoluteHumidityFromRH(ds):
    """ Calculate absolute humidity from relative humidity. """
    log.info(' Calculating absolute humidity from relative humidity')
    Ta,Ta_flag,a = qcutils.GetSeriesasMA(ds,"Ta")
    RH,RH_flag,a = qcutils.GetSeriesasMA(ds,"RH")
    Ah_new_flag = qcutils.MergeQCFlag([Ta_flag,RH_flag])
    Ah_new = mf.absolutehumidityfromRH(Ta,RH)
    if "Ah" in ds.series.keys():
        Ah,Ah_flag,Ah_attr = qcutils.GetSeriesasMA(ds,"Ah")
        index = numpy.where(numpy.ma.getmaskarray(Ah)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Ah)==True)[0]
        Ah[index] = Ah_new[index]
        Ah_flag[index] = Ah_new_flag[index]
        Ah_attr["long_name"] = Ah_attr["long_name"]+", merged with Ah calculated from RH"
        qcutils.CreateSeries(ds,"Ah",Ah,Flag=Ah_flag,Attr=Ah_attr)
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity',units='g/m3',standard_name='mass_concentration_of_water_vapor_in_air')
        qcutils.CreateSeries(ds,'Ah',Ah_new,Flag=Ah_new_flag,Attr=attr)

def AbsoluteHumidityFromq(ds):
    """ Calculate absolute humidity from specific humidity. """
    log.info(' Calculating absolute humidity from specific humidity')
    Ta,Ta_flag,a = qcutils.GetSeriesasMA(ds,"Ta")
    ps,ps_flag,a = qcutils.GetSeriesasMA(ds,"ps")
    q,q_flag,a = qcutils.GetSeriesasMA(ds,"q")
    Ah_new_flag = qcutils.MergeQCFlag([Ta_flag,ps_flag,q_flag])
    RH = mf.RHfromspecifichumidity(q,Ta,ps)
    Ah_new = mf.absolutehumidityfromRH(Ta,RH)
    if "Ah" in ds.series.keys():
        Ah,Ah_flag,Ah_attr = qcutils.GetSeriesasMA(ds,"Ah")
        index = numpy.where(numpy.ma.getmaskarray(Ah)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(Ah)==True)[0]
        Ah[index] = Ah_new[index]
        Ah_flag[index] = Ah_new_flag[index]
        Ah_attr["long_name"] = Ah_attr["long_name"]+", merged with Ah calculated from q"
        qcutils.CreateSeries(ds,"Ah",Ah,Flag=Ah_flag,Attr=Ah_attr)
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity',units='g/m3',standard_name='mass_concentration_of_water_vapor_in_air')
        qcutils.CreateSeries(ds,"Ah",Ah_new,Flag=Ah_new_flag,Attr=attr)

def RelativeHumidityFromq(ds):
    """ Calculate relative humidity from specific humidity. """
    log.info(' Calculating relative humidity from specific humidity')
    Ta,Ta_flag,a = qcutils.GetSeriesasMA(ds,"Ta")
    ps,ps_flag,a = qcutils.GetSeriesasMA(ds,"ps")
    q,q_flag,a = qcutils.GetSeriesasMA(ds,"q")
    RH_new_flag = qcutils.MergeQCFlag([Ta_flag,ps_flag,q_flag])
    RH_new = mf.RHfromspecifichumidity(q,Ta,ps)
    if "RH" in ds.series.keys():
        RH,RH_flag,RH_attr = qcutils.GetSeriesasMA(ds,"RH")
        index = numpy.where(numpy.ma.getmaskarray(RH)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(RH)==True)[0]
        RH[index] = RH_new[index]
        RH_flag[index] = RH_new_flag[index]
        RH_attr["long_name"] = RH_attr["long_name"]+", merged with RH calculated from q"
        qcutils.CreateSeries(ds,"RH",RH,Flag=RH_flag,Attr=RH_attr)
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='relative_humidity')
        qcutils.CreateSeries(ds,'RH',RH_new,Flag=RH_new_flag,Attr=attr)
    
def RelativeHumidityFromAh(ds):
    """ Calculate relative humidity from absolute humidity. """
    log.info(' Calculating relative humidity from absolute humidity')
    Ta,Ta_flag,a = qcutils.GetSeriesasMA(ds,"Ta")
    Ah,Ah_flag,a = qcutils.GetSeriesasMA(ds,"Ah")
    RH_new_flag = qcutils.MergeQCFlag([Ta_flag,Ah_flag])
    RH_new = mf.RHfromabsolutehumidity(Ah,Ta)     # relative humidity in units of percent
    if "RH" in ds.series.keys():
        RH,RH_flag,RH_attr = qcutils.GetSeriesasMA(ds,"RH")
        index = numpy.where(numpy.ma.getmaskarray(RH)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(RH)==True)[0]
        RH[index] = RH_new[index]
        RH_flag[index] = RH_new_flag[index]
        RH_attr["long_name"] = RH_attr["long_name"]+", merged with RH calculated from Ah"
        qcutils.CreateSeries(ds,"RH",RH,Flag=RH_flag,Attr=RH_attr)
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='relative_humidity')
        qcutils.CreateSeries(ds,"RH",RH_new,Flag=RH_new_flag,Attr=attr)

def smooth(x,window_len=11,window='hanning'):
    """
    Purpose:
        Smooth the data using a window with requested size.
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
    Input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    Output:
        the smoothed signal
    Example:
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
    See also: 
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string
    Note:
        1) length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        2) odd values for window_len return output with different length from input
    Source:
        Lifted from scipy Cookbook (http://wiki.scipy.org/Cookbook/SignalSmooth)
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    y=numpy.convolve(w/w.sum(),s,mode='valid')
#    return y
    return y[(window_len/2-1):-(window_len/2)]

def SpecificHumidityFromAh(ds):
    """ Calculate specific humidity from absolute humidity. """
    log.info(' Calculating specific humidity from absolute humidity')
    Ta,Ta_flag,a = qcutils.GetSeriesasMA(ds,"Ta")
    ps,ps_flag,a = qcutils.GetSeriesasMA(ds,"ps")
    Ah,Ah_flag,a = qcutils.GetSeriesasMA(ds,"Ah")
    q_new_flag = qcutils.MergeQCFlag([Ta_flag,ps_flag,Ah_flag])
    RH = mf.RHfromabsolutehumidity(Ah,Ta)
    q_new = mf.specifichumidityfromRH(RH, Ta, ps)
    if "q" in ds.series.keys():
        q,q_flag,q_attr = qcutils.GetSeriesasMA(ds,"q")
        index = numpy.where(numpy.ma.getmaskarray(q)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(q)==True)[0]
        q[index] = q_new[index]
        q_flag[index] = q_new_flag[index]
        q_attr["long_name"] = q_attr["long_name"]+", merged with q calculated from Ah"
        qcutils.CreateSeries(ds,"q",q,Flag=q_flag,Attr=q_attr)
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity',units='kg/kg',standard_name='specific_humidity')
        qcutils.CreateSeries(ds,'q',q_new,Flag=q_new_flag,Attr=attr)

def SpecificHumidityFromRH(ds):
    """ Calculate specific humidity from relative humidity."""
    log.info(' Calculating specific humidity from relative humidity')
    Ta,Ta_flag,a = qcutils.GetSeriesasMA(ds,"Ta")
    ps,ps_flag,a = qcutils.GetSeriesasMA(ds,"ps")
    RH,RH_flag,a = qcutils.GetSeriesasMA(ds,"RH")
    q_new_flag = qcutils.MergeQCFlag([Ta_flag,ps_flag,RH_flag])
    q_new = mf.specifichumidityfromRH(RH,Ta,ps)   # specific humidity in units of kg/kg
    if "q" in ds.series.keys():
        q,q_flag,q_attr = qcutils.GetSeriesasMA(ds,"q")
        index = numpy.where(numpy.ma.getmaskarray(q)==True)[0]
        #index = numpy.ma.where(numpy.ma.getmaskarray(q)==True)[0]
        q[index] = q_new[index]
        q_flag[index] = q_new_flag[index]
        q_attr["long_name"] = q_attr["long_name"]+", merged with q calculated from RH"
        qcutils.CreateSeries(ds,"q",q,Flag=q_flag,Attr=q_attr)
    else:        
        attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity',units='kg/kg',standard_name='specific_humidity')
        qcutils.CreateSeries(ds,"q",q_new,Flag=q_new_flag,Attr=attr)

def CalculateMeteorologicalVariables(ds,Ta_name='Ta',Tv_name='Tv_CSAT',ps_name='ps',
                                     q_name="q",Ah_name='Ah',RH_name='RH',Cc_name='Cc'):
    """
        Add time series of meteorological variables based on fundamental
        relationships (Stull 1988)

        Usage qcts.CalculateMeteorologicalVariables(ds,Ta_name,ps_name,Ah_name)
        ds: data structure
        Ta_name: data series name for air temperature
        ps_name: data series name for pressure
        Ah_name: data series name for absolute humidity
        q_name : data series name for specific humidity
        RH_name: data series for relative humidity

        Variables added:
            rhom: density of moist air, mf.densitymoistair(Ta,ps,Ah)
            Lv: latent heat of vapourisation, mf.Lv(Ta)
            q: specific humidity, mf.specifichumidity(mr)
                where mr (mixing ratio) = mf.mixingratio(ps,vp)
            Cpm: specific heat of moist air, mf.specificheatmoistair(q)
            VPD: vapour pressure deficit, VPD = esat - e
        """
    for item in [Ta_name,ps_name,Ah_name,Cc_name,q_name]:
        if item not in ds.series.keys():
            msg = " CalculateMeteorologicalVariables: series "
            msg = msg + item + " not found, returning ..."
            log.warning(msg)
            return
    log.info(' Adding standard met variables to database')
    # get the required data series
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_name)
    # use Tv_CSAT if it is in the data structure, otherwise use Ta
    if Tv_name not in ds.series.keys(): Tv_name = Ta_name
    Tv,f,a = qcutils.GetSeriesasMA(ds,Tv_name)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_name)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_name)
    Cc,f,a = qcutils.GetSeriesasMA(ds,Cc_name)
    Cc_units = a["units"]
    q,f,a = qcutils.GetSeriesasMA(ds,q_name)
    # do the calculations
    e = mf.vapourpressure(Ah,Ta)                  # vapour pressure from absolute humidity and temperature
    esat = mf.es(Ta)                              # saturation vapour pressure
    rhod = mf.densitydryair(Ta,ps,e)              # partial density of dry air
    rhom = mf.densitymoistair(Ta,ps,e)            # density of moist air
    rhow = mf.densitywatervapour(Ta,e)            # partial density of water vapour
    Lv = mf.Lv(Ta)                                # latent heat of vapourisation
    mr = mf.mixingratio(ps,e)                     # mixing ratio
    mrsat = mf.mixingratio(ps,esat)               # saturation mixing ratio
    qsat = mf.specifichumidity(mrsat)             # saturation specific humidity from saturation mixing ratio
    Cpd = mf.specificheatcapacitydryair(Tv)
    Cpw = mf.specificheatcapacitywatervapour(Ta,Ah)
    RhoCp = mf.densitytimesspecificheat(rhow,Cpw,rhod,Cpd)
    Cpm = mf.specificheatmoistair(q)              # specific heat of moist air
    VPD = esat - e                                # vapour pressure deficit
    SHD = qsat - q                                # specific humidity deficit
    if Cc_units=="mg/m3":
        c_ppm = mf.co2_ppmfrommgpm3(Cc,Ta,ps)     # CO2 concentration in units of umol/mol
    else:
        c_ppm = Cc
    h_ppt = mf.h2o_mmolpmolfromgpm3(Ah,Ta,ps)     # H2O concentration in units of mmol/mol
    # write the meteorological series to the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure',units='kPa',standard_name='water_vapor_partial_pressure_in_air')
    qcutils.CreateSeries(ds,'e',e,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Saturation vapour pressure',units='kPa')
    qcutils.CreateSeries(ds,'esat',esat,FList=[Ta_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Density of dry air',units='kg/m3')
    qcutils.CreateSeries(ds,'rhod',rhod,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Density of moist air',units='kg/m3',standard_name='air_density')
    qcutils.CreateSeries(ds,'rhom',rhom,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Partial density of water vapour',units='kg/m3')
    qcutils.CreateSeries(ds,'rhow',rhow,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Latent heat of vapourisation',units='J/kg')
    qcutils.CreateSeries(ds,'Lv',Lv,FList=[Ta_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity of dry air',units='J/kg-K')
    qcutils.CreateSeries(ds,'Cpd',Cpd,FList=[Tv_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity of water vapour',units='J/kg-K')
    qcutils.CreateSeries(ds,'Cpw',Cpw,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity of moist air',units='J/kg-K')
    qcutils.CreateSeries(ds,'Cpm',Cpm,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Product of air density and specific heat capacity',units='J/m3-K')
    qcutils.CreateSeries(ds,'RhoCp',RhoCp,FList=[Ta_name,Tv_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure deficit',units='kPa',standard_name='water_vapor_saturation_deficit_in_air')
    qcutils.CreateSeries(ds,'VPD',VPD,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity deficit',units='kg/kg')
    qcutils.CreateSeries(ds,'SHD',SHD,FList=[Ta_name,Ah_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration',units='umol/mol')
    qcutils.CreateSeries(ds,'C_ppm',c_ppm,FList=[Cc_name,Ta_name,ps_name],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='H2O concentration',units='mmol/mol')
    qcutils.CreateSeries(ds,'H_ppt',h_ppt,FList=[Ah_name,Ta_name,ps_name],Attr=attr)
    if 'CalculateMetVars' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CalculateMetVars'

def CalculateNetRadiation(cf,ds,Fn_out='Fn',Fsd_in='Fsd',Fsu_in='Fsu',Fld_in='Fld',Flu_in='Flu'):
    """
    Purpose:
     Calculate the net radiation from the 4 components of the surface
     radiation budget.
    Usage:
     qcts.CalculateNetRadiation(cf,ds,Fn_out,Fsd_in,Fsu_in,Fld_in,Flu_in)
        cf: control file
        ds: data structure
        Fn_out: output net radiation variable to ds.  Example: 'Fn_KZ'
        Fsd_in: input downwelling solar radiation in ds.  Example: 'Fsd'
        Fsu_in: input upwelling solar radiation in ds.  Example: 'Fsu'
        Fld_in: input downwelling longwave radiation in ds.  Example: 'Fld'
        Flu_in: input upwelling longwave radiation in ds.  Example: 'Flu'
    Side effects:
     Creates a new series in the data structure containing the net radiation.
    Author: PRI
    Date: Sometime early on
    """
    log.info(' Calculating net radiation from 4 components')
    if Fsd_in in ds.series.keys() and Fsu_in in ds.series.keys() and Fld_in in ds.series.keys() and Flu_in in ds.series.keys():
        Fsd,f,a = qcutils.GetSeriesasMA(ds,Fsd_in)
        Fsu,f,a = qcutils.GetSeriesasMA(ds,Fsu_in)
        Fld,f,a = qcutils.GetSeriesasMA(ds,Fld_in)
        Flu,f,a = qcutils.GetSeriesasMA(ds,Flu_in)
        Fn_calc = (Fsd - Fsu) + (Fld - Flu)
        if Fn_out not in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation using '+Fsd_in+','+Fsu_in+','+Fld_in+','+Flu_in,
                                 standard_name='surface_net_downwawrd_radiative_flux',units='W/m2')
            qcutils.CreateSeries(ds,Fn_out,Fn_calc,FList=[Fsd_in,Fsu_in,Fld_in,Flu_in],Attr=attr)
        else:
            Fn_exist,flag,attr = qcutils.GetSeriesasMA(ds,Fn_out)
            idx = numpy.where((numpy.ma.getmaskarray(Fn_exist)==True)&(numpy.ma.getmaskarray(Fn_calc)==False))[0]
            #idx = numpy.ma.where((numpy.ma.getmaskarray(Fn_exist)==True)&(numpy.ma.getmaskarray(Fn_calc)==False))[0]
            if len(idx)!=0:
                Fn_exist[idx] = Fn_calc[idx]
                flag[idx] = numpy.int32(20)
            qcutils.CreateSeries(ds,Fn_out,Fn_exist,Flag=flag,Attr=attr)
    else:
        nRecs = int(ds.globalattributes['nc_nrecs'])
        Fn = numpy.array([c.missing_value]*nRecs,dtype=numpy.float64)
        flag = numpy.ones(nRecs,dtype=numpy.int32)
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation (one or more components missing)',
                             standard_name='surface_net_downwawrd_radiative_flux',units='W/m2')
        qcutils.CreateSeries(ds,Fn_out,Fn,Flag=flag,Attr=attr)

def CheckCovarianceUnits(ds):
    """
    Purpose:
    Usage:
    Author: PRI
    Date: September 2015
    """
    log.info(' Checking covariance units')
    co2_list = ["UxC","UyC","UzC"]
    h2o_list = ["UxA","UyA","UzA","UxH","UyH","UzH"]
    for item in co2_list:
        if item not in ds.series.keys(): continue
        data,flag,attr = qcutils.GetSeriesasMA(ds,item)
        if "umol" in attr["units"]:
            Ta,f,a = qcutils.GetSeriesasMA(ds,"Ta")
            ps,f,a = qcutils.GetSeriesasMA(ds,"ps")
            data = mf.co2_mgpm3fromppm(data,Ta,ps)
            attr["units"] = "mg/m2/s"
            qcutils.CreateSeries(ds,item,data,Flag=flag,Attr=attr)
    for item in h2o_list:
        if item not in ds.series.keys(): continue
        data,flag,attr = qcutils.GetSeriesasMA(ds,item)
        if "mmol" in attr["units"]:
            Ta,f,a = qcutils.GetSeriesasMA(ds,"Ta")
            ps,f,a = qcutils.GetSeriesasMA(ds,"ps")
            data = mf.h2o_gpm3frommmolpmol(data,Ta,ps)
            attr["units"] = "g/m2/s"
            if "H" in item: item = item.replace("H","A")
            qcutils.CreateSeries(ds,item,data,Flag=flag,Attr=attr)

def CoordRotation2D(cf,ds):
    """
        2D coordinate rotation to force v = w = 0.  Based on Lee et al, Chapter
        3 of Handbook of Micrometeorology.  This routine does not do the third
        rotation to force v'w' = 0.
        
        Usage qcts.CoordRotation2D(ds)
        ds: data structure
        """
    # get the raw wind velocity components
    Ux,f,a = qcutils.GetSeriesasMA(ds,'Ux')          # longitudinal component in CSAT coordinate system
    Uy,f,a = qcutils.GetSeriesasMA(ds,'Uy')          # lateral component in CSAT coordinate system
    Uz,f,a = qcutils.GetSeriesasMA(ds,'Uz')          # vertical component in CSAT coordinate system
    # get the raw covariances
    UxUz,f,UxUz_a = qcutils.GetSeriesasMA(ds,'UxUz') # covariance(Ux,Uz)
    UyUz,f,UyUz_a = qcutils.GetSeriesasMA(ds,'UyUz') # covariance(Uy,Uz)
    UxUy,f,a = qcutils.GetSeriesasMA(ds,'UxUy')      # covariance(Ux,Uy)
    UyUy,f,a = qcutils.GetSeriesasMA(ds,'UyUy')      # variance(Uy)
    UxUx,f,a = qcutils.GetSeriesasMA(ds,'UxUx')      # variance(Ux)
    UzUz,f,a = qcutils.GetSeriesasMA(ds,'UzUz')      # variance(Ux)
    UzC,f,UzC_a = qcutils.GetSeriesasMA(ds,'UzC')    # covariance(Uz,C)
    UzA,f,UzA_a = qcutils.GetSeriesasMA(ds,'UzA')    # covariance(Uz,A)
    UzT,f,UzT_a = qcutils.GetSeriesasMA(ds,'UzT')    # covariance(Uz,T)
    UxC,f,a = qcutils.GetSeriesasMA(ds,'UxC')        # covariance(Ux,C)
    UyC,f,a = qcutils.GetSeriesasMA(ds,'UyC')        # covariance(Uy,C)
    UxA,f,a = qcutils.GetSeriesasMA(ds,'UxA')        # covariance(Ux,A)
    UyA,f,a = qcutils.GetSeriesasMA(ds,'UyA')        # covariance(Ux,A)
    UxT,f,a = qcutils.GetSeriesasMA(ds,'UxT')        # covariance(Ux,T)
    UyT,f,a = qcutils.GetSeriesasMA(ds,'UyT')        # covariance(Uy,T)
    nRecs = int(ds.globalattributes['nc_nrecs'])     # number of records
    # get the instrument heights
    fm_height = "not defined"
    if "height" in UxUz_a: fm_height = UxUz_a["height"]
    fc_height = "not defined"
    if "height" in UzC_a: fc_height = UzC_a["height"]
    fe_height = "not defined"
    if "height" in UzA_a: fe_height = UzA_a["height"]
    fh_height = "not defined"
    if "height" in UzT_a: fh_height = UzT_a["height"]
    # apply 2D coordinate rotation unless otherwise specified in control file
    rotate = True
    if ('Options' in cf) and ('2DCoordRotation' in cf['Options'].keys()):
        if not cf['Options'].as_bool('2DCoordRotation'): rotate = False
    if rotate:
        log.info(' Applying 2D coordinate rotation (components and covariances)')
        # get the 2D and 3D wind speeds
        ws2d = numpy.ma.sqrt(Ux**2 + Uy**2)
        ws3d = numpy.ma.sqrt(Ux**2 + Uy**2 + Uz**2)
        # get the sine and cosine of the angles through which to rotate
        #  - first we rotate about the Uz axis by eta to get v = 0
        #  - then we rotate about the v axis by theta to get w = 0
        ce = Ux/ws2d          # cos(eta)
        se = Uy/ws2d          # sin(eta)
        ct = ws2d/ws3d        # cos(theta)
        st = Uz/ws3d          # sin(theta)
        # get the rotation angles
        theta = numpy.rad2deg(numpy.arctan2(st,ct))
        eta = numpy.rad2deg(numpy.arctan2(se,ce))
        # do the wind velocity components first
        u = Ux*ct*ce + Uy*ct*se + Uz*st           # longitudinal component in natural wind coordinates
        v = Uy*ce - Ux*se                         # lateral component in natural wind coordinates
        w = Uz*ct - Ux*st*ce - Uy*st*se           # vertical component in natural wind coordinates
        # do the variances
        uu = UxUx*ct**2*ce**2 + UyUy*ct**2*se**2 + UzUz*st**2 + 2*UxUy*ct**2*ce*se + 2*UxUz*ct*st*ce + 2*UyUz*ct*st*se
        vv = UyUy*ce**2 + UxUx*se**2 - 2*UxUy*ce*se
        ww = UzUz*ct**2 + UxUx*st**2*ce**2 + UyUy*st**2*se**2 - 2*UxUz*ct*st*ce - 2*UyUz*ct*st*se + 2*UxUy*st**2*ce*se
        # now do the scalar covariances
        wT = UzT*ct - UxT*st*ce - UyT*st*se       # covariance(w,T) in natural wind coordinate system
        wA = UzA*ct - UxA*st*ce - UyA*st*se       # covariance(w,A) in natural wind coordinate system
        wC = UzC*ct - UxC*st*ce - UyC*st*se       # covariance(w,C) in natural wind coordinate system
        # now do the momentum covariances
        # full equations, Wesely PhD thesis via James Cleverly and EddyPro
        uw = UxUz*ce*(ct*ct-st*st) - 2*UxUy*ct*st*ce*se + UyUz*se*(ct*ct-st*st) - \
             UxUx*ct*st*ce*ce - UyUy*ct*st*se*se + UzUz*ct*st # covariance(w,x) in natural wind coordinate system
        uv = UxUy*ct*(ce*ce-se*se) + UyUz*st*ce - UxUz*st*se - \
             UxUx*ct*ce*se + UyUy*ct*ce*se                    # covariance(x,y) in natural wind coordinate system
        vw = UyUz*ct*ce - UxUz*ct*se - UxUy*st*(ce*ce-se*se) + \
             UxUx*st*ce*se - UyUy*st*ce*se                    # covariance(w,y) in natural wind coordinate system
    else:
        log.info(' 2D coordinate rotation disabled, using unrotated components and covariances')
        # dummy series for rotation angles
        theta = numpy.zeros(nRecs)
        eta = numpy.zeros(nRecs)
        # unrotated wind components
        u = Ux           # unrotated x xomponent
        v = Uy           # unrotated y xomponent
        w = Uz           # unrotated z xomponent
        # unrotated covariances
        wT = UzT       # unrotated  wT covariance
        wA = UzA       # unrotated  wA covariance
        wC = UzC       # unrotated  wC covariance
        uw = UxUz      # unrotated  uw covariance
        vw = UyUz      # unrotated  vw covariance
        uv = UxUy      # unrotated  uv covariance
        # unrotated variances
        uu = UxUx      # unrotated  u variance
        vv = UyUy      # unrotated  v variance
        ww = UzUz      # unrotated  w variance
    # store the rotated quantities in the nc object
    # default behaviour of CreateSeries is to use the maximum value of the QC flag for any series specified in FList
    attr = qcutils.MakeAttributeDictionary(long_name='Horizontal rotation angle',units='deg',height=fm_height)
    qcutils.CreateSeries(ds,'eta',eta,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical rotation angle',units='deg',height=fm_height)
    qcutils.CreateSeries(ds,'theta',theta,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal component of wind-speed in natural wind coordinates',
                                           units='m/s',height=fm_height)
    qcutils.CreateSeries(ds,'u',u,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Lateral component of wind-speed in natural wind coordinates',
                                           units='m/s',height=fm_height)
    qcutils.CreateSeries(ds,'v',v,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical component of wind-speed in natural wind coordinates',
                                           units='m/s',height=fm_height)
    qcutils.CreateSeries(ds,'w',w,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic heat flux, rotated to natural wind coordinates',
                                           units='mC/s',height=fh_height)
    qcutils.CreateSeries(ds,'wT',wT,FList=['Ux','Uy','Uz','UxT','UyT','UzT'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic vapour flux, rotated to natural wind coordinates',
                                           units='g/m2/s',height=fe_height)
    qcutils.CreateSeries(ds,'wA',wA,FList=['Ux','Uy','Uz','UxA','UyA','UzA'],Attr=attr)
    #ReplaceRotatedCovariance(cf,ds,'wA','UzA')
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic CO2 flux, rotated to natural wind coordinates',
                                           units='mg/m2/s',height=fc_height)
    qcutils.CreateSeries(ds,'wC',wC,FList=['Ux','Uy','Uz','UxC','UyC','UzC'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux X component, corrected to natural wind coordinates',
                                           units='m2/s2',height=fm_height)
    qcutils.CreateSeries(ds,'uw',uw,FList=['Ux','Uy','Uz','UxUz','UxUx','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Horizontal streamwise-crosswind covariance, rotated to natural wind coordinates',
                                           units='m2/s2',height=fm_height)
    qcutils.CreateSeries(ds,'uv',uv,FList=['Ux','Uy','Uz','UxUz','UxUx','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux Y component, corrected to natural wind coordinates',
                                           units='m2/s2',height=fm_height)
    qcutils.CreateSeries(ds,'vw',vw,FList=['Ux','Uy','Uz','UyUz','UxUy','UyUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Variance of streamwise windspeed, rotated to natural wind coordinates',
                                           units='m2/s2',height=fm_height)
    qcutils.CreateSeries(ds,'uu',uu,FList=['Ux','Uy','Uz','UxUx','UxUy','UxUz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Variance of crossstream windspeed, rotated to natural wind coordinates',
                                           units='m2/s2',height=fm_height)
    qcutils.CreateSeries(ds,'vv',vv,FList=['Ux','Uy','Uz','UyUy','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Variance of vertical windspeed, rotated to natural wind coordinates',
                                           units='m2/s2',height=fm_height)
    qcutils.CreateSeries(ds,'ww',ww,FList=['Ux','Uy','Uz','UzUz','UxUz','UyUz'],Attr=attr)
    # if RotateFlag is set, force the QC flag value from the maximum of the FList series to 11
    #if qcutils.cfkeycheck(cf,Base='General',ThisOne='RotateFlag') and cf['General']['RotateFlag'] == 'True':
        #keys = ['eta','theta','u','v','w','wT','wA','wC','uw','vw']
        #for ThisOne in keys:
            #testseries,f = qcutils.GetSeriesasMA(ds,ThisOne)
            #mask = numpy.ma.getmask(testseries)
            #index = numpy.where(mask.astype(int)==1)
            #ds.series[ThisOne]['Flag'][index] = numpy.int32(11)
    if 'CoordRotation2D' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CoordRotation2D'
    if qcutils.cfoptionskeylogical(cf,Key='RelaxRotation'):
        RotatedSeriesList = ['wT','wA','wC','uw','vw']
        NonRotatedSeriesList = ['UzT','UzA','UzC','UxUz','UyUz']
        for ThisOne, ThatOne in zip(RotatedSeriesList,NonRotatedSeriesList):
            ReplaceWhereMissing(ds.series[ThisOne],ds.series[ThisOne],ds.series[ThatOne],FlagValue=21)
        if 'RelaxRotation' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxRotation'

def CalculateFcStorage(cf,ds,Fc_out='Fc_storage',CO2_in='Cc'):
    """
    Calculate CO2 flux storage term in the air column beneath the CO2 instrument.  This
    routine assumes the air column between the sensor and the surface is well mixed.
    
    Usage qcts.CalculateFcStorage(cf,ds,Fc_out,CO2_in)
    cf: control file object    
    ds: data structure
    Fc_out: series label of the CO2 flux storage term
    CO2_in: series label of the CO2 concentration
    
    Parameters loaded from control file:
        zms: measurement height from surface, m
    """
    if 'Fc_storage' not in ds.series.keys():
        if qcutils.cfkeycheck(cf,Base='General',ThisOne='zms'):
            log.info(' Calculating Fc storage (single height)')
            nRecs = int(ds.globalattributes['nc_nrecs'])
            ts = int(ds.globalattributes['time_step'])
            zms = float(cf['General']['zms'])
            # get the input data
            Cc,Cc_flag,Cc_attr = qcutils.GetSeriesasMA(ds,CO2_in,si=0,ei=-1)
            Ta,f,a = qcutils.GetSeriesasMA(ds,'Ta',si=0,ei=-1)
            ps,f,a = qcutils.GetSeriesasMA(ds,'ps',si=0,ei=-1)
            # check the CO2 concentration units
            # if the units are mg/m3, convert CO2 concentration to umol/mol before taking the difference
            if Cc_attr['units']=='mg/m3': Cc = mf.co2_ppmfrommgpm3(Cc,Ta,ps)
            # calculate the change in CO2 concentration between time steps, CO2 concentration in umol/mol.
            dc = numpy.ma.ediff1d(Cc,to_begin=0)
            # convert the CO2 concentration difference from umol/mol to mg/m3
            dc = mf.co2_mgpm3fromppm(dc,Ta,ps)
            # calculate the time step in seconds
            dt=86400*numpy.ediff1d(ds.series['xlDateTime']['Data'],to_begin=float(ts)/1440)    # time step in seconds from the Excel datetime values
            # calculate the CO2 flux based on storage below the measurement height
            Fc_storage = zms*dc/dt
            Fc_storage_units = 'mg/m2/s'
            descr = 'Fc infered from CO2 storage using single point CO2 measurement'
            # make the output series attribute dictionary
            attr_out = qcutils.MakeAttributeDictionary(long_name=descr,units=Fc_storage_units)
            # put the storage flux in the data structure
            qcutils.CreateSeries(ds,Fc_out,Fc_storage,FList=[CO2_in],Attr=attr_out)
        else:
            log.error('CalculateFcStorage: zms expected in General section of control file but not found')
    else:
        log.info('CalculateFcStorage: Fc_storage found in data structure, not calculated')

def CorrectFcForStorage(cf,ds,Fc_out='Fc',Fc_in='Fc',Fc_storage_in='Fc_storage'):
    """
    Correct CO2 flux for storage in the air column beneath the CO2 instrument.
    
    Usage qcts.CorrectFcForStorage(cf,ds,Fc_out,Fc_in,Fc_storage_in)
    cf: control file object    
    ds: data structure
    Fc_out: series label of the corrected CO2 flux
    Fc_in: series label of the input CO2 flux
    Fc_storage: series label of the CO2 flux storage term

    """
    if not qcutils.cfoptionskeylogical(cf,Key="ApplyFcStorage"): return
    if (Fc_in not in ds.series.keys()) or (Fc_storage_in not in ds.series.keys()):
        msg = "CorrectFcForStorage: Fc or Fc_storage not found, skipping ..."
        log.warning(msg)
        return
    log.info(" ***!!! Applying Fc storage term !!!***")
    Fc_raw,Fc_flag,Fc_attr = qcutils.GetSeriesasMA(ds,Fc_in)
    Fc_storage,Fc_storage_flag,Fc_storage_attr = qcutils.GetSeriesasMA(ds,Fc_storage_in)
    if Fc_attr["units"]!=Fc_storage_attr["units"]:
        log.error("CorrectFcForStorage: units of Fc do not match those of storage term, storage not applied")
        return
    log.info(" Applying storage correction to Fc")
    Fc = Fc_raw + Fc_storage
    if qcutils.cfoptionskeylogical(cf,Key="RelaxFcStorage"):
        idx=numpy.where(numpy.ma.getmaskarray(Fc)==True)[0]
        Fc[idx]=Fc_raw[idx]
        log.info(" Replaced corrected Fc with "+str(len(idx))+" raw values")
    Fc_attr["long_name"] = Fc_attr["long_name"] + ", uncorrected"
    qcutils.CreateSeries(ds,"Fc_raw",Fc_raw,Flag=Fc_flag,Attr=Fc_attr)
    Fc_attr["long_name"] = Fc_attr["long_name"].replace(", uncorrected",", corrected for storage using supplied storage term")
    qcutils.CreateSeries(ds,Fc_out,Fc,FList=[Fc_in,Fc_storage_in],Attr=Fc_attr)
    if "CorrectFcForStorage" not in ds.globalattributes["Functions"]:
        ds.globalattributes["Functions"] = ds.globalattributes["Functions"]+", CorrectFcForStorage"

def CorrectIndividualFgForStorage(cf,ds):
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CFgArgs'):
        List = cf['FunctionArgs']['CFgArgs'].keys()
        for i in range(len(List)):
            CFgArgs = ast.literal_eval(cf['FunctionArgs']['CFgArgs'][str(i)])
            CorrectFgForStorage(cf,ds,Fg_out=CFgArgs[0],Fg_in=CFgArgs[1],Ts_in=CFgArgs[2],Sws_in=CFgArgs[3])
        return

def CorrectFgForStorage(cf,ds,Fg_out='Fg',Fg_in='Fg',Ts_in='Ts',Sws_in='Sws'):
    """
        Correct ground heat flux for storage in the layer above the heat flux plate
        
        Usage qcts.CorrectFgForStorage(cf,ds,Fg_out,Fg_in,Ts_in,Sws_in)
        ds: data structure
        Fg_out: output soil heat flux variable to ds.  Example: 'Fg'
        Fg_in: input soil heat flux in ds.  Example: 'Fg_Av'
        Ts_in: input soil temperature in ds.  Example: 'Ts'
        
        Parameters loaded from control file:
            FgDepth: Depth of soil heat flux plates, m
            BulkDensity: soil bulk density, kg/m3
            OrganicContent: soil organic content, fraction
            SwsDefault: default value of soil moisture content used when no sensors present
        """
    # check to see if the user wants to skip the correction
    if not qcutils.cfoptionskeylogical(cf,Key="CorrectFgForStorage",default=True):
        log.info(' CorrectFgForStorage: storage correction disabled in control file')
        return
    # check to see if there is a [Soil] section in the control file
    if 'Soil' not in cf.keys():
        # if there isn't, check to see if the soil information is in the netCDF global attributes
        if "FgDepth" in ds.globalattributes.keys():
            # if it is, read it into the control file object so we can use it later
            cf["Soil"] = {}
            cf["Soil"]["FgDepth"] = ds.globalattributes["FgDepth"]
            cf["Soil"]["BulkDensity"] = ds.globalattributes["BulkDensity"]
            cf["Soil"]["OrganicContent"] = ds.globalattributes["OrganicContent"]
            cf["Soil"]["SwsDefault"] = ds.globalattributes["SwsDefault"]
        else:
            # tell the user if we can't find the information needed
            log.warning(' CorrectFgForStorage: [Soil] section not found in control file or global attributes, Fg not corrected')
            return
    if Fg_in not in ds.series.keys() or Ts_in not in ds.series.keys():
        log.warning(' CorrectFgForStorage: '+Fg_in+' or '+Ts_in+' not found in data structure, Fg not corrected')
        return
    log.info(' Correcting soil heat flux for storage')
    # put the contents of the soil section into the global attributes
    for item in cf["Soil"].keys(): ds.globalattributes[item] = cf["Soil"][item]
    d = max(0.0,min(0.5,float(cf['Soil']['FgDepth'])))
    bd = max(1200.0,min(2500.0,float(cf['Soil']['BulkDensity'])))
    oc = max(0.0,min(1.0,float(cf['Soil']['OrganicContent'])))
    mc = 1.0 - oc
    Sws_default = min(1.0,max(0.0,float(cf['Soil']['SwsDefault'])))
    # get the data
    nRecs = int(ds.globalattributes["nc_nrecs"])
    Fg,Fg_flag,Fg_attr = qcutils.GetSeriesasMA(ds,Fg_in)
    Ts,Ts_flag,Ts_attr = qcutils.GetSeriesasMA(ds,Ts_in)
    Sws,Sws_flag,Sws_attr = qcutils.GetSeriesasMA(ds,Sws_in)
    iom = numpy.where(numpy.mod(Sws_flag,10)!=0)[0]
    if len(iom)!=0:
        log.warning('  CorrectFgForStorage: Sws_default used for '+str(len(iom))+' values')
        Sws[iom] = Sws_default
        Sws_flag[iom] = numpy.int32(22)
    # get the soil temperature difference from time step to time step
    dTs = numpy.ma.zeros(nRecs)
    dTs[1:] = numpy.ma.diff(Ts)
    # write the temperature difference into the data structure so we can use its flag later
    dTs_flag = numpy.zeros(nRecs,dtype=numpy.int32)
    index = numpy.where(numpy.ma.getmaskarray(dTs)==True)[0]
    #index = numpy.ma.where(numpy.ma.getmaskarray(dTs)==True)[0]
    dTs_flag[index] = numpy.int32(1)
    attr = qcutils.MakeAttributeDictionary(long_name='Change in soil temperature',units='C')
    qcutils.CreateSeries(ds,"dTs",dTs,Flag=dTs_flag,Attr=attr)
    # get the time difference
    dt = numpy.ma.zeros(nRecs)
    dt[1:] = numpy.diff(date2num(ds.series['DateTime']['Data']))*float(86400)
    dt[0] = dt[1]
    # calculate the specific heat capacity of the soil
    Cs = mc*bd*c.Cd + oc*bd*c.Co + Sws*c.rho_water*c.Cw
    # calculate the soil heat storage
    S = Cs*(dTs/dt)*d
    # apply the storage term
    Fg_out_data = Fg + S
    # work out the QC flag
    Fg_out_flag = numpy.zeros(nRecs,dtype=numpy.int32)
    for item in [Fg_flag,Ts_flag,Sws_flag]:
        Fg_out_flag = numpy.maximum(Fg_out_flag,item)
    # trap and re-instate flag values of 1 (data missing at L1)
    for item in [Fg_flag,Ts_flag,Sws_flag]:
        index = numpy.where(item==numpy.int32(1))[0]
        Fg_out_flag[index] = numpy.int32(1)
    # put the corrected soil heat flux into the data structure
    attr= qcutils.MakeAttributeDictionary(long_name='Soil heat flux corrected for storage',units='W/m2',standard_name='downward_heat_flux_in_soil')
    qcutils.CreateSeries(ds,Fg_out,Fg_out_data,Flag=Fg_out_flag,Attr=attr)
    # save the input (uncorrected) soil heat flux series, this will be used if the correction is relaxed
    attr = qcutils.MakeAttributeDictionary(long_name='Soil heat flux uncorrected for storage',units='W/m2')
    qcutils.CreateSeries(ds,'Fg_Av',Fg,Flag=Fg_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Soil heat flux storage',units='W/m2')
    qcutils.CreateSeries(ds,'S',S,Flag=Fg_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity',units='J/m3/K')
    qcutils.CreateSeries(ds,'Cs',Cs,Flag=Fg_flag,Attr=attr)
    if 'CorrectFgForStorage' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', CorrectFgForStorage'
    if qcutils.cfoptionskeylogical(cf,Key='RelaxFgStorage'):
        ReplaceWhereMissing(ds.series['Fg'],ds.series['Fg'],ds.series['Fg_Av'],FlagValue=20)
        if 'RelaxFgStorage' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxFgStorage'

def CorrectSWC(cf,ds):
    """
        Correct soil moisture data using calibration curve developed from
        collected soil samples.  To avoid unrealistic or unphysical extreme
        values upon extrapolation, exponential and logarithmic using ln
        functions are applied to small and large values, respectively.
        Threshold values where one model replaces the other is determined where
        the functions cross.  The logarithmic curve is constrained at with a
        point at which the soil measurement = field porosity and the sensor
        measurement is maximised under saturation at field capacity.
        
        Usage qcts.CorrectSWC(cf,ds)
        cf: control file
        ds: data structure
        
        Parameters loaded from control file:
            SWCempList: list of raw CS616 variables
            SWCoutList: list of corrected CS616 variables
            SWCattr:  list of meta-data attributes for corrected CS616 variables
            SWC_a0: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            SWC_a1: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            SWC_b0: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            SWC_b1: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            SWC_t: threshold parameter for switching from exponential to logarithmic model
            TDRempList: list of raw CS610 variables
            TDRoutList: list of corrected CS610 variables
            TDRattr:  list of meta-data attributes for corrected CS610 variables
            TDRlinList: list of deep TDR probes requiring post-hoc linear correction to match empirical samples
            TDR_a0: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            TDR_a1: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            TDR_b0: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            TDR_b1: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            TDR_t: threshold parameter for switching from exponential to logarithmic model
        """
    if not qcutils.cfoptionskeylogical(cf,Key='CorrectSWC'): return
    log.info(' Correcting soil moisture data ...')
    SWCempList = ast.literal_eval(cf['Soil']['empSWCin'])
    SWCoutList = ast.literal_eval(cf['Soil']['empSWCout'])
    SWCattr = ast.literal_eval(cf['Soil']['SWCattr'])
    if cf['Soil']['TDR']=='Yes':
        TDRempList = ast.literal_eval(cf['Soil']['empTDRin'])
        TDRoutList = ast.literal_eval(cf['Soil']['empTDRout'])
        TDRlinList = ast.literal_eval(cf['Soil']['linTDRin'])
        TDRattr = ast.literal_eval(cf['Soil']['TDRattr'])
        TDR_a0 = float(cf['Soil']['TDR_a0'])
        TDR_a1 = float(cf['Soil']['TDR_a1'])
        TDR_b0 = float(cf['Soil']['TDR_b0'])
        TDR_b1 = float(cf['Soil']['TDR_b1'])
        TDR_t = float(cf['Soil']['TDR_t'])
    SWC_a0 = float(cf['Soil']['SWC_a0'])
    SWC_a1 = float(cf['Soil']['SWC_a1'])
    SWC_b0 = float(cf['Soil']['SWC_b0'])
    SWC_b1 = float(cf['Soil']['SWC_b1'])
    SWC_t = float(cf['Soil']['SWC_t'])
    
    for i in range(len(SWCempList)):
        log.info('  Applying empirical correction to '+SWCempList[i])
        invar = SWCempList[i]
        outvar = SWCoutList[i]
        attr = SWCattr[i]
        Sws,f,a = qcutils.GetSeriesasMA(ds,invar)
        
        nRecs = len(Sws)
        
        Sws_out = numpy.ma.empty(nRecs,float)
        Sws_out.fill(c.missing_value)
        Sws_out.mask = numpy.ma.empty(nRecs,bool)
        Sws_out.mask.fill(True)
        
        index_high = numpy.ma.where((Sws.mask == False) & (Sws > SWC_t))[0]
        index_low = numpy.ma.where((Sws.mask == False) & (Sws < SWC_t))[0]
        
        Sws_out[index_low] = SWC_b0 * numpy.exp(SWC_b1 * Sws[index_low])
        Sws_out[index_high] = (SWC_a1 * numpy.log(Sws[index_high])) + SWC_a0
        
        attr = qcutils.MakeAttributeDictionary(long_name=attr,units='cm3 water/cm3 soil',standard_name='soil_moisture_content')
        qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)
    if cf['Soil']['TDR']=='Yes':
        for i in range(len(TDRempList)):
            log.info('  Applying empirical correction to '+TDRempList[i])
            invar = TDRempList[i]
            outvar = TDRoutList[i]
            attr = TDRattr[i]
            Sws,f,a = qcutils.GetSeriesasMA(ds,invar)
            
            nRecs = len(Sws)
            
            Sws_out = numpy.ma.empty(nRecs,float)
            Sws_out.fill(c.missing_value)
            Sws_out.mask = numpy.ma.empty(nRecs,bool)
            Sws_out.mask.fill(True)
            
            index_high = numpy.ma.where((Sws.mask == False) & (Sws > TDR_t))[0]
            index_low = numpy.ma.where((Sws.mask == False) & (Sws < TDR_t))[0]
            
            Sws_out[index_low] = TDR_b0 * numpy.exp(TDR_b1 * Sws[index_low])
            Sws_out[index_high] = (TDR_a1 * numpy.log(Sws[index_high])) + TDR_a0
            
            attr = qcutils.MakeAttributeDictionary(long_name=attr,units='cm3 water/cm3 soil',standard_name='soil_moisture_content')
            qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)

def CorrectWindDirection(cf,ds,Wd_in):
    """
        Correct wind direction for mis-aligned sensor direction.
        
        Usage qcts.CorrectWindDirection(cf,ds,Wd_in)
        cf: control file
        ds: data structure
        Wd_in: input/output wind direction variable in ds.  Example: 'Wd_CSAT'
        """
    log.info(' Correcting wind direction')
    Wd,f,a = qcutils.GetSeriesasMA(ds,Wd_in)
    ldt = ds.series['DateTime']['Data']
    KeyList = cf['Variables'][Wd_in]['Correction'].keys()
    for i in range(len(KeyList)):
        ItemList = ast.literal_eval(cf['Variables'][Wd_in]['Correction'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ItemList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ItemList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        Correction = float(ItemList[2])
        Wd[si:ei] = Wd[si:ei] + Correction
    Wd = numpy.mod(Wd,float(360))
    ds.series[Wd_in]['Data'] = numpy.ma.filled(Wd,float(c.missing_value))

def LowPassFilterSws(cf,ds,Sws_out='Sws_LP',Sws_in='Sws',npoles=5,co_ny=0.05):
    '''
    Create a series of daily averaged soil moisture data and then interpolate this
    back on to the time step of the data.  This result is a time series of soil
    moisture data that is less noisy than the data at the original time step but
    still resolves day-to-day changes and seasonal trends.
    '''
    b,a = butter(npoles,co_ny)
    Sws,f,a = qcutils.GetSeries(ds,Sws_in)
    Sws_LP = filtfilt(b,a,Sws)
    attr = qcutils.MakeAttributeDictionary(long_name=attr,units='cm3 water/cm3 soil',standard_name='soil_moisture_content')
    qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)
    

def do_attributes(cf,ds):
    """
        Import attriubes in xl2nc control file to netCDF dataset.  Included
        global and variable attributes.  Also attach flag definitions to global
        meta-data for reference.
        
        Usage qcts.do_attributes(cf,ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting the attributes given in control file')
    if 'Global' in cf.keys():
        for gattr in cf['Global'].keys():
            ds.globalattributes[gattr] = cf['Global'][gattr]
        ds.globalattributes['Flag00'] = 'Good data'
        ds.globalattributes['Flag10'] = 'Corrections: Apply Linear'
        ds.globalattributes['Flag20'] = 'GapFilling: Driver gap filled using ACCESS'
        ds.globalattributes['Flag30'] = 'GapFilling: Flux gap filled by ANN (SOLO)'
        ds.globalattributes['Flag40'] = 'GapFilling: Gap filled by climatology'
        ds.globalattributes['Flag50'] = 'GapFilling: Gap filled by interpolation'
        ds.globalattributes['Flag60'] = 'GapFilling: Flux gap filled using ratios'
        ds.globalattributes['Flag01'] = 'QA/QC: Missing value in L1 dataset'
        ds.globalattributes['Flag02'] = 'QA/QC: L2 Range Check'
        ds.globalattributes['Flag03'] = 'QA/QC: CSAT Diagnostic'
        ds.globalattributes['Flag04'] = 'QA/QC: LI7500 Diagnostic'
        ds.globalattributes['Flag05'] = 'QA/QC: L2 Diurnal SD Check'
        ds.globalattributes['Flag06'] = 'QA/QC: Excluded Dates'
        ds.globalattributes['Flag07'] = 'QA/QC: Excluded Hours'
        ds.globalattributes['Flag08'] = 'QA/QC: Missing value found with QC flag = 0'
        ds.globalattributes['Flag11'] = 'Corrections/Combinations: Coordinate Rotation (Ux, Uy, Uz, UxT, UyT, UzT, UxA, UyA, UzA, UxC, UyC, UzC, UxUz, UxUx, UxUy, UyUz, UxUy, UyUy)'
        ds.globalattributes['Flag12'] = 'Corrections/Combinations: Massman Frequency Attenuation Correction (Coord Rotation, Tv_CSAT, Ah_HMP, ps)'
        ds.globalattributes['Flag13'] = 'Corrections/Combinations: Virtual to Actual Fh (Coord Rotation, Massman, Ta_HMP)'
        ds.globalattributes['Flag14'] = 'Corrections/Combinations: WPL correction for flux effects on density measurements (Coord Rotation, Massman, Fhv to Fh, Cc_7500_Av)'
        ds.globalattributes['Flag15'] = 'Corrections/Combinations: Ta from Tv'
        ds.globalattributes['Flag16'] = 'Corrections/Combinations: L3 Range Check'
        ds.globalattributes['Flag17'] = 'Corrections/Combinations: L3 Diurnal SD Check'
        ds.globalattributes['Flag18'] = 'Corrections/Combinations: u* filter'
        ds.globalattributes['Flag19'] = 'Corrections/Combinations: Gap coordination'
        ds.globalattributes['Flag21'] = 'GapFilling: Used non-rotated covariance'
        ds.globalattributes['Flag31'] = 'GapFilling: Flux gap not filled by ANN'
        ds.globalattributes['Flag38'] = 'GapFilling: L4 Range Check'
        ds.globalattributes['Flag39'] = 'GapFilling: L4 Diurnal SD Check'
        # the following flags are used by James Cleverly's version but not
        # by the standard OzFlux version.
        #ds.globalattributes['Flag51'] = 'albedo: bad Fsd < threshold (290 W/m2 default) only if bad time flag (31) not set'
        #ds.globalattributes['Flag52'] = 'albedo: bad time flag (not midday 10.00 to 14.00)'
        #ds.globalattributes['Flag61'] = 'Penman-Monteith: bad rst (rst < 0) only if bad Uavg (35), bad Fe (33) and bad Fsd (34) flags not set'
        #ds.globalattributes['Flag62'] = 'Penman-Monteith: bad Fe < threshold (0 W/m2 default) only if bad Fsd (34) flag not set'
        #ds.globalattributes['Flag63'] = 'Penman-Monteith: bad Fsd < threshold (10 W/m2 default)'
        #ds.globalattributes['Flag64'] = 'Penman-Monteith: Uavg == 0 (undefined aerodynamic resistance under calm conditions) only if bad Fe (33) and bad Fsd (34) flags not set'
        #ds.globalattributes['Flag70'] = 'Partitioning Night: Re computed from exponential temperature response curves'
        #ds.globalattributes['Flag80'] = 'Partitioning Day: GPP/Re computed from light-response curves, GPP = Re - Fc'
        #ds.globalattributes['Flag81'] = 'Partitioning Day: GPP night mask'
        #ds.globalattributes['Flag82'] = 'Partitioning Day: Fc > Re, GPP = 0, Re = Fc'
    for ThisOne in ds.series.keys():
        if ThisOne in cf['Variables']:
            if 'Attr' in cf['Variables'][ThisOne].keys():
                ds.series[ThisOne]['Attr'] = {}
                for attr in cf['Variables'][ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = cf['Variables'][ThisOne]['Attr'][attr]
                if "missing_value" not in ds.series[ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr']["missing_value"] = numpy.int32(c.missing_value)

def DoFunctions(cf,ds):
    """
    Purpose:
     Evaluate functions used in the L1 control file.
    Usage:
    Author: PRI
    Date: September 2015
    """
    implemented_functions = [name for name,data in inspect.getmembers(qcfunc,inspect.isfunction)]
    for var in cf["Variables"].keys():
        if "Function" not in cf["Variables"][var].keys(): continue
        if "func" not in cf["Variables"][var]["Function"].keys():
            msg = " DoFunctions: 'func' keyword not found in [Functions] for "+var
            log.error(msg)
            continue
        function_string = cf["Variables"][var]["Function"]["func"]
        function_name = function_string.split("(")[0]
        if function_name not in implemented_functions:
            msg = " DoFunctions: Requested function "+function_name+" not imlemented, skipping ..."
            log.error(msg)
            continue
        function_args = function_string.split("(")[1].replace(")","").split(",")
        result = getattr(qcfunc,function_name)(ds,var,*function_args)
        msg = " Completed function for "+var
        log.info(msg)

def CalculateStandardDeviations(cf,ds):
    log.info(' Getting variances from standard deviations & vice versa')
    if 'AhAh' in ds.series.keys() and 'Ah_7500_Sd' not in ds.series.keys():
        AhAh,flag,attr = qcutils.GetSeriesasMA(ds,'AhAh')
        Ah_7500_Sd = numpy.ma.sqrt(AhAh)
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity from Li-7500, standard deviation',units='g/m3')
        qcutils.CreateSeries(ds,'Ah_7500_Sd',Ah_7500_Sd,Flag=flag,Attr=attr)
    if 'Ah_7500_Sd' in ds.series.keys() and 'AhAh' not in ds.series.keys():
        Ah_7500_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Ah_7500_Sd')
        AhAh = Ah_7500_Sd*Ah_7500_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity from Li-7500, variance',units='(g/m3)2')
        qcutils.CreateSeries(ds,'AhAh',AhAh,Flag=flag,Attr=attr)
    if 'CcCc' in ds.series.keys() and 'Cc_7500_Sd' not in ds.series.keys():
        CcCc,flag,attr = qcutils.GetSeriesasMA(ds,'CcCc')
        Cc_7500_Sd = numpy.ma.sqrt(CcCc)
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration from Li-7500, standard deviation',units='mg/m3')
        qcutils.CreateSeries(ds,'Cc_7500_Sd',Cc_7500_Sd,Flag=flag,Attr=attr)
    if 'Cc_7500_Sd' in ds.series.keys() and 'CcCc' not in ds.series.keys():
        Cc_7500_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Cc_7500_Sd')
        CcCc = Cc_7500_Sd*Cc_7500_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration from Li-7500, variance',units='(mg/m3)2')
        qcutils.CreateSeries(ds,'CcCc',CcCc,Flag=flag,Attr=attr)
    if 'Ux_Sd' in ds.series.keys() and 'UxUx' not in ds.series.keys():
        Ux_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Ux_Sd')
        UxUx = Ux_Sd*Ux_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UxUx',UxUx,Flag=flag,Attr=attr)
    if 'UxUx' in ds.series.keys() and 'Ux_Sd' not in ds.series.keys():
        UxUx,flag,attr = qcutils.GetSeriesasMA(ds,'UxUx')
        Ux_Sd = numpy.ma.sqrt(UxUx)
        attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal velocity component from CSAT, standard deviation',units='m/s')
        qcutils.CreateSeries(ds,'Ux_Sd',Ux_Sd,Flag=flag,Attr=attr)
    if 'Uy_Sd' in ds.series.keys() and 'UyUy' not in ds.series.keys():
        Uy_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Uy_Sd')
        UyUy = Uy_Sd*Uy_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Lateral velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UyUy',UyUy,Flag=flag,Attr=attr)
    if 'UyUy' in ds.series.keys() and 'Uy_Sd' not in ds.series.keys():
        UyUy,flag,attr = qcutils.GetSeriesasMA(ds,'UyUy')
        Uy_Sd = numpy.ma.sqrt(UyUy)
        attr = qcutils.MakeAttributeDictionary(long_name='Lateral velocity component from CSAT, standard deviation',units='m/s')
        qcutils.CreateSeries(ds,'Uy_Sd',Uy_Sd,Flag=flag,Attr=attr)
    if 'Uz_Sd' in ds.series.keys() and 'UzUz' not in ds.series.keys():
        Uz_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Uz_Sd')
        UzUz = Uz_Sd*Uz_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Vertical velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UzUz',UzUz,Flag=flag,Attr=attr)
    if 'UzUz' in ds.series.keys() and 'Uz_Sd' not in ds.series.keys():
        UzUz,flag,attr = qcutils.GetSeriesasMA(ds,'UzUz')
        Uz_Sd = numpy.ma.sqrt(UzUz)
        attr = qcutils.MakeAttributeDictionary(long_name='Vertical velocity component from CSAT, standard deviation',units='m/s')
        qcutils.CreateSeries(ds,'Uz_Sd',Uz_Sd,Flag=flag,Attr=attr)

def do_solo(cf,ds4,Fc_in='Fc',Fe_in='Fe',Fh_in='Fh',Fc_out='Fc',Fe_out='Fe',Fh_out='Fh'):
    ''' duplicate gapfilled fluxes for graphing comparison'''
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='SOLOvars'):
        invars = ast.literal_eval(cf['FunctionArgs']['SOLOvars'])
        Fc_in = invars[0]
        Fe_in = invars[1]
        Fh_in = invars[2]
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='SOLOplot'):
        outvars = ast.literal_eval(cf['FunctionArgs']['SOLOplot'])
        Fc_out = outvars[0]
        Fe_out = outvars[1]
        Fh_out = outvars[2]
    # add relevant meteorological values to L3 data
    log.info(' Adding standard met variables to database')
    CalculateMeteorologicalVariables(ds4)
    ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', CalculateMetVars'
    if Fe_in in ds4.series.keys():
        Fe,flag,attr = qcutils.GetSeriesasMA(ds4,Fe_in)
        attr = qcutils.MakeAttributeDictionary(long_name='ANN gapfilled Latent Heat Flux',units='W/m2',standard_name='surface_upward_latent_heat_flux')
        qcutils.CreateSeries(ds4,Fe_out,Fe,Flag=flag,Attr=attr)
    if Fc_in in ds4.series.keys():
        Fc,flag,attr = qcutils.GetSeriesasMA(ds4,Fc_in)
        attr = qcutils.MakeAttributeDictionary(long_name='ANN gapfilled Carbon Flux',units='mg/m2/s')
        qcutils.CreateSeries(ds4,Fc_out,Fc,Flag=flag,Attr=attr)
    if Fh_in in ds4.series.keys():
        Fh,flag,attr = qcutils.GetSeriesasMA(ds4,Fh_in)
        attr = qcutils.MakeAttributeDictionary(long_name='ANN gapfilled Sensible Heat Flux',units='W/m2',standard_name='surface_upward_sensible_heat_flux')
        qcutils.CreateSeries(ds4,Fh_out,Fh,Flag=flag,Attr=attr)

def Fc_WPL(cf,ds,Fc_wpl_out='Fc',Fc_raw_in='Fc',Fh_in='Fh',Fe_in='Fe',Ta_in='Ta',Ah_in='Ah',Cc_in='Cc',ps_in='ps'):
    """
        Apply Webb, Pearman and Leuning correction to carbon flux.  This
        correction is necessary to account for flux effects on density
        measurements.  Original formulation: Campbell Scientific
        
        Usage qcts.Fc_WPL(ds,Fc_wpl_out,Fc_raw_in,Fh_in,Fe_raw_in,Ta_in,Ah_in,Cc_in,ps_in)
        ds: data structure
        Fc_wpl_out: output corrected carbon flux variable to ds.  Example: 'Fc'
        Fc_raw_in: input carbon flux in ds.  Example: 'Fc'
        Fh_in: input sensible heat flux in ds.  Example: 'Fh'
        Fe_raw_in: input uncorrected latent heat flux in ds.  Example: 'Fe_raw'
        Ta_in: input air temperature in ds.  Example: 'Ta'
        Ah_in: input absolute humidity in ds.  Example: 'Ah'
        Cc_in: input co2 density in ds.  Example: 'Cc'
        ps_in: input atmospheric pressure in ds.  Example: 'ps'
        
        Used for fluxes that are raw or rotated.
        
        Pre-requisite: CalculateFluxes, CalculateFluxes_Unrotated or CalculateFluxesRM
        Pre-requisite: FhvtoFh
        Pre-requisite: Fe_WPL
        
        Accepts meteorological constants or variables
        """
    if 'DisableFcWPL' in cf['Options'] and cf['Options'].as_bool('DisableFcWPL'):
        log.warning(" WPL correction for Fc disabled in control file")
        return
    log.info(' Applying WPL correction to Fc')
    Fc_raw,Fc_raw_flag,Fc_raw_attr = qcutils.GetSeriesasMA(ds,Fc_raw_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    Fe,f,a = qcutils.GetSeriesasMA(ds,Fe_in)
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    TaK = Ta+c.C2K                                # air temperature from C to K
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    Ah = Ah*c.g2kg                                # absolute humidity from g/m3 to kg/m3
    Cc,f,a = qcutils.GetSeriesasMA(ds,Cc_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    rhod,f,a = qcutils.GetSeriesasMA(ds,'rhod')
    RhoCp,f,a = qcutils.GetSeriesasMA(ds,'RhoCp')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    sigma = Ah/rhod
    co2_wpl_Fe = (c.mu/(1+c.mu*sigma))*(Cc/rhod)*(Fe/Lv)
    co2_wpl_Fh = (Cc/TaK)*(Fh/RhoCp)
    Fc_wpl_data = Fc_raw+co2_wpl_Fe+co2_wpl_Fh
    Fc_wpl_flag = numpy.zeros(len(Fc_wpl_data))
    index = numpy.where(numpy.ma.getmaskarray(Fc_wpl_data)==True)[0]
    Fc_wpl_flag[index] = numpy.int32(14)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL corrected Fc',units='mg/m2/s')
    if "height" in Fc_raw_attr: attr["height"] = Fc_raw_attr["height"]
    qcutils.CreateSeries(ds,Fc_wpl_out,Fc_wpl_data,Flag=Fc_wpl_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL correction to Fc due to Fe',units='mg/m2/s')
    if "height" in Fc_raw_attr: attr["height"] = Fc_raw_attr["height"]
    qcutils.CreateSeries(ds,'co2_wpl_Fe',co2_wpl_Fe,Flag=Fc_wpl_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL correction to Fc due to Fh',units='mg/m2/s')
    if "height" in Fc_raw_attr: attr["height"] = Fc_raw_attr["height"]
    qcutils.CreateSeries(ds,'co2_wpl_Fh',co2_wpl_Fh,Flag=Fc_wpl_flag,Attr=attr)

def Fe_WPL(cf,ds,Fe_wpl_out='Fe',Fe_raw_in='Fe',Fh_in='Fh',Ta_in='Ta',Ah_in='Ah',ps_in='ps'):
    """
        Apply Webb, Pearman and Leuning correction to vapour flux.  This
        correction is necessary to account for flux effects on density
        measurements.  Original formulation: Campbell Scientific
        
        Usage qcts.Fe_WPL(ds,Fe_wpl_out,Fe_raw_in,Fh_in,Ta_in,Ah_in,ps_in)
        ds: data structure
        Fe_wpl_out: output corrected water vapour flux variable to ds.  Example: 'Fe'
        Fe_raw_in: input water vapour flux in ds.  Example: 'Fe'
        Fh_in: input sensible heat flux in ds.  Example: 'Fh'
        Ta_in: input air temperature in ds.  Example: 'Ta'
        Ah_in: input absolute humidity in ds.  Example: 'Ah'
        ps_in: input atmospheric pressure in ds.  Example: 'ps'
        
        Used for fluxes that are raw or rotated.
        
        Pre-requisite: CalculateFluxes, CalculateFluxes_Unrotated or CalculateFluxesRM
        Pre-requisite: FhvtoFh
        
        Accepts meteorological constants or variables
        """
    if 'DisableFeWPL' in cf['Options'] and cf['Options'].as_bool('DisableFeWPL'):
        log.warning(" WPL correction for Fe disabled in control file")
        return
    log.info(' Applying WPL correction to Fe')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='EWPL'):
        Eargs = ast.literal_eval(cf['FunctionArgs']['EWPL'])
        Fe_wpl_out = Eargs[0]
        Fe_raw_in = Eargs[1]
        Fh_in = Eargs[2]
        Ta_in = Eargs[3]
        Ah_in = Eargs[4]
        ps_in = Eargs[5]
    Fe_raw,Fe_raw_flag,Fe_raw_attr = qcutils.GetSeriesasMA(ds,Fe_raw_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    TaK = Ta + c.C2K                              # air temperature from C to K
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    rhod,f,a = qcutils.GetSeriesasMA(ds,'rhod')     # density dry air
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')     # density moist air
    RhoCp,f,a = qcutils.GetSeriesasMA(ds,'RhoCp')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    Ah = Ah*c.g2kg                                # absolute humidity from g/m3 to kg/m3
    sigma = Ah/rhod
    h2o_wpl_Fe = c.mu*sigma*Fe_raw
    h2o_wpl_Fh = (1+c.mu*sigma)*Ah*Lv*(Fh/RhoCp)/TaK
    Fe_wpl_data = Fe_raw+h2o_wpl_Fe+h2o_wpl_Fh
    Fe_wpl_flag = numpy.zeros(len(Fe_wpl_data))
    mask = numpy.ma.getmask(Fe_wpl_data)
    index = numpy.where(numpy.ma.getmaskarray(Fe_wpl_data)==True)[0]
    Fe_wpl_flag[index] = numpy.int32(14)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL corrected Fe',
                                           standard_name='surface_upward_latent_heat_flux',
                                           units='W/m2')
    if "height" in Fe_raw_attr: attr["height"] = Fe_raw_attr["height"]
    qcutils.CreateSeries(ds,Fe_wpl_out,Fe_wpl_data,Flag=Fe_wpl_flag,Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Fe (uncorrected for WPL)',units='W/m2')
    if "height" in Fe_raw_attr: attr["height"] = Fe_raw_attr["height"]
    qcutils.CreateSeries(ds,'Fe_raw',Fe_raw,Flag=Fe_raw_flag,Attr=attr)
    if qcutils.cfoptionskeylogical(cf,Key='RelaxFeWPL'):
        ReplaceWhereMissing(ds.series['Fe'],ds.series['Fe'],ds.series['Fe_raw'],FlagValue=20)
        if 'RelaxFeWPL' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxFeWPL'

def FhvtoFh(cf,ds,Fh_out='Fh',Fhv_in='Fhv',Tv_in='Tv_CSAT',q_in='q',wA_in='wA',wT_in='wT'):
    '''
    Convert the virtual heat flux to the sensible heat flux.
    USEAGE:
     qcts.FhvtoFh_EP(cf,ds,Fhv_in='Fhv',RhoCp_in='RhoCp',Tv_in='Tv_CSAT',wA_in='wA',rhom_in='rhom',q_in='q',wT_in='wT')
    INPUT:
     All inputs are read from the data structure.
      Fhv_in   - label of the virtual heat flux series, default is 'Fhv'
      RhoCp_in - label of the RhoCp series, default is 'RhoCp'
      Tv_in    - label of the virtual temperature series, default is 'Tv_CSAT'
      wA_in    - label of the wA covariance series, default is 'wA'
      rhom_in  - label of the moist air density series, default is 'rhom'
      q_in     - label of the specific humidity series, default is 'q'
      wT_in    - label of the wT covariance series, default is 'wT'
    OUTPUT:
     All outputs are written to the data structure.
      Fh_out   - label of sensible heat flux, default is 'Fh'
    '''
    log.info(' Converting virtual Fh to Fh')
    # get the input series
    Fhv,f,a = qcutils.GetSeriesasMA(ds,Fhv_in)              # get the virtual heat flux
    Tv,f,a = qcutils.GetSeriesasMA(ds,Tv_in)                # get the virtual temperature, C
    TvK = Tv + c.C2K                                        # convert from C to K
    wA,f,a = qcutils.GetSeriesasMA(ds,wA_in)                # get the wA covariance, g/m2/s
    wA = wA * c.g2kg                                        # convert from g/m2/s to kg/m2/s
    q,f,a = qcutils.GetSeriesasMA(ds,q_in)                  # get the specific humidity, kg/kg
    wT,f,wT_a = qcutils.GetSeriesasMA(ds,wT_in)             # get the wT covariance, mK/s
    # get the utility series
    RhoCp,f,a = qcutils.GetSeriesasMA(ds,'RhoCp')           # get rho*Cp
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')             # get the moist air density, kg/m3
    # define local constants
    alpha = 0.51
    # do the conversion
    Fh = Fhv - RhoCp*alpha*TvK*wA/rhom - RhoCp*alpha*q*wT
    # put the calculated sensible heat flux into the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Sensible heat flux from virtual heat flux',
                                           units='W/m2',standard_name='surface_upward_sensible_heat_flux')
    if "height" in wT_a: attr["height"] = wT_a["height"]
    qcutils.CreateSeries(ds,Fh_out,Fh,FList=[Fhv_in,Tv_in,wA_in,q_in,wT_in],Attr=attr)
    if 'FhvtoFh' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', FhvtoFh'
    if qcutils.cfoptionskeylogical(cf,Key='RelaxFhvtoFh'):
        ReplaceWhereMissing(ds.series['Fh'],ds.series['Fh'],ds.series['Fhv'],FlagValue=20)
        if 'RelaxFhvtoFh' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', RelaxFhvtoFh'

def FilterUstar(cf,ds,ustar_in='ustar',ustar_out='ustar_filtered'):
    """
    Filter ustar for low turbulence periods.  The filtering is done by checking the
    friction velocity for each time period.  If ustar is less than or equal to the
    threshold specified in the control file then ustar is set to missing.  If
    the ustar is greater than the threshold, no action is taken.  Filtering is not
    done "in place", a new series is created with the label given in the control file.
    The QC flag is set to 18 to indicate the missing low ustar values.
    
    Usage: qcts.FilterUstar(cf,ds)
    cf: control file object
    ds: data structure object
    """
    if ustar_out not in cf['Variables'].keys(): return
    if 'ustar_threshold' in cf['Variables'][ustar_out].keys():
        log.info(' Filtering ustar to remove values below threshold')
        ustar_threshold = float(cf['Variables'][ustar_out]['ustar_threshold'])
        ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,ustar_in)
        index = numpy.ma.where(ustar<=ustar_threshold)[0]
        ustar = numpy.ma.masked_where(ustar<=ustar_threshold,ustar)
        ustar_flag[index] = 18
        descr = 'ustar filtered for low turbulence conditions (<'+str(ustar_threshold)+')'
        units = qcutils.GetUnitsFromds(ds, ustar_in)
        attr = qcutils.MakeAttributeDictionary(long_name=descr,units=units)
        qcutils.CreateSeries(ds,ustar_out,ustar,Flag=ustar_flag,Attr=attr)
    else:
        log.error(' ustar threshold (ustar_threshold) not found in '+ustar_out+' section of control file')

def get_averages(Data):
    """
        Get daily averages on days when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_averages(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Av = c.missing_value
    elif Num == 48:
        Av = numpy.ma.mean(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Av = numpy.ma.mean(Data[li])
        else:
            Av = c.missing_value
    return Num, Av

def get_laggedcorrelation(x_in,y_in,maxlags):
    """
    Calculate the lagged cross-correlation between 2 1D arrays.
    Taken from the matplotlib.pyplot.xcorr source code.
    PRI added handling of masked arrays.
    """
    lags = numpy.arange(-maxlags,maxlags+1)
    mask = numpy.ma.mask_or(x_in.mask,y_in.mask,copy=True,shrink=False)
    x = numpy.ma.array(x_in,mask=mask,copy=True)
    y = numpy.ma.array(y_in,mask=mask,copy=True)
    x = numpy.ma.compressed(x)
    y = numpy.ma.compressed(y)
    corr = numpy.correlate(x, y, mode=2)
    corr/= numpy.sqrt(numpy.dot(x,x) * numpy.dot(y,y))
    if maxlags is None: maxlags = len(x) - 1
    if maxlags >= len(x) or maxlags < 1:
        raise ValueError('qcts.get_laggedcorrelation: maxlags must be None or strictly positive < %d'%len(x))
    corr = corr[len(x)-1-maxlags:len(x)+maxlags]
    return lags,corr

def get_minmax(Data):
    """
        Get daily minima and maxima on days when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are sample size (Num), minimum (Min) and maximum (Max)
        
        Usage qcts.get_minmax(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Min = c.missing_value
        Max = c.missing_value
    elif Num == 48:
        Min = numpy.ma.min(Data[li])
        Max = numpy.ma.max(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Min = numpy.ma.min(Data[li])
            Max = numpy.ma.max(Data[li])
        else:
            Min = c.missing_value
            Max = c.missing_value
    return Num, Min, Max

def get_nightsums(Data):
    """
        Get nightly sums and averages on nights when no 30-min observations are missing.
        Nights with missing observations return a value of c.missing_value
        Values returned are sample size (Num), sums (Sum) and average (Av)
        
        Usage qcts.get_nightsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(Data.mask == False)[0]
    Num = numpy.size(li)
    if Num == 0:
        Sum = c.missing_value
        Av = c.missing_value
    else:
        x = 0
        for i in range(len(Data)):
            if Data.mask[i] == True:
                x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
            Av = numpy.ma.mean(Data[li])
        else:
            Sum = c.missing_value
            Av = c.missing_value
    
    return Num, Sum, Av

def get_soilaverages(Data):
    """
        Get daily averages of soil water content on days when 15 or fewer 30-min observations are missing.
        Days with 16 or more missing observations return a value of c.missing_value
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_soilaverages(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num > 33:
        Av = numpy.ma.mean(Data[li])
    else:
        Av = c.missing_value
    return Num, Av

def get_subsums(Data):
    """
        Get separate daily sums of positive and negative fluxes when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are positive and negative sample sizes (PosNum and NegNum) and sums (SumPos and SumNeg)
        
        Usage qcts.get_subsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 48:
        pi = numpy.ma.where(Data[li]>0)
        ni = numpy.ma.where(Data[li]<0)
        PosNum = numpy.size(pi)
        NegNum = numpy.size(ni)
        if PosNum > 0:
            SumPos = numpy.ma.sum(Data[pi])
        else:
            SumPos = 0
        if NegNum > 0:
            SumNeg = numpy.ma.sum(Data[ni])
        else:
            SumNeg = 0
    else:
        pi = numpy.ma.where(Data[li]>0)
        ni = numpy.ma.where(Data[li]<0)
        PosNum = numpy.size(pi)
        NegNum = numpy.size(ni)
        SumPos = c.missing_value
        SumNeg = c.missing_value
    return PosNum, NegNum, SumPos, SumNeg

def get_sums(Data):
    """
        Get daily sums when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are sample size (Num) and sum (Sum)
        
        Usage qcts.get_sums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-float(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Sum = c.missing_value
    elif Num == 48:
        Sum = numpy.ma.sum(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
        else:
            Sum = c.missing_value
    return Num, Sum

def get_qcflag(ds):
    """
        Set up flags during ingest of L1 data.
        Identifies missing observations as c.missing_value and sets flag value 1
        
        Usage qcts.get_qcflag(ds)
        ds: data structure
        """
    log.info(' Setting up the QC flags')
    nRecs = len(ds.series['xlDateTime']['Data'])
    for ThisOne in ds.series.keys():
        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        index = numpy.where(ds.series[ThisOne]['Data']==c.missing_value)[0]
        ds.series[ThisOne]['Flag'][index] = numpy.int32(1)

def get_synthetic_fsd(ds):
    """
    Purpose:
     Calculates a time series of synthetic downwelling shortwave radiation.  The
     solar altitude is also output.
    Useage:
     qcts.get_synthetic_fsd(ds)
    Author: PRI
    Date: Sometime in 2014
    """
    log.info(' Calculating synthetic Fsd')
    # get the latitude and longitude
    lat = float(ds.globalattributes["latitude"])
    lon = float(ds.globalattributes["longitude"])
    # get the UTC time from the local time
    ldt_UTC = qcutils.get_UTCfromlocaltime(ds)
    # get the solar altitude
    alt_solar = [pysolar.GetAltitude(lat,lon,dt) for dt in ldt_UTC]
    # get the synthetic downwelling shortwave radiation
    Fsd_syn = [pysolar.GetRadiationDirect(dt,alt) for dt,alt in zip(ldt_UTC,alt_solar)]
    Fsd_syn = numpy.ma.array(Fsd_syn)
    # get the QC flag
    nRecs = len(Fsd_syn)
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    # add the synthetic downwelling shortwave radiation to the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Synthetic downwelling shortwave radiation',units='W/m2',
                                           standard_name='surface_downwelling_shortwave_flux_in_air')
    qcutils.CreateSeries(ds,"Fsd_syn",Fsd_syn,Flag=flag,Attr=attr)
    # add the solar altitude to the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Solar altitude',units='deg',
                                           standard_name='not defined')
    qcutils.CreateSeries(ds,"solar_altitude",alt_solar,Flag=flag,Attr=attr)

def InvertSign(ds,ThisOne):
    log.info(' Inverting sign of '+ThisOne)
    index = numpy.where(abs(ds.series[ThisOne]['Data']-float(c.missing_value))>c.eps)[0]
    ds.series[ThisOne]['Data'][index] = float(-1)*ds.series[ThisOne]['Data'][index]

def InterpolateOverMissing(ds,series='',maxlen=1000):
    """
    Purpose:
     Interpolate over periods of missing data.  Uses linear interpolation.
    Usage:
     qcts.InterpolateOverMissing(ds,series=ThisOne,maxlen=3)
     where ds is the data structure
           ThisOne is a series label
           maxlen is the maximum gap length (hours) to be filled by interpolation
    Side effects:
     Fills gaps.
    Author: PRI
    Date: September 2014
    """
    # check that series is in the data structure
    if series not in ds.series.keys():
        log.error("InterpolateOverMissing: series "+series+" not found in data structure")
        return
    # convert the Python datetime to a number
    DateNum = date2num(ds.series['DateTime']['Data'])
    # get the data
    data_org,flag_org,attr_org = qcutils.GetSeries(ds,series)
    # number of records
    nRecs = len(data_org)
    # index of good values
    iog = numpy.where(abs(data_org-float(c.missing_value))>c.eps)[0]
    # index of missing values
    iom = numpy.where(abs(data_org-float(c.missing_value))<=c.eps)[0]
    # return if there is not enough data to use
    if len(iog)<2:
        log.info(' InterpolateOverMissing: Less than 2 good points available for series '+str(series))
        return
    # linear interpolation function
    f = interpolate.interp1d(DateNum[iog],data_org[iog],bounds_error=False,fill_value=float(c.missing_value))
    # interpolate over the whole time series
    data_int = f(DateNum).astype(numpy.float64)
    # copy the original flag
    flag_int = numpy.copy(flag_org)
    # index of interpolates that are not equal to the missing value
    index = numpy.where(abs(data_int-float(c.missing_value))>c.eps)[0]
    # set the flag for these points
    if len(index)!=0:
        flag_int[index] = numpy.int32(50)
    # restore the original good data
    data_int[iog] = data_org[iog]
    flag_int[iog] = flag_org[iog]
    # now replace data in contiguous blocks of length > min with missing data
    # first, a conditional index, 0 where data is good, 1 where it is missing
    cond_ind = numpy.zeros(nRecs,dtype=numpy.int32)
    cond_ind[iom] = 1
    cond_bool = (cond_ind==1)
    # start and stop indices of contiguous blocks
    for start, stop in qcutils.contiguous_regions(cond_bool):
        # code to handle minimum segment length goes here
        duration = stop - start
        if duration>maxlen:
            #data_int[start:stop+1] = numpy.float(c.missing_value)
            #flag_int[start:stop+1] = flag_org[start:stop+1]
            data_int[start:stop] = numpy.float64(c.missing_value)
            flag_int[start:stop] = flag_org[start:stop]
    # put data_int back into the data structure
    attr_int = dict(attr_org)
    qcutils.CreateSeries(ds,series,data_int,Flag=flag_int,Attr=attr_int)
    if 'InterpolateOverMissing2' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', InterpolateOverMissing2'

def MassmanStandard(cf,ds,Ta_in='Ta',Ah_in='Ah',ps_in='ps',ustar_in='ustar',ustar_out='ustar',L_in='L',L_out ='L',uw_out='uw',vw_out='vw',wT_out='wT',wA_out='wA',wC_out='wC'):
    """
       Massman corrections.
       The steps involved are as follows:
        1) calculate ustar and L using rotated but otherwise uncorrected covariances
       """
    if not qcutils.cfoptionskeylogical(cf,Key='MassmanCorrection'): return
    if 'Massman' not in cf:
        log.info(' Massman section not found in control file, no corrections applied')
        return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MassmanVars'):
        MArgs = ast.literal_eval(cf['FunctionArgs']['MassmanVars'])
        Ta_in = MArgs[0]
        Ah_in = MArgs[1]
        ps_in = MArgs[2]
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MassmanOuts'):
        MOut = ast.literal_eval(cf['FunctionArgs']['MassmanOuts'])
        ustar_in = MOut[0]
        ustar_out = MOut[1]
        L_in = MOut[2]
        L_out = MOut[3]
        uw_out = MOut[4]
        vw_out = MOut[5]
        wT_out = MOut[6]
        wA_out = MOut[7]
        wC_out = MOut[8]
    log.info(' Correcting for flux loss from spectral attenuation')
    zmd = float(cf['Massman']['zmd'])             # z-d for site
    angle = float(cf['Massman']['angle'])         # CSAT3-IRGA separation angle
    CSATarm = float(cf['Massman']['CSATarm'])     # CSAT3 mounting distance
    IRGAarm = float(cf['Massman']['IRGAarm'])     # IRGA mounting distance
    lLat = numpy.ma.sin(numpy.deg2rad(angle)) * IRGAarm
    lLong = CSATarm - (numpy.ma.cos(numpy.deg2rad(angle)) * IRGAarm)
    # *** Massman_1stpass starts here ***
    #  The code for the first and second passes is very similar.  It would be useful to make them the
    #  same and put into a loop to reduce the nu,ber of lines in this function.
    # calculate ustar and Monin-Obukhov length from rotated but otherwise uncorrected covariances
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    nRecs = numpy.size(Ta)
    u,f,a = qcutils.GetSeriesasMA(ds,'u')
    uw,f,a = qcutils.GetSeriesasMA(ds,'uw')
    vw,f,a = qcutils.GetSeriesasMA(ds,'vw')
    wT,f,a = qcutils.GetSeriesasMA(ds,'wT')
    wC,f,a = qcutils.GetSeriesasMA(ds,'wC')
    wA,f,a = qcutils.GetSeriesasMA(ds,'wA')
    if ustar_in not in ds.series.keys():
        ustarm = numpy.ma.sqrt(numpy.ma.sqrt(uw ** 2 + vw ** 2))
    else:
        ustarm,f,a = qcutils.GetSeriesasMA(ds,ustar_in)
    if L_in not in ds.series.keys():
        Lm = mf.molen(Ta, Ah, ps, ustarm, wT, fluxtype='kinematic')
    else:
        Lm,f,a = qcutils.GetSeriesasMA(ds,Lm_in)
    # now calculate z on L
    zoLm = zmd / Lm
    # start calculating the correction coefficients for approximate corrections
    #  create nxMom, nxScalar and alpha series with their unstable values by default
    nxMom, nxScalar, alpha = qcutils.nxMom_nxScalar_alpha(zoLm)
    # now calculate the fxMom and fxScalar coefficients
    fxMom = nxMom * u / zmd
    fxScalar = nxScalar * u / zmd
    # compute spectral filters
    tao_eMom = ((c.lwVert / (5.7 * u)) ** 2) + ((c.lwHor / (2.8 * u)) ** 2)
    tao_ewT = ((c.lwVert / (8.4 * u)) ** 2) + ((c.lTv / (4.0 * u)) ** 2)
    tao_ewIRGA = ((c.lwVert / (8.4 * u)) ** 2) + ((c.lIRGA / (4.0 * u)) ** 2) \
                 + ((lLat / (1.1 * u)) ** 2) + ((lLong / (1.05 * u)) ** 2)
    tao_b = c.Tb / 2.8
    # calculate coefficients
    bMom = qcutils.bp(fxMom,tao_b)
    bScalar = qcutils.bp(fxScalar,tao_b)
    pMom = qcutils.bp(fxMom,tao_eMom)
    pwT = qcutils.bp(fxScalar,tao_ewT)
    # calculate corrections for momentum and scalars
    rMom = qcutils.r(bMom, pMom, alpha)        # I suspect that rMom and rwT are the same functions
    rwT = qcutils.r(bScalar, pwT, alpha)
    # determine approximately-true Massman fluxes
    uwm = uw / rMom
    vwm = vw / rMom
    wTm = wT / rwT
    # *** Massman_1stpass ends here ***
    # *** Massman_2ndpass starts here ***
    # we have calculated the first pass corrected momentum and temperature covariances, now we use
    # these to calculate the final corrections
    #  first, get the 2nd pass corrected friction velocity and Monin-Obukhov length
    ustarm = numpy.ma.sqrt(numpy.ma.sqrt(uwm ** 2 + vwm ** 2))
    Lm = mf.molen(Ta, Ah, ps, ustarm, wTm, fluxtype='kinematic')
    zoLm = zmd / Lm
    nxMom, nxScalar, alpha = qcutils.nxMom_nxScalar_alpha(zoLm)
    fxMom = nxMom * (u / zmd)
    fxScalar = nxScalar * (u / zmd)
    # calculate coefficients
    bMom = qcutils.bp(fxMom,tao_b)
    bScalar = qcutils.bp(fxScalar,tao_b)
    pMom = qcutils.bp(fxMom,tao_eMom)
    pwT = qcutils.bp(fxScalar,tao_ewT)
    pwIRGA = qcutils.bp(fxScalar,tao_ewIRGA)
    # calculate corrections for momentum and scalars
    rMom = qcutils.r(bMom, pMom, alpha)
    rwT = qcutils.r(bScalar, pwT, alpha)
    rwIRGA = qcutils.r(bScalar, pwIRGA, alpha)
    # determine true fluxes
    uwM = uw / rMom
    vwM = vw / rMom
    wTM = wT / rwT
    wCM = wC / rwIRGA
    wAM = wA / rwIRGA
    ustarM = numpy.ma.sqrt(numpy.ma.sqrt(uwM ** 2 + vwM ** 2))
    LM = mf.molen(Ta, Ah, ps, ustarM, wTM, fluxtype='kinematic')
    # write the 2nd pass Massman corrected covariances to the data structure
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true ustar',units='m/s')
    qcutils.CreateSeries(ds,ustar_out,ustarM,FList=['uw','vw'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Obukhov Length',units='m')
    qcutils.CreateSeries(ds,L_out,LM,FList=[Ta_in,Ah_in,ps_in,'wT'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(uw)',units='m2/s2')
    qcutils.CreateSeries(ds,uw_out,uwM,FList=['uw',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(vw)',units='m2/s2')
    qcutils.CreateSeries(ds,vw_out,vwM,FList=['vw',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(wT)',units='mC/s')
    qcutils.CreateSeries(ds,wT_out,wTM,FList=['wT',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(wA)',units='g/m2/s')
    qcutils.CreateSeries(ds,wA_out,wAM,FList=['wA',L_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Massman true Cov(wC)',units='mg/m2/s')
    qcutils.CreateSeries(ds,wC_out,wCM,FList=['wC',L_out],Attr=attr)
    # *** Massman_2ndpass ends here ***
    
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='MassmanFlag') and cf['General']['MassmanFlag'] == 'True':
        keys = [ustar_out,L_out,uw_out,vw_out,wT_out,wA_out,wC_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(int)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(12)

def MergeSeriesUsingDict(ds,merge_order=""):
    """ Merge series as defined in the ds.merge dictionary."""
    # check that ds has a "merge" attribute
    if "merge" not in dir(ds): raise Exception("MergeSeriesUsingDict: No merge dictionary in ds")
    if merge_order not in ds.merge.keys():
        msg = "MergeSeriesUsingDict: merge_order ("+merge_order+") not found in merge dictionary"
        log.info(msg)
        return
    # loop over the entries in ds.merge
    for target in ds.merge[merge_order].keys():
        srclist = ds.merge[merge_order][target]["source"]
        log.info(' Merging '+str(srclist)+' ==> '+target)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series '+srclist[0]+' not found')
            continue
        data = ds.series[srclist[0]]['Data'].copy()
        flag1 = ds.series[srclist[0]]['Flag'].copy()
        flag2 = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        tmplist = list(srclist)
        tmplist.remove(tmplist[0])
        for label in tmplist:
            if label in ds.series.keys():
                SeriesNameString = SeriesNameString+', '+label
                index = numpy.where(numpy.mod(flag1,10)==0)[0]         # find the elements with flag = 0, 10, 20 etc
                flag2[index] = 0                                        # set them all to 0
                if label=="Fg":
                    index = numpy.where(flag2==22)[0]
                    if len(index)!=0: flag2[index] = 0
                index = numpy.where(flag2!=0)[0]                        # index of flag values other than 0,10,20,30 ...
                data[index] = ds.series[label]['Data'][index].copy()  # replace bad primary with good secondary
                flag1[index] = ds.series[label]['Flag'][index].copy()
            else:
                log.error(" MergeSeries: secondary input series "+label+" not found")
        attr["long_name"] = attr["long_name"]+", merged from " + SeriesNameString
        qcutils.CreateSeries(ds,target,data,Flag=flag1,Attr=attr)
    del ds.merge[merge_order]

def MergeHumidities(cf,ds):
    if "Ah" not in cf["Variables"] and "RH" not in cf["Variables"] and "q" not in cf["Variables"]:
        log.error(" MergeHumidities: No humidities found in control file, returning ...")
        return
    if "Ah" in cf["Variables"]:
        MergeSeries(cf,ds,'Ah',[0,10])
    if "RH" in cf["Variables"]:
        MergeSeries(cf,ds,'RH',[0,10])
    if "q" in cf["Variables"]:
        MergeSeries(cf,ds,'q',[0,10])

def MergeSeries(cf,ds,series,okflags):
    """
        Merge two series of data to produce one series containing the best data from both.
        Calling syntax is: MergeSeries(cf,ds,series,okflags)
         where ds is the data structure containing all series
               series (str) is the label of the destination series
               okflags (list) is a list of QC flag values for which the data is considered acceptable
        If the QC flag for Primary is in okflags, the value from Primary is placed in destination.
        If the QC flag for Primary is not in okflags but the QC flag for Secondary is, the value
        from Secondary is placed in Destination.
        """
    # check to see if the series is specified in the control file
    section = qcutils.get_cfsection(cf,series=series)
    if len(section)==0: return
    # check to see if the entry for series in the control file has the MergeSeries key
    if 'MergeSeries' not in cf[section][series].keys(): return
    # check to see if the series has already been merged
    if series in ds.mergeserieslist: return
    # now get the source list and the standard name
    srclist, standardname = qcutils.GetMergeSeriesKeys(cf,series,section=section)
    nSeries = len(srclist)
    if nSeries==0:
        log.info(' MergeSeries: no input series specified for '+str(series))
        return
    if nSeries==1:
        log.info(' Merging '+str(srclist)+'==>'+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found for'+str(series))
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
    else:
        log.info(' Merging '+str(srclist)+'==>'+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found')
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            if ThisOne in ds.series.keys():
                SeriesNameString = SeriesNameString+', '+ThisOne
                indx1 = numpy.zeros(numpy.size(data),dtype=numpy.int)
                indx2 = numpy.zeros(numpy.size(data),dtype=numpy.int)
                for okflag in okflags:
                    index = numpy.where((flag==okflag))[0]                             # index of acceptable primary values
                    indx1[index] = 1                                                   # set primary index to 1 when primary good
                    index = numpy.where((ds.series[ThisOne]['Flag']==okflag))[0]       # same process for secondary
                    indx2[index] = 1
                index = numpy.where((indx1!=1)&(indx2==1))[0]           # index where primary bad but secondary good
                data[index] = ds.series[ThisOne]['Data'][index]         # replace bad primary with good secondary
                flag[index] = ds.series[ThisOne]['Flag'][index]
            else:
                log.error('  MergeSeries: secondary input series'+ThisOne+'not found')
    ds.mergeserieslist.append(series)
    #attr = qcutils.MakeAttributeDictionary(long_name='Merged from '+SeriesNameString,
                             #standard_name=standardname,units=SeriesUnitString)
    attr["long_name"] = attr["long_name"]+", merged from " + SeriesNameString
    qcutils.CreateSeries(ds,series,data,Flag=flag,Attr=attr)

def PT100(ds,T_out,R_in,m):
    log.info(' Calculating temperature from PT100 resistance')
    R,f,a = qcutils.GetSeriesasMA(ds,R_in)
    R = m*R
    T = (-c.PT100_alpha+numpy.sqrt(c.PT100_alpha**2-4*c.PT100_beta*(-R/100+1)))/(2*c.PT100_beta)
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated PT100 temperature using '+str(R_in),units='degC')
    qcutils.CreateSeries(ds,T_out,T,FList=[R_in],Attr=attr)

def ReplaceRotatedCovariance(cf,ds,rot_cov_label,non_cov_label):
    log.info(' Replacing missing '+rot_cov_label+' when '+non_cov_label+' is good')
    cr_data,cr_flag,cr_attr = qcutils.GetSeriesasMA(ds,rot_cov_label)
    cn_data,cn_flag,cn_attr = qcutils.GetSeriesasMA(ds,non_cov_label)
    index = numpy.where((numpy.ma.getmaskarray(cr_data)==True)&
                           (numpy.ma.getmaskarray(cn_data)==False))[0]
    #index = numpy.ma.where((numpy.ma.getmaskarray(cr_data)==True)&
                           #(numpy.ma.getmaskarray(cn_data)==False))[0]
    if len(index)!=0:
        ds.series[rot_cov_label]['Data'][index] = cn_data[index]
        ds.series[rot_cov_label]['Flag'][index] = numpy.int32(20)
    return

def ReplaceOnDiff(cf,ds,series=''):
    # Gap fill using data from alternate sites specified in the control file
    ts = ds.globalattributes['time_step']
    if len(series)!=0:
        ds_alt = {}                     # create a dictionary for the data from alternate sites
        open_ncfiles = []               # create an empty list of open netCDF files
        for ThisOne in series:          # loop over variables in the series list
            # has ReplaceOnDiff been specified for this series?
            if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'ReplaceOnDiff'):
                # loop over all entries in the ReplaceOnDiff section
                for Alt in cf['Variables'][ThisOne]['ReplaceOnDiff'].keys():
                    if 'FileName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                        alt_filename = cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['FileName']
                        if 'AltVarName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            alt_varname = cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['AltVarName']
                        else:
                            alt_varname = ThisOne
                        if alt_filename not in open_ncfiles:
                            n = len(open_ncfiles)
                            open_ncfiles.append(alt_filename)
                            ds_alt[n] = qcio.nc_read_series_file(alt_filename)
                        else:
                            n = open_ncfiles.index(alt_filename)
                        if 'Transform' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            AltDateTime = ds_alt[n].series['DateTime']['Data']
                            AltSeriesData = ds_alt[n].series[alt_varname]['Data']
                            TList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Transform'])
                            for TListEntry in TList:
                                qcts.TransformAlternate(TListEntry,AltDateTime,AltSeriesData,ts=ts)
                        if 'Range' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            RList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Range'])
                            for RListEntry in RList:
                                qcts.ReplaceWhenDiffExceedsRange(ds.series['DateTime']['Data'],ds.series[ThisOne],
                                                                 ds.series[ThisOne],ds_alt[n].series[alt_varname],
                                                                 RListEntry)
                    elif 'AltVarName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                        alt_varname = ThisOne
                        if 'Range' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            RList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Range'])
                            for RListEntry in RList:
                                qcts.ReplaceWhenDiffExceedsRange(ds.series['DateTime']['Data'],ds.series[ThisOne],
                                                                 ds.series[ThisOne],ds.series[alt_varname],
                                                                 RListEntry)
                    else:
                        log.error('ReplaceOnDiff: Neither AltFileName nor AltVarName given in control file')
    else:
        log.error('ReplaceOnDiff: No input series specified')

def ReplaceWhereMissing(Destination,Primary,Secondary,FlagOffset=None,FlagValue=None):
    #print time.strftime('%X')+' Merging series '+Primary+' and '+Secondary+' into '+Destination
    p_data = Primary['Data'].copy()
    p_flag = Primary['Flag'].copy()
    s_data = Secondary['Data'].copy()
    s_flag = Secondary['Flag'].copy()
    if numpy.size(p_data)>numpy.size(s_data):
        p_data = p_data[0:numpy.size(s_data)]
    if numpy.size(s_data)>numpy.size(p_data):
        s_data = s_data[0:numpy.size(p_data)]
    index = numpy.where((abs(p_data-float(c.missing_value))<c.eps)&
                        (abs(s_data-float(c.missing_value))>c.eps))[0]
    p_data[index] = s_data[index]
    if FlagValue is None and FlagOffset is not None:
        p_flag[index] = s_flag[index] + numpy.int32(FlagOffset)
    elif FlagValue is not None and FlagOffset is None:
        p_flag[index] = numpy.int32(FlagValue)
    else:
        p_flag[index] = s_flag[index]
    Destination['Data'] = Primary['Data'].copy()
    Destination['Flag'] = Primary['Flag'].copy()
    Destination['Data'][0:len(p_data)] = p_data
    Destination['Flag'][0:len(p_flag)] = p_flag
    Destination['Attr']['long_name'] = 'Merged from original and alternate'
    Destination['Attr']['units'] = Primary['Attr']['units']

def ReplaceWhenDiffExceedsRange(DateTime,Destination,Primary,Secondary,RList):
    #print time.strftime('%X')+' Replacing '+Primary+' with '+Secondary+' when difference exceeds threshold'
    # get the primary data series
    p_data = numpy.ma.array(Primary['Data'])
    p_flag = Primary['Flag'].copy()
    # get the secondary data series
    s_data = numpy.ma.array(Secondary['Data'])
    s_flag = Secondary['Flag'].copy()
    # truncate the longest series if the sizes do not match
    if numpy.size(p_data)!=numpy.size(s_data):
        log.warning(' ReplaceWhenDiffExceedsRange: Series lengths differ, longest will be truncated')
        if numpy.size(p_data)>numpy.size(s_data):
            p_data = p_data[0:numpy.size(s_data)]
        if numpy.size(s_data)>numpy.size(p_data):
            s_data = s_data[0:numpy.size(p_data)]
    # get the difference between the two data series
    d_data = p_data-s_data
    # normalise the difference if requested
    if RList[3]=='s':
        d_data = (p_data-s_data)/s_data
    elif RList[3]=='p':
        d_data = (p_data-s_data)/p_data
    #si = qcutils.GetDateIndex(DateTime,RList[0],0)
    #ei = qcutils.GetDateIndex(DateTime,RList[1],0)
    Range = RList[2]
    Upper = float(Range[0])
    Lower = float(Range[1])
    index = numpy.ma.where((abs(d_data)<Lower)|(abs(d_data)>Upper))
    p_data[index] = s_data[index]
    p_flag[index] = 35
    Destination['Data'] = numpy.ma.filled(p_data,float(c.missing_value))
    Destination['Flag'] = p_flag.copy()
    Destination['Attr']['long_name'] = 'Replaced original with alternate when difference exceeded threshold'
    Destination['Attr']['units'] = Primary['Attr']['units']

def savitzky_golay(y, window_size, order, deriv=0):
    ''' Apply Savitsky-Golay low-pass filter to data.'''
    try:
        window_size = numpy.abs(numpy.int(window_size))
        order = numpy.abs(numpy.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = numpy.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))
    return numpy.convolve( m, y, mode='valid')

def Square(Series):
    tmp = numpy.array([c.missing_value]*numpy.size(Series),Series.dtype)
    index = numpy.where(Series!=float(c.missing_value))[0]
    tmp[index] = Series[index] ** 2
    return tmp

def SquareRoot(Series):
    tmp = numpy.array([c.missing_value]*numpy.size(Series),Series.dtype)
    index = numpy.where(Series!=float(c.missing_value))[0]
    tmp[index] = Series[index] ** .5
    return tmp

def TaFromTv(cf,ds,Ta_out='Ta_CSAT',Tv_in='Tv_CSAT',Ah_in='Ah',RH_in='RH',q_in='q',ps_in='ps'):
    # Calculate the air temperature from the virtual temperature, the
    # absolute humidity and the pressure.
    # NOTE: the virtual temperature is used in place of the air temperature
    #       to calculate the vapour pressure from the absolute humidity, the
    #       approximation involved here is of the order of 1%.
    log.info(' Calculating Ta from Tv')
    # check to see if we have enough data to proceed
    if Tv_in not in ds.series.keys():
        log.error(" TaFromTv: sonic virtual temperature ("+str(Tv_in)+") not found in data structure")
        return
    if Ah_in not in ds.series.keys() and RH_in not in ds.series.keys() and q_in not in ds.series.keys():
        labstr = str(Ah_in)+","+str(RH_in)+","+str(q_in)
        log.error(" TaFromTv: no humidity data ("+labstr+") found in data structure")
        return
    if ps_in not in ds.series.keys():
        log.error(" TaFromTv: pressure ("+str(ps_in)+") not found in data structure")
        return
    # we seem to have enough to continue
    Tv,f,a = qcutils.GetSeriesasMA(ds,Tv_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    if Ah_in in ds.series.keys():
        Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
        vp = mf.vapourpressure(Ah,Tv)
        mr = mf.mixingratio(ps,vp)
        q = mf.specifichumidity(mr)
    elif RH_in in ds.series.keys():
        RH,f,a = qcutils.GetSeriesasMA(ds,RH_in)
        q = mf.specifichumidityfromRH(RH,Tv,ps)
    elif q_in in ds.series.keys():
        q,f,a = qcutils.GetSeriesasMA(ds,q_in)
    Ta_data = mf.tafromtv(Tv,q)
    nRecs = int(ds.globalattributes['nc_nrecs'])
    Ta_flag = numpy.zeros(nRecs,numpy.int32)
    mask = numpy.ma.getmask(Ta_data)
    index = numpy.where(mask.astype(numpy.int32)==1)
    Ta_flag[index] = 15
    attr = qcutils.MakeAttributeDictionary(long_name='Ta calculated from Tv using '+Tv_in,units='C',standard_name='air_temperature')
    qcutils.CreateSeries(ds,Ta_out,Ta_data,Flag=Ta_flag,Attr=attr)
    if 'TaFromTv' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+', TaFromTv'

def TransformAlternate(TList,DateTime,Series,ts=30):
    # Apply polynomial transform to data series being used as replacement data for gap filling
    #print time.strftime('%X')+' Applying polynomial transform to '+ThisOne
    si = qcutils.GetDateIndex(DateTime,TList[0],ts=ts,default=0,match='exact')
    ei = qcutils.GetDateIndex(DateTime,TList[1],ts=ts,default=-1,match='exact')
    Series = numpy.ma.masked_where(abs(Series-float(c.missing_value))<c.eps,Series)
    Series[si:ei] = qcutils.polyval(TList[2],Series[si:ei])
    Series = numpy.ma.filled(Series,float(c.missing_value))
