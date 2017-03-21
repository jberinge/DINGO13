import ast
import copy
import constants as c
import datetime
import numpy
import time
import qcts
import qcutils
import logging

log = logging.getLogger('qc.ck')

def cliptorange(data, lower, upper):
    data = rangecheckserieslower(data,lower)
    data = rangecheckseriesupper(data,upper)
    return data

def rangecheckserieslower(data,lower):
    if lower is None:
        log.info(' rangecheckserieslower: no lower bound set')
        return data
    if numpy.ma.isMA(data):
        data = numpy.ma.masked_where(data<lower,data)
    else:
        index = numpy.where((abs(data-numpy.float64(c.missing_value))>c.eps)&(data<lower))[0]
        data[index] = numpy.float64(c.missing_value)
    return data

def rangecheckseriesupper(data,upper):
    if upper is None:
        log.info(' rangecheckserieslower: no upper bound set')
        return data
    if numpy.ma.isMA(data):
        data = numpy.ma.masked_where(data>upper,data)
    else:
        index = numpy.where((abs(data-numpy.float64(c.missing_value))>c.eps)&(data>upper))[0]
        data[index] = numpy.float64(c.missing_value)
    return data

def CoordinateFluxGaps(cf,ds,Fc_in='Fc',Fe_in='Fe',Fh_in='Fh'):
    if not qcutils.cfoptionskeylogical(cf,Key='CoordinateFluxGaps'): return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='gapsvars'):
        vars = ast.literal_eval(cf['FunctionArgs']['gapsvars'])
        Fc_in = vars[0]
        Fe_in = vars[1]
        Fh_in = vars[2]
    Fc,f,a = qcutils.GetSeriesasMA(ds,Fc_in)
    Fe,f,a = qcutils.GetSeriesasMA(ds,Fe_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    # April 2015 PRI - changed numpy.ma.where to numpy.where
    index = numpy.where((numpy.ma.getmaskarray(Fc)==True)|
                        (numpy.ma.getmaskarray(Fe)==True)|
                        (numpy.ma.getmaskarray(Fh)==True))[0]
    #index = numpy.ma.where((numpy.ma.getmaskarray(Fc)==True)|
                           #(numpy.ma.getmaskarray(Fe)==True)|
                           #(numpy.ma.getmaskarray(Fh)==True))[0]
    # the following for ... in loop is not necessary
    for i in range(len(index)):
        j = index[i]
        if Fc.mask[j]==False:
            Fc.mask[j]=True
            Fc[j] = numpy.float64(c.missing_value)
            ds.series[Fc_in]['Flag'][j] = numpy.int32(19)
        if Fe.mask[j]==False:
            Fe.mask[j]=True
            Fe[j] = numpy.float64(c.missing_value)
            ds.series[Fe_in]['Flag'][j] = numpy.int32(19)           
        if Fh.mask[j]==False:
            Fh.mask[j]=True
            Fh[j] = numpy.float64(c.missing_value)
            ds.series[Fh_in]['Flag'][j] = numpy.int32(19)
    ds.series[Fc_in]['Data']=numpy.ma.filled(Fc,float(c.missing_value))
    ds.series[Fe_in]['Data']=numpy.ma.filled(Fe,float(c.missing_value))
    ds.series[Fh_in]['Data']=numpy.ma.filled(Fh,float(c.missing_value))
    log.info(' Finished gap co-ordination')

def CreateNewSeries(cf,ds):
    '''Create a new series using the MergeSeries or AverageSeries instructions.'''
    log.info(' Checking for new series to create')
    for ThisOne in cf['Variables'].keys():
        if 'MergeSeries' in cf['Variables'][ThisOne].keys():
            qcts.MergeSeries(cf,ds,ThisOne,[0,10,20,30,40,50])
        if 'AverageSeries' in cf['Variables'][ThisOne].keys():
            qcts.AverageSeriesByElements(cf,ds,ThisOne)

def do_IRGAcheck(cf,ds):
    """
    Purpose:
     Decide which IRGA check routine to use depending on the setting
     of the "irga_type" key in the [Options] section of the control
     file.  The default is Li7500.
    Usage:
    Author: PRI
    Date: September 2015
    """
    irga_list = ["li7500","li7500a","ec155"]
    # get the IRGA type from the control file
    irga_type = qcutils.get_keyvaluefromcf(cf,["Options"],"irga_type", default="li7500")
    # remove any hyphens or spaces
    for item in ["-"," "]:
        if item in irga_type: irga_type = irga_type.replace(item,"")
    # check the IRGA type against the list of suppprted devices
    if irga_type.lower() not in irga_list:
        msg = " Unrecognised IRGA type "+irga_type+" given in control file, IRGA checks skipped ..."
        log.error(msg)
        return
    # do the IRGA checks
    if irga_type.lower()=="li7500" or irga_type.lower()=="li7500a":
        ds.globalattributes["irga_type"] = irga_type
        do_7500check(cf,ds)
    elif irga_type.lower()=="ec155":
        ds.globalattributes["irga_type"] = irga_type
        do_EC155check(cf,ds)
    else:
        msg = " Unsupported IRGA type "+irga_type+", contact the devloper ..."
        log.error(msg)
        return

def do_7500check(cf,ds):
    '''Rejects data values for series specified in LI75List for times when the Diag_7500
       flag is non-zero.  If the Diag_7500 flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       LI75Lisat.  Additional checks are done for AGC_7500 (the LI-7500 AGC value),
       Ah_7500_Sd (standard deviation of absolute humidity) and Cc_7500_Sd (standard
       deviation of CO2 concentration).'''
    if "Diag_7500" not in ds.series.keys():
        msg = " Diag_7500 not found in data, skipping 7500 checks ..."
        log.warning(msg)
        return
    log.info(' Doing the 7500 check')
    LI75List = ['Ah_7500_Av','Cc_7500_Av','Ah_7500_Sd','Cc_7500_Sd',
                'UzA','UxA','UyA','UzC','UxC','UyC']
    index = numpy.where(ds.series['Diag_7500']['Flag']!=0)
    log.info('  7500Check: Diag_7500 ' + str(numpy.size(index)))
    LI75_dependents = []
    for item in ['AGC_7500','Ah_7500_Sd','Cc_7500_Sd','AhAh','CcCc']:
        if item in ds.series.keys(): LI75_dependents.append(item)
    if "Ah_7500_Sd" and "AhAh" in LI75_dependents: LI75_dependents.remove("AhAh")
    if "Cc_7500_Sd" and "CcCc" in LI75_dependents: LI75_dependents.remove("CcCc")
    for item in LI75_dependents:
        if item in ds.series.keys():
            index = numpy.where(ds.series[item]['Flag']!=0)
            log.info('  7500Check: '+item+' rejected '+str(numpy.size(index))+' points')
            ds.series['Diag_7500']['Flag'] = ds.series['Diag_7500']['Flag'] + ds.series[item]['Flag']
    index = numpy.where((ds.series['Diag_7500']['Flag']!=0))
    log.info('  7500Check: Total ' + str(numpy.size(index)))
    for ThisOne in LI75List:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            log.error('  qcck.do_7500check: series '+str(ThisOne)+' in LI75List not found in ds.series')
    if '7500Check' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',7500Check'

def do_EC155check(cf,ds):
    """
    Purpose:
    Usage:
    Author: PRI
    Date: September 2015
    """
    # check to see if we have a Diag_IRGA series to work with
    if "Diag_IRGA" not in ds.series.keys():
        msg = " Diag_IRGA not found in data, skipping IRGA checks ..."
        log.warning(msg)
        return
    # seems OK to continue
    log.info(' Doing the EC155 check')
    # list of series that depend on IRGA data quality
    EC155_list = ['H2O_IRGA_Av','CO2_IRGA_Av','H2O_IRGA_Sd','CO2_IRGA_Sd',
                 'UzH','UxH','UyH','UzC','UxC','UyC']
    index = numpy.where(ds.series['Diag_IRGA']['Flag']!=0)
    log.info('  EC155Check: Diag_IRGA rejects ' + str(numpy.size(index)))
    EC155_dependents = []
    for item in ['Signal_H2O','Signal_CO2','H2O_IRGA_Sd','CO2_IRGA_Sd']:
        if item in ds.series.keys(): EC155_dependents.append(item)
    for item in EC155_dependents:
        index = numpy.where(ds.series[item]['Flag']!=0)
        log.info('  EC155Check: '+item+' rejected '+str(numpy.size(index))+' points')
        ds.series['Diag_IRGA']['Flag'] = ds.series['Diag_IRGA']['Flag'] + ds.series[item]['Flag']
    index = numpy.where((ds.series['Diag_IRGA']['Flag']!=0))
    log.info('  EC155Check: Total rejected ' + str(numpy.size(index)))
    for ThisOne in EC155_list:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            log.error(' do_EC155check: series '+str(ThisOne)+' in EC155 list not found in data structure')
    if 'EC155Check' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',EC155Check'

def CoordinateAh7500AndFcGaps(cf,ds,Fcvar='Fc'):
    '''Cleans up Ah_7500_Av based upon Fc gaps to for QA check on Ah_7500_Av v Ah_HMP.'''
    if not qcutils.cfoptionskeylogical(cf,Key='CoordinateAh7500&FcGaps'): return
    log.info(' Doing the Ah_7500 check')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='AhcheckFc'):
        Fclist = ast.literal_eval(cf['FunctionArgs']['AhcheckFc'])
        Fcvar = Fclist[0]
    
    # index1  Index of bad Ah_7500_Av observations
    index1 = numpy.where((ds.series['Ah_7500_Av']['Flag']!=0) & (ds.series['Ah_7500_Av']['Flag']!=10))
    
    # index2  Index of bad Fc observations
    index2 = numpy.where((ds.series[Fcvar]['Flag']!=0) & (ds.series[Fcvar]['Flag']!=10))
    
    ds.series['Ah_7500_Av']['Data'][index2] = numpy.float64(c.missing_value)
    ds.series['Ah_7500_Av']['Flag'][index2] = ds.series[Fcvar]['Flag'][index2]
    ds.series['Ah_7500_Av']['Flag'][index1] = ds.series['Ah_7500_Av']['Flag'][index1]
    if 'CoordinateAh7500AndFcGaps' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',CoordinateAh7500AndFcGaps'

def do_CSATcheck(cf,ds):
    '''Rejects data values for series specified in CSATList for times when the Diag_CSAT
       flag is non-zero.  If the Diag_CSAT flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       CSATList.'''
    if "Diag_CSAT" not in ds.series.keys():
        msg = " Diag_CSAT not found in data, skipping CSAT checks ..."
        log.warning(msg)
        return
    log.info(' Doing the CSAT check')
    CSAT_all = ['Ux','Uy','Uz',
                'Ws_CSAT','Wd_CSAT','Wd_CSAT_Compass','Tv_CSAT',
                'UzT','UxT','UyT','UzA','UxA','UyA','UzC','UxC','UyC',
                'UxUz','UyUz','UxUy','UxUx','UyUy','UzUz']
    CSAT_list = []
    for item in CSAT_all:
        if item in ds.series.keys(): CSAT_list.append(item)
    index = numpy.where(ds.series['Diag_CSAT']['Flag']!=0)
    log.info('  CSATCheck: Diag_CSAT ' + str(numpy.size(index)))
    for ThisOne in CSAT_list:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(3)
        else:
            log.error('  qcck.do_CSATcheck: series '+str(ThisOne)+' in CSATList not found in ds.series')
    if 'CSATCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',CSATCheck'

def do_dependencycheck(cf,ds,section='',series='',code=23,mode="quiet"):
    if len(section)==0 and len(series)==0: return
    if len(section)==0: section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if "DependencyCheck" not in cf[section][series].keys(): return
    if "Source" not in cf[section][series]["DependencyCheck"]:
        msg = " DependencyCheck: keyword Source not found for series "+series+", skipping ..."
        log.error(msg)
        return
    if mode=="verbose":
        msg = " Doing DependencyCheck for "+series
        log.info(msg)
    # get the precursor source list from the control file
    source_list = ast.literal_eval(cf[section][series]["DependencyCheck"]["Source"])
    # get the data
    dependent_data,dependent_flag,dependent_attr = qcutils.GetSeriesasMA(ds,series)
    # loop over the precursor source list
    for item in source_list:
        # check the precursor is in the data structure
        if item not in ds.series.keys():
            msg = " DependencyCheck: "+series+" precursor series "+item+" not found, skipping ..."
            continue
        # get the precursor data
        precursor_data,precursor_flag,precursor_attr = qcutils.GetSeriesasMA(ds,item)
        # mask the dependent data where the precurso is masked
        dependent_data = numpy.ma.masked_where(numpy.ma.getmaskarray(precursor_data)==True,dependent_data)
        # get an index of masked precursor data
        index = numpy.ma.where(numpy.ma.getmaskarray(precursor_data)==True)[0]
        # set the dependent QC flag
        dependent_flag[index] = numpy.int32(code)
    # put the data back into the data structure
    if series=="Fc":
        pass
    dependent_attr["DependencyCheck_source"] = str(source_list)
    qcutils.CreateSeries(ds,series,dependent_data,Flag=dependent_flag,Attr=dependent_attr)
    if 'do_dependencychecks' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_dependencychecks'

def do_diurnalcheck(cf,ds,section='',series='',code=5):
    if 'DiurnalCheck' not in cf[section][series].keys(): return
    if 'NumSd' not in cf[section][series]['DiurnalCheck'].keys(): return
    dt = float(ds.globalattributes['time_step'])
    n = int((60./dt) + 0.5)             #Number of timesteps per hour
    nInts = int((1440.0/dt)+0.5)        #Number of timesteps per day
    Av = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Sd = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    NSd = numpy.array(eval(cf[section][series]['DiurnalCheck']['NumSd']),dtype=float)
    for m in range(1,13):
        mindex = numpy.where(ds.series['Month']['Data']==m)[0]
        if len(mindex)!=0:
            lHdh = ds.series['Hdh']['Data'][mindex]
            l2ds = ds.series[series]['Data'][mindex]
            for i in range(nInts):
                li = numpy.where(abs(lHdh-(float(i)/float(n))<c.eps)&(l2ds!=float(c.missing_value)))
                if numpy.size(li)!=0:
                    Av[i] = numpy.mean(l2ds[li])
                    Sd[i] = numpy.std(l2ds[li])
                else:
                    Av[i] = float(c.missing_value)
                    Sd[i] = float(c.missing_value)
            Lwr = Av - NSd[m-1]*Sd
            Upr = Av + NSd[m-1]*Sd
            hindex = numpy.array(n*lHdh,int)
            index = numpy.where(((l2ds!=float(c.missing_value))&(l2ds<Lwr[hindex]))|
                                ((l2ds!=float(c.missing_value))&(l2ds>Upr[hindex])))[0] + mindex[0]
            ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['diurnalcheck_numsd'] = cf[section][series]['DiurnalCheck']['NumSd']
    if 'DiurnalCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',DiurnalCheck'

def do_excludedates(cf,ds,section='',series='',code=6):
    if 'ExcludeDates' not in cf[section][series].keys(): return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeDates'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeDateList = ast.literal_eval(cf[section][series]['ExcludeDates'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeDateList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeDateList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        ds.series[series]['Data'][si:ei] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][si:ei] = numpy.int32(code)
        ds.series[series]['Attr']['ExcludeDates_'+str(i)] = cf[section][series]['ExcludeDates'][str(i)]
    if 'ExcludeDates' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',ExcludeDates'

def do_excludehours(cf,ds,section='',series='',code=7):
    if 'ExcludeHours' not in cf[section][series].keys(): return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeHours'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeHourList = ast.literal_eval(cf[section][series]['ExcludeHours'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeHourList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeHourList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        for j in range(len(ExcludeHourList[2])):
            ExHr = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').hour
            ExMn = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').minute
            index = numpy.where((ds.series['Hour']['Data'][si:ei]==ExHr)&
                                (ds.series['Minute']['Data'][si:ei]==ExMn))[0] + si
            ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['ExcludeHours_'+str(i)] = cf[section][series]['ExcludeHours'][str(i)]
    if 'ExcludeHours' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',ExcludeHours'

def do_linear(cf,ds):
    level = ds.globalattributes['nc_level']
    for ThisOne in cf['Variables'].keys():
        if qcutils.haskey(cf,ThisOne,'Linear'):
            qcts.ApplyLinear(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'Drift'):
            qcts.ApplyLinearDrift(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'LocalDrift'):
            qcts.ApplyLinearDriftLocal(cf,ds,ThisOne)
    if 'do_linear' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_linear'

def do_rangecheck(cf,ds,section='',series='',code=2):
    '''Applies a range check to data series listed in the control file.  Data values that
       are less than the lower limit or greater than the upper limit are replaced with
       c.missing_value and the corresponding QC flag element is set to 2.'''
    if 'RangeCheck' not in cf[section][series].keys(): return
    if 'Lower' in cf[section][series]['RangeCheck'].keys():
        lwr = numpy.array(eval(cf[section][series]['RangeCheck']['Lower']))
        valid_lower = numpy.min(lwr)
        lwr = lwr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(c.missing_value))>c.eps)&
                                (ds.series[series]['Data']<lwr))
        ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_lower'] = cf[section][series]['RangeCheck']['Lower']
    if 'Upper' in cf[section][series]['RangeCheck'].keys():
        upr = numpy.array(eval(cf[section][series]['RangeCheck']['Upper']))
        valid_upper = numpy.min(upr)
        upr = upr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(c.missing_value))>c.eps)&
                                (ds.series[series]['Data']>upr))
        ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_upper'] = cf[section][series]['RangeCheck']['Upper']
        ds.series[series]['Attr']['valid_range'] = str(valid_lower)+','+str(valid_upper)
    if 'RangeCheck' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',RangeCheck'

def do_qcchecks(cf,ds,mode="verbose"):
    if "nc_level" in ds.globalattributes:
        level = str(ds.globalattributes["nc_level"])
        if mode!="quiet": log.info(" Doing the QC checks at level "+str(level))
    else:
        if mode!="quiet": log.info(" Doing the QC checks")
    # get the series list from the control file
    series_list = []
    for item in ["Variables","Drivers","Fluxes"]:
        if item in cf:
            section = item
            series_list = cf[item].keys()
    if len(series_list)==0:
        msg = " do_qcchecks: Variables, Drivers or Fluxes section not found in control file, skipping QC checks ..."
        log.warning(msg)
        return
    # loop over the series specified in the control file
    # first time for general QC checks
    for series in series_list:
        # check the series is in the data structure
        if series not in ds.series.keys():
            if mode!="quiet":
                msg = " do_qcchecks: series "+series+" not found in data structure, skipping ..."
                log.warning(msg)
            continue
        # if so, do the QC checks
        do_qcchecks_oneseries(cf,ds,section=section,series=series)
    # loop over the series in the control file
    # second time for dependencies
    for series in series_list:
        # check the series is in the data structure
        if series not in ds.series.keys():
            if mode!="quiet":
                msg = " do_qcchecks: series "+series+" not found in data structure, skipping ..."
                log.warning(msg)
            continue
        # if so, do dependency check
        do_dependencycheck(cf,ds,section=section,series=series,code=23,mode="quiet")

def do_qcchecks_oneseries(cf,ds,section='',series=''):
    if len(section)==0:
        section = qcutils.get_cfsection(cf,series=series,mode='quiet')
        if len(section)==0: return
    # do the range check
    do_rangecheck(cf,ds,section=section,series=series,code=2)
    # do the diurnal check
    do_diurnalcheck(cf,ds,section=section,series=series,code=5)
    # do exclude dates
    do_excludedates(cf,ds,section=section,series=series,code=6)
    # do exclude hours
    do_excludehours(cf,ds,section=section,series=series,code=7)
    if 'do_qcchecks' not in ds.globalattributes['Functions']:
        ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_qcchecks'
