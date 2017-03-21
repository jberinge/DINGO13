# Python modules
from configobj import ConfigObj
from collections import OrderedDict
import ast
import copy
import csv
import datetime
import dateutil
import logging
import netCDF4
import numpy
import ntpath
import os
import pdb
import platform
import sys
import time
import Tkinter, tkFileDialog
import xlrd
import xlwt
import xlsxwriter
# OzFluxQC modules
import cfg
import constants as c
import meteorologicalfunctions as mf
import qcck
import qcts
import qcutils

log = logging.getLogger('qc.io')

class DataStructure(object):
    def __init__(self):
        self.series = {}
        self.globalattributes = {}
        self.globalattributes["Functions"] = ""
        self.mergeserieslist = []
        self.averageserieslist = []
        self.returncodes = {}

def convert_v27tov28():
    """ Convert V2.7 (1D) netCDF files to V2.8 (3D). """
    # get the file names
    ncV27name = get_filename_dialog(path="../Sites")
    ncV28name = ncV27name.replace(".nc","_V28.nc")
    # read the V2.7 file
    ds = nc_read_series(ncV27name)
    # add the "time_zone" global attribute if it is not present
    if "time_zone" not in ds.globalattributes.keys():
        for gattr in ["site_name","SiteName"]:
            if gattr in ds.globalattributes.keys():
                time_zone,found = qcutils.get_timezone(ds.globalattributes[gattr],prompt="yes")
        ds.globalattributes["time_zone"] = time_zone
    # add the "missing_value" attribute if it is not present
    for ThisOne in ds.series.keys():
        if "missing_value" not in ds.series[ThisOne]["Attr"].keys():
            ds.series[ThisOne]["Attr"]["missing_value"] = numpy.int32(c.missing_value)
    # write the V2.8 file
    ncFile = nc_open_write(ncV28name,nctype='NETCDF4')
    nc_write_series(ncFile, ds)

def copy_datastructure(cf,ds_in):
    '''
    Return a copy of a data structure based on the following rules:
     1) if the netCDF file at the "copy_to" level does not exist
        then copy the existing data structure at the "input" level
        to create a new data structure at the "output" level.
    '''
    # assumptions that need to be checked are:
    #  - the start datetime of the two sets of data are the same
    #  - the end datetime of the L3 data is the same or after the
    #    end datetime of the the L4 data
    #    - if the end datetimes are the same then we are just re-processing something
    #    - if the end datetime for the L3 data is after the end date of the L4 data
    #      then more data has been added to this year and the user wants to gap fill
    #      the new data
    # modificatons to be made:
    #  - check the modification datetime of the L3 and L4 files:
    #     - if the L3 file is newer than the L4 file the disregard the "UseExistingOutFile" setting
    # get the output (L4) file name
    ct_filename = cf['Files']['file_path']+cf['Files']['out_filename']
    # if the L4 file does not exist then create the L4 data structure as a copy
    # of the L3 data structure
    if not os.path.exists(ct_filename):
        ds_out = copy.deepcopy(ds_in)
    # if the L4 file does exist ...
    if os.path.exists(ct_filename):
        # check to see if the user wants to use it
        if qcutils.get_keyvaluefromcf(cf,["Options"],"UseExistingOutFile",default="No")!='Yes':
            # if the user doesn't want to use the existing L4 data then create
            # the L4 data structure as a copy of the L3 data structure
            ds_out = copy.deepcopy(ds_in)
        else:
            # the user wants to use the data from an existing L4 file
            # get the netCDF file name at the "input" level
            outfilename = get_outfilenamefromcf(cf)
            # read the netCDF file at the "input" level
            ds_file = nc_read_series(outfilename)
            dt_file = ds_file.series['DateTime']['Data']
            sd_file = str(dt_file[0])
            ed_file = str(dt_file[-1])
            # create a copy of the data
            ds_out = copy.deepcopy(ds_in)
            dt_out = ds_out.series['DateTime']['Data']
            ts = ds_out.globalattributes['time_step']
            sd_out = str(dt_out[0])
            ed_out = str(dt_out[-1])
            # get the start and end indices based on the start and end dates
            si = qcutils.GetDateIndex(dt_out,sd_file,ts=ts,default=0,match='exact')
            ei = qcutils.GetDateIndex(dt_out,ed_file,ts=ts,default=-1,match='exact')
            # now replace parts of ds_out with the data read from file
            for ThisOne in ds_file.series.keys():
                # check to see if the L4 series exists in the L3 data
                if ThisOne in ds_out.series.keys():
                    # ds_out is the copy of the L3 data, now fill it with the L4 data read from file
                    ds_out.series[ThisOne]['Data'][si:ei+1] = ds_file.series[ThisOne]['Data']
                    ds_out.series[ThisOne]['Flag'][si:ei+1] = ds_file.series[ThisOne]['Flag']
                else:
                    # if it doesn't, create the series and put the data into it
                    ds_out.series[ThisOne] = {}
                    ds_out.series[ThisOne] = ds_file.series[ThisOne].copy()
                    # check to see if we have to append data to make the copy of the L4 data now
                    # in the L3 data structure the same length as the existing L3 data
                    nRecs_file = int(ds_file.globalattributes['nc_nrecs'])
                    nRecs_out = int(ds_out.globalattributes['nc_nrecs'])
                    if nRecs_file < nRecs_out:
                        # there is more data at L3 than at L4
                        # append missing data to make the series the same length
                        nRecs_append = nRecs_out - nRecs_file
                        data = numpy.array([c.missing_value]*nRecs_append,dtype=numpy.float64)
                        flag = numpy.ones(nRecs_append,dtype=numpy.int32)
                        ds_out.series[ThisOne]['Data'] = numpy.concatenate((ds_out.series[ThisOne]['Data'],data))
                        ds_out.series[ThisOne]['Flag'] = numpy.concatenate((ds_out.series[ThisOne]['Flag'],flag))
                    elif nRecs_file > nRecs_out:
                        # tell the user something is wrong
                        log.error('copy_datastructure: L3 file contains less data than L4 file')
                        # return an empty dictionary
                        ds_out = {}
                    else:
                        # nRecs_file and nRecs_out are equal so we do not need to do anything
                        pass
    return ds_out

def csv_read_series(cf):
    """
    Purpose:
     Reads a CSV file and returns the data in a data structure.
     The CSV file must conform to the following rules:
      1) the first line of the CSV file is a header line that contains
         the variable names for each column
      2) the second and all subsequent lines must contain data
      3) the first colume must contain a string representation
         of the datetime
      4) missing data in the CSV file is represented by a blank
         or by "NA".
    Usage:
     ds = csv_read_series(cf)
     where cf is a control file
           ds is a data structure
    Author: PRI
    Date: September 2015
    """
    ds = DataStructure()
    # get the filename
    csv_name = get_infilenamefromcf(cf)
    if len(csv_name)==0:
        log.error(' in_filename not found in control file')
        return ds
    if not os.path.exists(csv_name):
        log.error(' Input file '+csv_name+' specified in control file not found')
        return ds
    # get the header and first data rows
    FirstDataRow = int(qcutils.get_keyvaluefromcf(cf,["Files"],"in_firstdatarow")) - 1
    header_row = int(qcutils.get_keyvaluefromcf(cf,["Files"],"in_headerrow")) - 1
    # get the CSV workbook object
    log.info(" Opening and reading CSV file "+csv_name)
    # get the header line from the CSV file
    csv_file = open(csv_name,'rb')
    csv_reader = csv.reader(csv_file, delimiter=',')
    header = csv_reader.next()
    # construct the dtype list for use in numpy.genfromtxt
    dtype = [(header[0],'object')]+[(a,'f') for a in header[1:]]    
    # converter function for datetimes
    convert_datetime = lambda x: dateutil.parser.parse(x)
    # read the CSV file
    missing_values = "NA,NAN"
    data = numpy.genfromtxt(csv_name,delimiter=",",skip_header=0,names=True,converters={0:convert_datetime},
                            dtype=dtype,missing_values="NA",filling_values=float(-9999))
    # set some global attributes
    nRecs = len(data[header[0]])
    ds.globalattributes["nc_nrecs"] = str(nRecs)
    ds.globalattributes['featureType'] = 'timeseries'
    ds.globalattributes['csv_filename'] = csv_name
    ds.globalattributes['xl_datemode'] = str(0)
    s = os.stat(csv_name)
    t = time.localtime(s.st_mtime)
    ds.globalattributes['csv_moddatetime'] = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    # get the datetime
    ds.series["DateTime"] = {}
    ds.series["DateTime"]["Data"] = list(data["TIMESTAMP"])
    ds.series["DateTime"]["Flag"] = numpy.zeros(nRecs,dtype=numpy.int32)
    ds.series["DateTime"]["Attr"] = {}
    ds.series["DateTime"]["Attr"]["long_name"] = "Datetime in local timezone"
    ds.series["DateTime"]["Attr"]["units"] = "None"
    # get the variables
    var_list = cf["Variables"].keys()
    if "xlDateTime" in var_list: var_list.remove("xlDateTime")
    for var in var_list:
        if "csv" in cf["Variables"][var]:
            csv_var_name = cf["Variables"][var]["csv"]["name"]
        elif "xl" in cf["Variables"][var]:
            csv_var_name = cf["Variables"][var]["xl"]["name"]
        else:
            log.error(" No csv or xl section found in control file for "+var+", skipping ...")
            continue
        csv_var_name = csv_var_name.replace(".","")
        if csv_var_name not in data.dtype.names:
            log.warning("Requested variable "+csv_var_name+" not found in CSV file, skipping ...")
            continue
        ds.series[var] = {}
        ds.series[var]["Data"] = data[csv_var_name]
        ds.series[var]["Flag"] = numpy.zeros(nRecs,dtype=numpy.int32)
    return ds

def nc_2xls(ncfilename,outputlist=None):
    # read the netCDF file
    ds = nc_read_series(ncfilename,checktimestep=False)
    nRecs = int(ds.globalattributes["nc_nrecs"])
    nCols = len(ds.series.keys())
    if outputlist!=None: nCols = len(outputlist)
    # xlwt seems to only handle 225 columns
    if nRecs<65535 and nCols<220:
        # write the variables to the Excel 97/2003 file
        xlfilename= ncfilename.replace('.nc','.xls')
        xl_write_series(ds,xlfilename,outputlist=outputlist)
    else:
        # write the variables to the Excel 2010 file
        xlsxfilename= ncfilename.replace('.nc','.xlsx')
        xlsx_write_series(ds,xlsxfilename,outputlist=outputlist)

def read_eddypro_full(csvname):
    ds = DataStructure()
    csvfile = open(csvname,'rb')
    csvreader = csv.reader(csvfile)
    n = 0
    adatetime = []
    us_data_list = []
    us_flag_list = []
    Fh_data_list = []
    Fh_flag_list = []
    Fe_data_list = []
    Fe_flag_list = []
    Fc_data_list = []
    Fc_flag_list = []
    for row in csvreader:
        if n==0:
            header=row
        elif n==1:
            varlist=row
            us_data_col = varlist.index('u*')
            us_flag_col = varlist.index('qc_Tau')
            Fh_data_col = varlist.index('H')
            Fh_flag_col = varlist.index('qc_H')
            Fe_data_col = varlist.index('LE')
            Fe_flag_col = varlist.index('qc_LE')
            Fc_data_col = varlist.index('co2_flux')
            Fc_flag_col = varlist.index('qc_co2_flux')
        elif n==2:
            unitlist=row
        else:
            adatetime.append(datetime.datetime.strptime(row[1]+' '+row[2],'%Y-%m-%d %H:%M'))
            us_data_list.append(float(row[us_data_col]))
            us_flag_list.append(float(row[us_flag_col]))
            Fh_data_list.append(float(row[Fh_data_col]))
            Fh_flag_list.append(float(row[Fh_flag_col]))
            Fe_data_list.append(float(row[Fe_data_col]))
            Fe_flag_list.append(float(row[Fe_flag_col]))
            Fc_data_list.append(float(row[Fc_data_col]))
            Fc_flag_list.append(float(row[Fc_flag_col]))
        n = n + 1
    nRecs = len(adatetime)
    ds.series['DateTime'] = {}
    ds.series['DateTime']['Data'] = adatetime
    qcutils.round_datetime(ds,mode="nearest_timestep")
    ds.series['ustar'] = {}
    ds.series['ustar']['Data'] = numpy.array(us_data_list,dtype=numpy.float64)
    ds.series['ustar']['Flag'] = numpy.array(us_flag_list,dtype=numpy.int32)
    ds.series['Fh'] = {}
    ds.series['Fh']['Data'] = numpy.array(Fh_data_list,dtype=numpy.float64)
    ds.series['Fh']['Flag'] = numpy.array(Fh_flag_list,dtype=numpy.int32)
    ds.series['Fe'] = {}
    ds.series['Fe']['Data'] = numpy.array(Fe_data_list,dtype=numpy.float64)
    ds.series['Fe']['Flag'] = numpy.array(Fe_flag_list,dtype=numpy.int32)
    ds.series['Fc'] = {}
    ds.series['Fc']['Data'] = numpy.array(Fc_data_list,dtype=numpy.float64)
    ds.series['Fc']['Flag'] = numpy.array(Fc_flag_list,dtype=numpy.int32)
    ds.globalattributes["nc_nrecs"] = nRecs
    return ds

def smap_datetodatadictionary(ds,data_dict,nperday,ndays,si,ei):
    ldt = ds.series["DateTime"]["Data"][si:ei+1]
    # do the months
    month_1d,f,a = qcutils.GetSeries(ds,"Month",si=si,ei=ei)
    data_dict["Mo"] = {}
    data_dict["Mo"]["data"] = numpy.reshape(month_1d,[ndays,nperday])[:,0]
    data_dict["Mo"]["fmt"] = "0"
    # do the days
    data_dict["Day"] = {}
    day_1d,f,a = qcutils.GetSeries(ds,"Day",si=si,ei=ei)
    data_dict["Day"]["data"] = numpy.reshape(day_1d,[ndays,nperday])[:,0]
    data_dict["Day"]["fmt"] = "0"
    # day of the year
    data_dict["DOY"] = {}
    doy_1d = numpy.array([item.timetuple().tm_yday for item in ldt])
    data_dict["DOY"]["data"] = numpy.reshape(doy_1d,[ndays,nperday])[:,0]
    data_dict["DOY"]["fmt"] = "0"

def smap_docarbonfluxes(cf,ds,smap_label,si,ei):
    ncname = cf["Variables"][smap_label]["ncname"]
    data,flag,attr = qcutils.GetSeriesasMA(ds,ncname,si=si,ei=ei)
    data = data*12.01*1800/1E6
    data = numpy.ma.filled(data,float(-9999))
    return data,flag

def smap_donetshortwave(ds,smap_label,si,ei):
    ts = int(ds.globalattributes["time_step"])
    # do the net shortwave radiation
    Fsd,Fsd_flag,a = qcutils.GetSeriesasMA(ds,"Fsd",si=si,ei=ei)
    Fsu,Fsu_flag,a = qcutils.GetSeriesasMA(ds,"Fsu",si=si,ei=ei)
    # get the net shortwave radiation and convert to MJ/m2/day at the same time
    Fnsw = ((Fsd - Fsu)*ts*60)/1E6
    # now get the QC flag
    Fnsw_flag = Fsd_flag+Fsu_flag
    Fnsw = numpy.ma.filled(Fnsw,float(-9999))
    return Fnsw,Fnsw_flag

def smap_dopressure(ds,smap_label,si,ei):
    ps,ps_flag,attr = qcutils.GetSeriesasMA(ds,"ps",si=si,ei=ei)
    ps = ps/float(1000)
    ps = numpy.ma.filled(ps,float(-9999))
    return ps,ps_flag

def smap_doshortwave(ds,smap_label,si,ei):
    ts = int(ds.globalattributes["time_step"])
    Fsd,Fsd_flag,a = qcutils.GetSeriesasMA(ds,"Fsd",si=si,ei=ei)
    Fsd = (Fsd*ts*60)/1E6
    Fsd = numpy.ma.filled(Fsd,float(-9999))
    return Fsd,Fsd_flag

def smap_parseformat(fmt):
    if "." in fmt:
        numdec = len(fmt) - (fmt.index(".") + 1)
        strfmt = "{0:."+str(numdec)+"f}"
    else:
        strfmt = "{0:d}"
    return strfmt

def smap_qclabel(smap_label):
    if "_f" in smap_label:
        smap_qc_label=smap_label.replace("_f","_qc")
    else:
        smap_qc_label=smap_label+"_qc"
    return smap_qc_label

def smap_updatedatadictionary(cfvars,data_dict,data,flag,smap_label,nperday,ndays):
    data_dict[smap_label] = {}
    if cfvars[smap_label]["daily"].lower()=="sum":
        data_dict[smap_label]["data"] = numpy.ma.sum(numpy.ma.reshape(data,[ndays,nperday]),axis=1)
    elif cfvars[smap_label]["daily"].lower()=="average":
        data_dict[smap_label]["data"] = numpy.ma.average(numpy.ma.reshape(data,[ndays,nperday]),axis=1)
    elif cfvars[smap_label]["daily"].lower()=="skip":
        data_dict[smap_label]["data"] = numpy.reshape(data,[ndays,nperday])[:,0]
    else:
        print "smap_updatedatadictionary: unrecognised option for daily ("+str(cfvars[smap_label]["daily"])+")"
    data_dict[smap_label]["fmt"] = cfvars[smap_label]["format"]
    if cfvars[smap_label]["genqc"]=="True":
        smap_qc_label = smap_qclabel(smap_label)
        data_dict[smap_qc_label] = {}
        index_0 = numpy.where(flag==0)[0]
        index_not0 = numpy.where(flag>0)[0]
        flag[index_0] = numpy.int32(1)
        flag[index_not0] = numpy.int32(0)
        data_dict[smap_qc_label]["data"] = numpy.ma.sum(numpy.ma.reshape(flag,[ndays,nperday]),axis=1)/float(nperday)
        data_dict[smap_qc_label]["fmt"] = "0.00"

def smap_write_csv(cf):
    cfvars = cf["Variables"]
    smap_list = cfvars.keys()
    ncFileName = get_infilenamefromcf(cf)
    csvFileName_base = get_outfilenamefromcf(cf)    
    # read the netCDF file
    ds = nc_read_series(ncFileName)
    ts = int(ds.globalattributes["time_step"])
    nRecs = int(ds.globalattributes["nc_nrecs"])
    nperhr = int(float(60)/ts+0.5)
    nperday = int(float(24)*nperhr+0.5)
    dt = ds.series["DateTime"]["Data"]
    # get a list of years in the data file
    year_list = range(dt[0].year,dt[-1].year+1)
    years = numpy.array([item.year for item in dt])
    # loop over years in the data file
    data_dict = OrderedDict()
    for year in year_list:
        csvFileName = csvFileName_base+"_"+str(year)+"_SMAP.csv"
        # open the csv file
        csvfile = open(csvFileName,'wb')
        # write the header lines
        writer = smap_writeheaders(cf,csvfile)
        # get the start and end datetime
        year_index = numpy.where(years==year)[0]
        # add the last record from this year
        year_index = numpy.append(year_index,year_index[-1]+1)
        sdate = dt[max([0,year_index[0]])]
        edate = dt[min([year_index[-1],nRecs-1])]
        si = qcutils.GetDateIndex(dt,str(sdate),ts=ts,default=0,match="startnextday")
        ei = qcutils.GetDateIndex(dt,str(edate),ts=ts,default=nRecs-1,match="endpreviousday")
        data_dict["DateTime"] = dt[si:ei+1]
        log.info(" Writing "+str(data_dict["DateTime"][0])+" to "+ str(data_dict["DateTime"][-1]))
        ndays = len(data_dict["DateTime"])/nperday
        # put the month, day and DOY into the data dictionary
        smap_datetodatadictionary(ds,data_dict,nperday,ndays,si,ei)
        # first column in SMAP csv file will be the SMAP ID number
        smap_id = numpy.array([cf["General"]["SMAP_ID"]]*ndays)
        # loop over the data required, massage units if necessary and put the data into a dictionary for later use
        smap_list = ["Rn_f","Rs_f","PAR_f","Ta","VPD","Ts_f","PREC","SWC","NEE","GPP","Reco","PRESS","SNOWD"]
        for smap_label in smap_list:
            if smap_label in ["Mo","Day","DOY"]: continue
            if smap_label=="Rn_f":
                data,flag = smap_donetshortwave(ds,smap_label,si,ei)
            elif smap_label=="Rs_f":
                data,flag = smap_doshortwave(ds,smap_label,si,ei)
            elif smap_label=="PAR_f" or smap_label=="SNOWD":
                data = numpy.array([-9999]*len(data_dict["DateTime"]))
                flag = numpy.array([1]*len(data_dict["DateTime"]))
                cfvars[smap_label]["daily"] = "skip"
            elif smap_label=="PRESS":
                data,flag = smap_dopressure(ds,smap_label,si,ei)
            elif smap_label in ["GPP","NEE","Reco"]:
                data,flag = smap_docarbonfluxes(cf,ds,smap_label,si,ei)
            else:
                data,flag,attr = qcutils.GetSeries(ds,cfvars[smap_label]["ncname"],si=si,ei=ei)
            smap_updatedatadictionary(cfvars,data_dict,data,flag,smap_label,nperday,ndays)
        # now loop over the days and write the data out
        for i in range(ndays):
            data_list = []
            data_list.append(smap_id[i])
            for smap_label in data_dict.keys():
                if smap_label=="DateTime": continue
                strfmt = smap_parseformat(data_dict[smap_label]["fmt"])
                if "d" in strfmt:
                    data_list.append(strfmt.format(int(round(data_dict[smap_label]["data"][i]))))
                else:
                    data_list.append(strfmt.format(data_dict[smap_label]["data"][i]))
            writer.writerow(data_list)
        csvfile.close()

def smap_writeheaders(cf,csvfile):
    writer = csv.writer(csvfile)
    # write the header lines to the csv file
    series_list = cf["Variables"].keys()
    for item in cf["General"]:
        if item in ["SMAP_ID"]: continue
        writer.writerow([item,str(cf['General'][item])])
    # write the units and variable name header lines to the csv file
    units_list = ["-","-","-","-"]
    row_list = ['ID','Mo','Day','DOY']
    for smap_label in series_list:
        row_list.append(smap_label)
        units_list.append(cf["Variables"][smap_label]["units"])
        if cf["Variables"][smap_label]["genqc"]=="True":
            smap_qc_label = smap_qclabel(smap_label)
            row_list.append(smap_qc_label)
            units_list.append("-")
    writer.writerow(units_list)
    writer.writerow(row_list)
    return writer

def xl2nc(cf,InLevel):
    # get the data series from the Excel file
    in_filename = get_infilenamefromcf(cf)
    if not qcutils.file_exists(in_filename,mode="quiet"):
        msg = " Input file "+in_filename+" not found ..."
        log.error(msg)
        return
    file_name,file_extension = os.path.splitext(in_filename)
    if "csv" in file_extension.lower():
        ds = csv_read_series(cf)
        # get a series of Excel datetime from the Python datetime objects
        qcutils.get_xldatefromdatetime(ds)
    else:
        ds = xl_read_series(cf)
        # get a series of Python datetime objects from the Excel datetime
        qcutils.get_datetimefromxldate(ds)
    #if len(ds.series.keys())==0: return 1
    # get the netCDF attributes from the control file
    qcts.do_attributes(cf,ds)
    # round the Python datetime to the nearest second
    qcutils.round_datetime(ds,mode="nearest_second")
    #check for gaps in the Python datetime series and fix if present
    fixtimestepmethod = qcutils.get_keyvaluefromcf(cf,["options"],"FixTimeStepMethod",default="round")
    if qcutils.CheckTimeStep(ds): qcutils.FixTimeStep(ds,fixtimestepmethod=fixtimestepmethod)
    # recalculate the Excel datetime
    qcutils.get_xldatefromdatetime(ds)
    # get the Year, Month, Day etc from the Python datetime
    qcutils.get_ymdhmsfromdatetime(ds)
    # write the processing level to a global attribute
    ds.globalattributes['nc_level'] = str(InLevel)
    # get the start and end date from the datetime series unless they were
    # given in the control file
    if 'start_date' not in ds.globalattributes.keys():
        ds.globalattributes['start_date'] = str(ds.series['DateTime']['Data'][0])
    if 'end_date' not in ds.globalattributes.keys():
        ds.globalattributes['end_date'] = str(ds.series['DateTime']['Data'][-1])
    # calculate variances from standard deviations and vice versa
    qcts.CalculateStandardDeviations(cf,ds)
    # create new variables using user defined functions
    qcts.DoFunctions(cf,ds)
    # create a series of synthetic downwelling shortwave radiation
    qcts.get_synthetic_fsd(ds)
    # write the data to the netCDF file
    outfilename = get_outfilenamefromcf(cf)
    ncFile = nc_open_write(outfilename)
    nc_write_series(ncFile,ds)
    return 0

def fn_write_csv(cf):
    # get the file names
    ncFileName = get_infilenamefromcf(cf)
    csvFileName = get_outfilenamefromcf(cf)
    # open the csv file
    csvfile = open(csvFileName,'wb')
    writer = csv.writer(csvfile)
    # read the netCDF file
    ds = nc_read_series(ncFileName)
    # Tumbarumba doesn't have RH in the netCDF files
    if "RH" not in ds.series.keys():
        Ah,f,a = qcutils.GetSeriesasMA(ds,'Ah')
        Ta,f,a = qcutils.GetSeriesasMA(ds,'Ta')
        RH = mf.RHfromabsolutehumidity(Ah, Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='relative_humidity')
        qcutils.CreateSeries(ds,"RH",RH,FList=['Ta','Ah'],Attr=attr)
    ts = int(ds.globalattributes["time_step"])
    ts_delta = datetime.timedelta(minutes=ts)
    # get the datetime series
    dt = ds.series["DateTime"]["Data"]
    # check the start datetime of the series and adjust if necessary
    start_datetime = dateutil.parser.parse(str(cf["General"]["start_datetime"]))
    if dt[0]<start_datetime:
        # requested start_datetime is after the start of the file
        log.info(" Truncating start of file")
        si = qcutils.GetDateIndex(dt,str(start_datetime),ts=ts,match="exact")
        for thisone in ds.series.keys():
            ds.series[thisone]["Data"] = ds.series[thisone]["Data"][si:]
            ds.series[thisone]["Flag"] = ds.series[thisone]["Flag"][si:]
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
    elif dt[0]>start_datetime:
        # requested start_datetime is before the start of the file
        log.info(" Padding start of file")
        dt_patched = [ldt for ldt in qcutils.perdelta(start_datetime, dt[0]-ts_delta, ts_delta)]
        data_patched = numpy.ones(len(dt_patched))*float(c.missing_value)
        flag_patched = numpy.ones(len(dt_patched))
        # list of series in the data structure
        series_list = ds.series.keys()
        # ds.series["DateTime"]["Data"] is a list not a numpy array so we must treat it differently
        ds.series["DateTime"]["Data"] = dt_patched+ds.series["DateTime"]["Data"]
        ds.series["DateTime"]["Flag"] = numpy.concatenate((flag_patched,ds.series["DateTime"]["Flag"]))
        series_list.remove("DateTime")
        for thisone in series_list:
            ds.series[thisone]["Data"] = numpy.concatenate((data_patched,ds.series[thisone]["Data"]))
            ds.series[thisone]["Flag"] = numpy.concatenate((flag_patched,ds.series[thisone]["Flag"]))
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
        # refresh the year, month, day etc arrays now that we have padded the datetime series
        qcutils.get_ymdhmsfromdatetime(ds)
    # now check the end datetime of the file
    end_datetime = dateutil.parser.parse(str(cf["General"]["end_datetime"]))
    if dt[-1]>end_datetime:
        # requested end_datetime is before the end of the file
        msg = " Truncating end of file "+dt[-1].strftime("%Y-%m-%d %H:%M")+" "+end_datetime.strftime("%Y-%m-%d %H:%M")
        log.info(msg)
        ei = qcutils.GetDateIndex(dt,str(end_datetime),ts=ts,match="exact")
        for thisone in ds.series.keys():
            ds.series[thisone]["Data"] = ds.series[thisone]["Data"][:ei+1]
            ds.series[thisone]["Flag"] = ds.series[thisone]["Flag"][:ei+1]
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
    elif dt[-1]<end_datetime:
        # requested end_datetime is before the requested end date
        msg = " Padding end of file "+dt[-1].strftime("%Y-%m-%d %H:%M")+" "+end_datetime.strftime("%Y-%m-%d %H:%M")
        log.info(msg)
        dt_patched = [ldt for ldt in qcutils.perdelta(dt[-1]+ts_delta, end_datetime, ts_delta)]
        data_patched = numpy.ones(len(dt_patched))*float(c.missing_value)
        flag_patched = numpy.ones(len(dt_patched))
        # list of series in the data structure
        series_list = ds.series.keys()
        # ds.series["DateTime"]["Data"] is a list not a numpy array so we must treat it differently
        ds.series["DateTime"]["Data"] = ds.series["DateTime"]["Data"]+dt_patched
        ds.series["DateTime"]["Flag"] = numpy.concatenate((ds.series["DateTime"]["Flag"],flag_patched))
        series_list.remove("DateTime")
        for thisone in series_list:
            ds.series[thisone]["Data"] = numpy.concatenate((ds.series[thisone]["Data"],data_patched))
            ds.series[thisone]["Flag"] = numpy.concatenate((ds.series[thisone]["Flag"],flag_patched))
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
        # refresh the year, month, day etc arrays now that we have padded the datetime series
        qcutils.get_ymdhmsfromdatetime(ds)
    if ts==30:
        nRecs_year = 17520
        nRecs_leapyear = 17568
    elif ts==60:
        nRecs_year = 8760
        nRecs_leapyear = 8784
    else:
        log.error(" Unrecognised time step ("+str(ts)+")")
        return
    if (int(ds.globalattributes["nc_nrecs"])!=nRecs_year) & (int(ds.globalattributes["nc_nrecs"])!=nRecs_leapyear):
        log.error(" Number of records in file does not equal "+str(nRecs_year)+" or "+str(nRecs_leapyear))
        msg = str(len(ds.series["DateTime"]["Data"]))+" "+str(ds.series["DateTime"]["Data"][0])
        msg = msg+" "+str(ds.series["DateTime"]["Data"][-1])
        log.error(msg)
        return
    # get the date and time data
    Day,flag,attr = qcutils.GetSeries(ds,'Day')
    Month,flag,attr = qcutils.GetSeries(ds,'Month')
    Year,flag,attr = qcutils.GetSeries(ds,'Year')
    Hour,flag,attr = qcutils.GetSeries(ds,'Hour')
    Minute,flag,attr = qcutils.GetSeries(ds,'Minute')
    # get the data
    data = {}
    series_list = cf["Variables"].keys()
    for series in series_list:
        ncname = cf["Variables"][series]["ncname"]
        if ncname not in ds.series.keys():
            log.error("Series "+ncname+" not in netCDF file, skipping ...")
            series_list.remove(series)
            continue
        data[series] = ds.series[ncname]
        fmt = cf["Variables"][series]["format"]
        if "." in fmt:
            numdec = len(fmt) - (fmt.index(".") + 1)
            strfmt = "{0:."+str(numdec)+"f}"
        else:
            strfmt = "{0:d}"
        data[series]["fmt"] = strfmt
    #adjust units if required
    for series in series_list:
        if series=="FC" and data[series]["Attr"]["units"]=='mg/m2/s':
            data[series]["Data"] = mf.Fc_umolpm2psfrommgpm2ps(data[series]["Data"])
            data[series]["Attr"]["units"] = "umol/m2/s"
        if series=="CO2" and data[series]["Attr"]["units"]=='mg/m3':
            CO2 = data["CO2"]["Data"]
            TA = data["TA"]["Data"]
            PA = data["PA"]["Data"]
            data[series]["Data"] = mf.co2_ppmfrommgpm3(CO2,TA,PA)
            data[series]["Attr"]["units"] = "umol/mol"
        if series=="H2O" and data[series]["Attr"]["units"]=='g/m3':
            H2O = data["H2O"]["Data"]
            TA = data["TA"]["Data"]
            PA = data["PA"]["Data"]
            data[series]["Data"] = mf.h2o_mmolpmolfromgpm3(H2O,TA,PA)
            data[series]["Attr"]["units"] = "mmol/mol"
        if series=="RH" and data[series]["Attr"]["units"] in ["fraction","frac"]:
            data[series]["Data"] = float(100)*data[series]["Data"]
            data[series]["Attr"]["units"] = "%"
    # write the general information to csv file
    for item in cf["General"]:
        writer.writerow([item,str(cf['General'][item])])
    # write the variable names to the csv file
    row_list = ['DateTime','Year','Month','Day','HHMM']
    for item in series_list:
        row_list.append(item)
    writer.writerow(row_list)
    # write the units line to the csv file
    units_list = ["-","-","-","-","-"]
    for item in series_list:
        units_list.append(data[item]["Attr"]["units"])
    writer.writerow(units_list)
    # now write the data
    for i in range(len(Year)):
        # get the datetime string
        dtstr = '%02d/%02d/%d %02d:%02d'%(Day[i],Month[i],Year[i],Hour[i],Minute[i])
        hrmn = '%02d%02d'%(Hour[i],Minute[i])
        dttup = datetime.datetime(Year[i],Month[i],Day[i],Hour[i],Minute[i]).timetuple()
        doy = float(dttup.tm_yday) + float(dttup.tm_hour)/24 + float(dttup.tm_min)/1440
        data_list = [dtstr,'%d'%(Year[i]),'%02d'%(Month[i]),'%02d'%(Day[i]),hrmn]
        for series in series_list:
            strfmt = data[series]["fmt"]
            if "d" in strfmt:
                data_list.append(strfmt.format(int(round(data[series]["Data"][i]))))
            else:
                data_list.append(strfmt.format(data[series]["Data"][i]))
        writer.writerow(data_list)
    # close the csv file
    csvfile.close()
    return

def get_controlfilecontents(ControlFileName,mode="verbose"):
    if mode!="quiet": log.info(' Processing the control file ')
    if len(ControlFileName)!=0:
        cf = ConfigObj(ControlFileName)
        cf['controlfile_name'] = ControlFileName
    else:
        cf = ConfigObj()
    if "Files" in cf:
        if "plot_path" not in cf["Files"].keys():
            cf["Files"]["plot_path"] = "plots/"
    return cf

def get_controlfilename(path='.',title='Choose a control file'):
    log.info(' Choosing the control file ')
    root = Tkinter.Tk(); root.withdraw()
    name = tkFileDialog.askopenfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    return name

def get_ncdtype(Series):
    sd = Series.dtype.name
    dt = 'f'
    if sd=='float64': dt = 'd'
    if sd=='int32': dt = 'i'
    if sd=='int64': dt = 'l'
    return dt

def get_filename_dialog(path='.',title='Choose a file'):
    """
    Purpose:
     Put up a file open dialog and let the user browse to open a file
    Usage:
     fname = qcio.get_filename_dialog(path=<path_to_file>,title=<tile>)
     where path  - the path to the file location, optional
           title - the title for the file open dialog, optional
    Returns:
     fname - the full file name including path, string
    Author: PRI
    Date: Back in the day
    """
    root = Tkinter.Tk(); root.withdraw()
    FileName = tkFileDialog.askopenfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    return str(FileName)

def get_infilenamefromcf(cf):
    path = qcutils.get_keyvaluefromcf(cf,["Files"],"file_path",default="")
    name = qcutils.get_keyvaluefromcf(cf,["Files"],"in_filename",default="")
    return str(path)+str(name)

def get_outfilenamefromcf(cf):
    path = qcutils.get_keyvaluefromcf(cf,["Files"],"file_path",default="")
    name = qcutils.get_keyvaluefromcf(cf,["Files"],"out_filename",default="")
    return str(path)+str(name)

def get_outputlistfromcf(cf,filetype):
    try:
        outputlist = ast.literal_eval(cf['Output'][filetype])
    except:
        #log.info('get_outputlistfromcf: Unable to get output list from Output section in control file')
        outputlist = None
    return outputlist

def get_seriesstats(cf,ds):
    # open an Excel file for the flag statistics
    level = ds.globalattributes['nc_level']
    out_filename = get_outfilenamefromcf(cf)
    xl_filename = out_filename.replace('.nc','_FlagStats.xls')
    log.info(' Writing flag stats to Excel file '+xl_filename)
    xlFile = xlwt.Workbook()
    xlFlagSheet = xlFile.add_sheet('Flag')
    # get the flag statistics
    #xlRow = 0
    #xlCol = 0
    #xlFlagSheet.write(xlRow,xlCol,'0:')
    #xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag00'])
    #xlFlagSheet.write(xlRow,xlCol+2,'1:')
    #xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag01'])
    #xlFlagSheet.write(xlRow,xlCol+4,'2:')
    #xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag02'])
    #xlFlagSheet.write(xlRow,xlCol+6,'3:')
    #xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag03'])
    #xlFlagSheet.write(xlRow,xlCol+8,'4:')
    #xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag04'])
    #xlFlagSheet.write(xlRow,xlCol+10,'5:')
    #xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag05'])
    #xlFlagSheet.write(xlRow,xlCol+12,'6:')
    #xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag06'])
    #xlFlagSheet.write(xlRow,xlCol+14,'7:')
    #xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag07'])
    #xlRow = xlRow + 1
    #xlFlagSheet.write(xlRow,xlCol,'10:')
    #xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag10'])
    #xlFlagSheet.write(xlRow,xlCol+2,'11:')
    #xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag11'])
    #xlFlagSheet.write(xlRow,xlCol+4,'12:')
    #xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag12'])
    #xlFlagSheet.write(xlRow,xlCol+6,'13:')
    #xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag13'])
    #xlFlagSheet.write(xlRow,xlCol+8,'14:')
    #xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag14'])
    #xlFlagSheet.write(xlRow,xlCol+10,'15:')
    #xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag15'])
    #xlFlagSheet.write(xlRow,xlCol+12,'16:')
    #xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag16'])
    #xlFlagSheet.write(xlRow,xlCol+14,'17:')
    #xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag17'])
    #xlFlagSheet.write(xlRow,xlCol+16,'18:')
    #xlFlagSheet.write(xlRow,xlCol+17,ds.globalattributes['Flag18'])
    #xlFlagSheet.write(xlRow,xlCol+18,'19:')
    #xlFlagSheet.write(xlRow,xlCol+19,ds.globalattributes['Flag19'])
    #xlRow = xlRow + 1
    #xlFlagSheet.write(xlRow,xlCol,'30:')
    #xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag30'])
    #xlFlagSheet.write(xlRow,xlCol+2,'31:')
    #xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag31'])
    #xlFlagSheet.write(xlRow,xlCol+4,'32:')
    #xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag32'])
    #xlFlagSheet.write(xlRow,xlCol+6,'33:')
    #xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag33'])
    #xlFlagSheet.write(xlRow,xlCol+8,'34:')
    #xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag34'])
    #xlFlagSheet.write(xlRow,xlCol+10,'35:')
    #xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag35'])
    #xlFlagSheet.write(xlRow,xlCol+12,'36:')
    #xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag36'])
    #xlFlagSheet.write(xlRow,xlCol+14,'37:')
    #xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag37'])
    #xlFlagSheet.write(xlRow,xlCol+16,'38:')
    #xlFlagSheet.write(xlRow,xlCol+17,ds.globalattributes['Flag38'])
    #xlFlagSheet.write(xlRow,xlCol+18,'39:')
    #xlFlagSheet.write(xlRow,xlCol+19,ds.globalattributes['Flag39'])
    bins = numpy.arange(-0.5,23.5)
    xlRow = 5
    xlCol = 1
    for Value in bins[:len(bins)-1]:
        xlFlagSheet.write(xlRow,xlCol,int(Value+0.5))
        xlCol = xlCol + 1
    xlRow = xlRow + 1
    xlCol = 0
    dsVarNames = ds.series.keys()
    dsVarNames.sort(key=unicode.lower)
    for ThisOne in dsVarNames:
        data,flag,attr = qcutils.GetSeries(ds, ThisOne)
        hist, bin_edges = numpy.histogram(flag, bins=bins)
        xlFlagSheet.write(xlRow,xlCol,ThisOne)
        xlCol = xlCol + 1
        for Value in hist:
            xlFlagSheet.write(xlRow,xlCol,float(Value))
            xlCol = xlCol + 1
        xlCol = 0
        xlRow = xlRow + 1
    xlFile.save(xl_filename)

def load_controlfile(path='.',title='Choose a control file'):
    ''' 
    Returns a control file object.
    USAGE: cf = load_controlfile([path=<some_path_to_a_controlfile>])
    The "path" keyword is optional.
    '''
    name = get_controlfilename(path=path,title=title)
    cf = get_controlfilecontents(name)
    return cf

def nc_concatenate(cf):
    # initialise logicals
    TimeGap = False
    # get an instance of the data structure
    ds = DataStructure()
    # get the input file list
    InFile_list = cf['Files']['In'].keys()
    # read in the first file
    ncFileName = cf['Files']['In'][InFile_list[0]]
    log.info(' Reading data from '+ncFileName)
    fixtimestepmethod = qcutils.get_keyvaluefromcf(cf,["Options"],"FixTimeStepMethod",default="round")
    ds_n = nc_read_series(ncFileName,fixtimestepmethod=fixtimestepmethod)
    if len(ds_n.series.keys())==0:
        log.error(' An error occurred reading netCDF file: '+ncFileName)
        return
    # fill the global attributes
    for ThisOne in ds_n.globalattributes.keys():
        ds.globalattributes[ThisOne] = ds_n.globalattributes[ThisOne]
    # check that we have 'Ws' and 'Wd' series
    if "Ws" not in ds_n.series.keys():
        if "Ws_CSAT" in ds_n.series.keys():
            msg = "Ws not found, copying series Ws_CSAT to Ws"
            log.info(msg)
            ds_n.series["Ws"] = ds_n.series["Ws_CSAT"].copy()
        else:
            msg = "Both Ws and Ws_CSAT missing from file"
            log.warning(msg)
    if "Wd" not in ds_n.series.keys():
        if "Wd_CSAT" in ds_n.series.keys():
            msg = "Wd not found, copying series Wd_CSAT to Wd"
            log.info(msg)
            ds_n.series["Wd"] = ds_n.series["Wd_CSAT"].copy()
        else:
            msg = "Both Wd and Wd_CSAT missing from file"
            log.warning(msg)
    # fill the variables
    for ThisOne in ds_n.series.keys():
        if ThisOne=="Fc":
            Fc,flag,attr = qcutils.GetSeriesasMA(ds_n, ThisOne)
            if attr['units']=='mg/m2/s':
                log.info("Converting Fc to umol/m2/s")
                Fc = mf.Fc_umolpm2psfrommgpm2ps(Fc)
                attr['units'] = 'umol/m2/s'
                attr['standard_name'] = 'surface_upward_mole_flux_of_carbon_dioxide'
                qcutils.CreateSeries(ds_n,ThisOne,Fc,Flag=flag,Attr=attr)
        ds.series[ThisOne] = {}
        ds.series[ThisOne]['Data'] = ds_n.series[ThisOne]['Data']
        ds.series[ThisOne]['Flag'] = ds_n.series[ThisOne]['Flag']
        ds.series[ThisOne]['Attr'] = {}
        for attr in ds_n.series[ThisOne]['Attr'].keys():
            ds.series[ThisOne]['Attr'][attr] = ds_n.series[ThisOne]['Attr'][attr]
    ts = int(ds.globalattributes['time_step'])
    # loop over the remaining files given in the control file
    for n in InFile_list[1:]:
        ncFileName = cf['Files']['In'][InFile_list[int(n)]]
        log.info(' Reading data from '+ncFileName)
        #print 'ncconcat: reading data from '+ncFileName
        ds_n = nc_read_series(ncFileName,fixtimestepmethod=fixtimestepmethod)
        if len(ds.series.keys())==0:
            log.error(' An error occurred reading the netCDF file: '+ncFileName)
            return
        dt_n = ds_n.series['DateTime']['Data']
        dt = ds.series['DateTime']['Data']
        nRecs_n = len(ds_n.series['xlDateTime']['Data'])
        nRecs = len(ds.series['xlDateTime']['Data'])
        # check that we have 'Ws' and 'Wd' series
        if "Ws" not in ds_n.series.keys():
            if "Ws_CSAT" in ds_n.series.keys():
                msg = "Ws not found, copying series Ws_CSAT to Ws"
                log.info(msg)
                ds_n.series["Ws"] = ds_n.series["Ws_CSAT"].copy()
            else:
                msg = "Both Ws and Ws_CSAT missing from file"
                log.warning(msg)
        if "Wd" not in ds_n.series.keys():
            if "Wd_CSAT" in ds_n.series.keys():
                msg = "Wd not found, copying series Wd_CSAT to Wd"
                log.info(msg)
                ds_n.series["Wd"] = ds_n.series["Wd_CSAT"].copy()
            else:
                msg = "Both Wd and Wd_CSAT missing from file"
                log.warning(msg)
        #print ds.series['DateTime']['Data'][-1],ds_n.series['DateTime']['Data'][-1]
        #print dt[-1],dt[-1]+datetime.timedelta(minutes=ts),dt_n[0]
        if dt_n[0]<dt[-1]+datetime.timedelta(minutes=ts):
            log.info(' Overlapping times detected in consecutive files')
            si = qcutils.GetDateIndex(dt_n,str(dt[-1]),ts=ts)+1
            ei = -1
        if dt_n[0]==dt[-1]+datetime.timedelta(minutes=ts):
            log.info(' Start and end times OK in consecutive files')
            si = 0; ei = -1
        if dt_n[0]>dt[-1]+datetime.timedelta(minutes=ts):
            log.info(' Gap between start and end times in consecutive files')
            si = 0; ei = -1
            #TimeGap = True
        # loop over the data series in the concatenated file
        for ThisOne in ds.series.keys():
            if ThisOne=="Fc":
                Fc,flag,attr = qcutils.GetSeriesasMA(ds_n,ThisOne)
                if attr['units']=='mg/m2/s':
                    log.info(" Converting Fc to umol/m2/s")
                    Fc = mf.Fc_umolpm2psfrommgpm2ps(Fc)
                    attr['units'] = 'umol/m2/s'
                    attr['standard_name'] = 'surface_upward_mole_flux_of_carbon_dioxide'
                    qcutils.CreateSeries(ds_n,ThisOne,Fc,Flag=flag,Attr=attr)
            # does this series exist in the file being added to the concatenated file
            if ThisOne in ds_n.series.keys():
                # if so, then append this series to the concatenated series
                if type(ds.series[ThisOne]["Data"]) is list:
                    ds.series[ThisOne]['Data'] = ds.series[ThisOne]['Data']+ds_n.series[ThisOne]['Data'][si:ei]
                else:
                    ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
                ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
            else:
                # if not, then create a dummy series and concatenate that
                ds_n.series[ThisOne] = {}
                ds_n.series[ThisOne]['Data'] = numpy.array([c.missing_value]*nRecs_n,dtype=numpy.float64)
                ds_n.series[ThisOne]['Flag'] = numpy.array([1]*nRecs_n,dtype=numpy.int32)
                ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
                ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
        # and now loop over the series in the file being concatenated
        for ThisOne in ds_n.series.keys():
            # does this series exist in the concatenated data
            if ThisOne not in ds.series.keys():
                # if not then add it
                ds.series[ThisOne] = {}
                ds.series[ThisOne]['Data'] = numpy.array([c.missing_value]*nRecs,dtype=numpy.float64)
                ds.series[ThisOne]['Flag'] = numpy.array([1]*nRecs,dtype=numpy.int32)
                ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
                ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
                ds.series[ThisOne]['Attr'] = {}
                for attr in ds_n.series[ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = ds_n.series[ThisOne]['Attr'][attr]
    # update the number of records
    ds.globalattributes["nc_nrecs"] = len(ds.series["DateTime"]["Data"])
    # now sort out any time gaps
    if qcutils.CheckTimeStep(ds):
        fixtimestepmethod = qcutils.get_keyvaluefromcf(cf,["Options"],"FixTimeStepMethod",default="round")
        qcutils.FixTimeStep(ds,fixtimestepmethod=fixtimestepmethod)
        # update the Excel datetime from the Python datetime
        qcutils.get_xldatefromdatetime(ds)
        # update the Year, Month, Day etc from the Python datetime
        qcutils.get_ymdhmsfromdatetime(ds)
    # if requested, fill any small gaps by interpolation
    # get a list of series in ds excluding the QC flags
    series_list = [item for item in ds.series.keys() if "_QCFlag" not in item]
    # remove the datetime variables, these will have no gaps
    datetime_list = ["xlDateTime","DateTime","Year","Month","Day","Hour","Minute","Second","Hdh","Ddd"]
    for item in datetime_list:
        if item in series_list: series_list.remove(item)
    # loop over the non-datetime data series in ds and interpolate
    # get the maximum gap length (in hours) from the control file
    maxlen = int(qcutils.get_keyvaluefromcf(cf,["Options"],"MaxGapInterpolate",default=3))
    # now loop over the series and do the interpolation
    log.info(" Interpolating over fixed time gaps ("+str(maxlen)+" hour max)")
    for item in series_list:
        qcts.InterpolateOverMissing(ds,series=item,maxlen=maxlen)
    # make sure we have all of the humidities
    qcts.CalculateHumidities(ds)
    # and make sure we have all of the meteorological variables
    qcts.CalculateMeteorologicalVariables(ds)
    # re-calculate the synthetic Fsd
    qcts.get_synthetic_fsd(ds)
    # re-apply the quality control checks (range, diurnal and rules)
    qcck.do_qcchecks(cf,ds)
    # update the coverage statistics
    qcutils.get_coverage_individual(ds)
    # write the netCDF file
    ncFileName = qcutils.get_keyvaluefromcf(cf,["Files","Out"],"ncFileName",default="out.nc")
    log.info(' Writing data to '+ncFileName)
    ncFile = nc_open_write(ncFileName)
    ndims = qcutils.get_keyvaluefromcf(cf,["Options"],"NumberOfDimensions", default=3)
    nc_write_series(ncFile,ds,ndims=ndims)

def nc_split():
    split_info = {}
    split_gui = Tkinter.Toplevel()
    split_gui.wm_title("Split netCDF file")
    split_gui.grid()
    # first row contains the input file name selection
    nrow = 0
    split_gui.infilenameLabel = Tkinter.Label(split_gui,text="Input file")
    split_gui.infilenameLabel.grid(row=nrow,column=0,columnspan=1)
    split_gui.infilename = Tkinter.StringVar()
    split_gui.infilename.set("")
    split_gui.infilenameEntry = Tkinter.Entry(split_gui,textvariable=split_gui.infilename,width=30)
    split_gui.infilenameEntry.grid(row=nrow,column=1)
    split_gui.infilenameBrowse = Tkinter.Button(split_gui,text="Browse",command=lambda:ncsplit_infilename_browse(split_gui))
    split_gui.infilenameBrowse.grid(row=nrow,column=2)
    # second row has the start date entry
    nrow = nrow + 1
    split_gui.startLabel = Tkinter.Label(split_gui, text="Start date (YYYY-MM-DD)")
    split_gui.startLabel.grid(row=nrow,column=0,columnspan=1)
    split_gui.startEntry = Tkinter.Entry(split_gui)
    split_gui.startEntry.grid(row=nrow,column=1,columnspan=1)
    # third row has the end date entry
    nrow = nrow + 1
    split_gui.endLabel = Tkinter.Label(split_gui, text="End date   (YYYY-MM-DD)")
    split_gui.endLabel.grid(row=nrow,column=0,columnspan=1)
    split_gui.endEntry = Tkinter.Entry(split_gui)
    split_gui.endEntry.grid(row=nrow,column=1,columnspan=1)
    # fourth row contains the output file name selection
    nrow = nrow + 1
    split_gui.outfilenameLabel = Tkinter.Label(split_gui,text="Output file")
    split_gui.outfilenameLabel.grid(row=nrow,column=0,columnspan=1)
    split_gui.outfilename = Tkinter.StringVar()
    split_gui.outfilename.set("")
    split_gui.outfilenameEntry = Tkinter.Entry(split_gui,textvariable=split_gui.outfilename,width=30)
    split_gui.outfilenameEntry.grid(row=nrow,column=1)
    split_gui.outfilenameBrowse = Tkinter.Button(split_gui,text="Browse",command=lambda:ncsplit_outfilename_browse(split_gui))
    split_gui.outfilenameBrowse.grid(row=nrow,column=2)
    # action buttons on the bottom row
    nrow = nrow + 1
    split_gui.doneButton = Tkinter.Button(split_gui,text="Done",command=lambda:ncsplit_done(split_gui))
    split_gui.doneButton.grid(row=nrow,column=0,columnspan=1)
    split_gui.runButton = Tkinter.Button(split_gui,text="Run",command=lambda:ncsplit_run(split_gui))
    split_gui.runButton.grid(row=nrow,column=1,columnspan=1)
    # progress message area
    nrow = nrow + 1
    split_gui.progress_row = nrow
    split_gui.progress = Tkinter.Label(split_gui, text='Waiting for input ...')
    split_gui.progress.grid(row=nrow,column=0,columnspan=6,sticky="W")
    # event loop for GUI
    split_gui.wait_window(split_gui)

def ncsplit_infilename_browse(split_gui):
    root = Tkinter.Tk(); root.withdraw()
    filename = tkFileDialog.askopenfilename(parent=root,title="Choose an input netCDF file",
                                            initialdir="../Sites")
    root.destroy()
    split_gui.inpathname = ntpath.split(filename)[0]+"/"
    split_gui.infilename.set(ntpath.split(filename)[1])

def ncsplit_outfilename_browse(split_gui):
    root = Tkinter.Tk(); root.withdraw()
    filename = tkFileDialog.asksaveasfilename(parent=root,initialdir=split_gui.inpathname,
                                          title="Choose an output netCDF file")
    root.destroy()
    split_gui.outpathname = ntpath.split(filename)[0]+"/"
    split_gui.outfilename.set(ntpath.split(filename)[1])

def ncsplit_done(split_gui):
    split_gui.destroy()

def ncsplit_run(split_gui):
    msg = " Splitting "+split_gui.infilename.get()
    ncsplit_progress(split_gui,msg)
    infilename = split_gui.inpathname+split_gui.infilename.get()
    outfilename = split_gui.inpathname+split_gui.outfilename.get()
    msg = " Splitting "+infilename
    log.info(msg)
    msg = " Output to "+outfilename
    log.info(msg)
    startdate = str(split_gui.startEntry.get())
    enddate = str(split_gui.endEntry.get())
    # read the input file into the input data structure
    ds_in = nc_read_series(infilename)
    ts = int(ds_in.globalattributes["time_step"])
    ldt_in = ds_in.series["DateTime"]["Data"]
    ldt_in_flag = ds_in.series["DateTime"]["Flag"]
    # create the output data structure
    ds_out = DataStructure()
    # copy the global attributes
    for item in ds_in.globalattributes.keys():
        ds_out.globalattributes[item] = ds_in.globalattributes[item]
    # get the indices of the start and end datetimes
    si = qcutils.GetDateIndex(ldt_in,startdate,ts=ts,default=0,match="exact")
    ei = qcutils.GetDateIndex(ldt_in,enddate,ts=ts,default=len(ldt_in),match="exact")
    # get a list of the series in ds_in
    series_list = [item for item in ds_in.series.keys() if "_QCFlag" not in item]
    # remove the Python datetime series
    for item in ["DateTime","DateTime_UTC"]:
        if item in series_list: series_list.remove(item)
    # loop over the series
    for item in series_list:
        data,flag,attr = qcutils.GetSeriesasMA(ds_in,item,si=si,ei=ei)
        qcutils.CreateSeries(ds_out,item,data,Flag=flag,Attr=attr)
    # deal with the Python datetime series
    ldt_out = ldt_in[si:ei+1]
    ldt_out_flag = ldt_in_flag[si:ei+1]
    ds_out.series["DateTime"] = {}
    ds_out.series["DateTime"]["Data"] = ldt_out
    ds_out.series["DateTime"]["Flag"] = ldt_out_flag
    ds_out.series["DateTime"]["Attr"]= ds_in.series["DateTime"]["Attr"]
    # update the number of records global attribute
    ds_out.globalattributes["nc_nrecs"] = len(ldt_out)
    # update the start and end datetime global attribute
    ds_out.globalattributes["start_date"] = str(ldt_out[0])
    ds_out.globalattributes["end_date"] = str(ldt_out[-1])
    # write the output data structure to a netCDF file
    ncFile = nc_open_write(outfilename)
    nc_write_series(ncFile, ds_out)
    msg = " Finished splitting "+split_gui.infilename.get()
    ncsplit_progress(split_gui,msg)
    log.info(msg)

def ncsplit_progress(split_gui,text):
    """ Update progress message in nc split GUI."""
    split_gui.progress.destroy()
    split_gui.progress = Tkinter.Label(split_gui, text=text)
    split_gui.progress.grid(row=9,column=0,columnspan=6,sticky="W")
    split_gui.update()

def nc_read_series(ncFullName,checktimestep=True,fixtimestepmethod=""):
    """
    Purpose:
     Reads a netCDF file and returns the meta-data and data in a DataStructure.
     The returned data structure is an instance of qcio.DataStructure().
     The data structure consists of:
      1) ds.globalattributes
         A dictionary containing the global attributes of the netCDF file.
      2) ds.series
         A dictionary containing the variable data, meta-data and QC flag
         Each variable dictionary in ds.series contains;
         a) ds.series[variable]["Data"]
            A 1D numpy float64 array containing the variable data, missing
            data value is -9999.
         b) ds.series[variable]["Flag"]
            A 1D numpy int32 array containing the QC flag data for this variable.
         c) ds.series[variable]["Attr"]
            A dictionary containing the variable attributes.
    Usage:
     nc_name = qcio.get_filename_dialog(path="../Sites/Whroo/Data/Processed/")
     ds = qcio.nc_read_series(nc_name)
     where nc_name is the full name of the netCDF file to be read
           ds is the returned data structure
    Side effects:
     This routine checks the time step of the data read from the netCDF file
     against the value of the global attribute "time_step", see qcutils.CheckTimeStep.
     If a problem is found with the time step (duplicate records, non-integral
     time steps or gaps) then qcutils.FixTimeStep is called to repair the time step.
     Fixing non-integral timne steps requires some user input.  The options are to
     quit ([Q]), interpolate ([I], not implemented yet) or round ([R]).  Quitting
     causes the script to exit and return to the command prompt.  Interpolation
     is not implemented yet but will interpolate the data from the original time
     step to a regular time step.  Rounding will round any non-itegral time steps
     to the nearest time step.
    Author: PRI
    Date: Back in the day
    """
    log.info(" Reading netCDF file "+ntpath.split(ncFullName)[1])
    netCDF4.default_encoding = 'latin-1'
    ds = DataStructure()
    # check to see if the requested file exists, return empty ds if it doesn't
    if not qcutils.file_exists(ncFullName,mode="quiet"):
        log.error(' netCDF file '+ncFullName+' not found')
        raise Exception("nc_read_series: file not found")
    # file probably exists, so let's read it
    ncFile = netCDF4.Dataset(ncFullName,'r')
    # now deal with the global attributes
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
        if "time_step" in ds.globalattributes: c.ts = ds.globalattributes["time_step"]
    # get a list of the variables in the netCDF file (not their QC flags)
    varlist = [x for x in ncFile.variables.keys() if "_QCFlag" not in x]
    for ThisOne in varlist:
        # skip variables that do not have time as a dimension
        dimlist = [x.lower() for x in ncFile.variables[ThisOne].dimensions]
        if "time" not in dimlist: continue
        # create the series in the data structure
        ds.series[unicode(ThisOne)] = {}
        # get the data and the QC flag
        data,flag,attr = nc_read_var(ncFile,ThisOne)
        ds.series[ThisOne]["Data"] = data
        ds.series[ThisOne]["Flag"] = flag
        ds.series[ThisOne]["Attr"] = attr
    ncFile.close()
    # make sure all values of -9999 have non-zero QC flag
    qcutils.CheckQCFlags(ds)
    # get a series of Python datetime objects
    if "time" in ds.series.keys():
        time,f,a = qcutils.GetSeries(ds,"time")
        qcutils.get_datetimefromnctime(ds,time,a["units"])
    else:
        qcutils.get_datetimefromymdhms(ds)
    # round the Python datetime to the nearest second
    qcutils.round_datetime(ds,mode="nearest_second")
    # check the time step and fix it required
    if checktimestep:
        if qcutils.CheckTimeStep(ds):
            qcutils.FixTimeStep(ds,fixtimestepmethod=fixtimestepmethod)
            # update the Excel datetime from the Python datetime
            qcutils.get_xldatefromdatetime(ds)
            # update the Year, Month, Day etc from the Python datetime
            qcutils.get_ymdhmsfromdatetime(ds)
    # tell the user when the data starts and ends
    ldt = ds.series["DateTime"]["Data"]
    msg = " Got data from "+ldt[0].strftime("%Y-%m-%d %H:%M:%S")+" to "+ldt[-1].strftime("%Y-%m-%d %H:%M:%S")
    log.info(msg)
    return ds

def nc_read_todf(ncFullName,var_data=[]):
    """
    Purpose:
     Read an OzFlux netCDF file and return the data in an Pandas data frame.
    Usage:
     df = qcio.nc_read_todf(ncFullName)
      where ncFullName is the full name of the netCDF file.
    Side effects:
     Returns a Pandas data frame
    Author: PRI using code originally written by Ian McHugh
    Date: August 2014
    """
    log.info(" Reading netCDF file "+ncFullName+" to Pandas data frame")
    netCDF4.default_encoding = 'latin-1'
    # check to see if the requested file exists, return empty ds if it doesn't
    if not qcutils.file_exists(ncFullName,mode="quiet"):
        log.error(' netCDF file '+ncFullName+' not found')
        raise Exception("nc_read_todf: file not found")
    # file probably exists, so let's read it
    ncFile = netCDF4.Dataset(ncFullName,"r")
    # now deal with the global attributes
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        gattr = {}
        for attr in gattrlist:
            gattr[attr] = getattr(ncFile,attr)
            if "time_step" in gattr.keys(): c.ts = gattr["time_step"]
    # get a list of Python datetimes from the xlDatetime
    # this may be better replaced with one of the standard OzFluxQC routines
    dates_list=[datetime.datetime(*xlrd.xldate_as_tuple(elem,0)) for elem in ncFile.variables['xlDateTime']]
    # get a list of variables to read from the netCDF file
    if len(var_data)==0:
        # get the variable list from the netCDF file contents
        var_data = ncFile.variables.keys()
    else:
        # add the QC flags to the list entered as an argument
        var_flag = []
        for var in var_data: var_flag.append(var+"_QCFlag")
        var_list = var_data+var_flag
    # read the variables and attributes from the netCDF file
    # create dictionaries to hold the data and the variable attributes
    d = {}
    vattr = {}
    for item in var_list:
        d[item] = ncFile.variables[item][:]
        vattrlist = ncFile.variables[item].ncattrs()
        if len(vattrlist)!=0:
            vattr[item] = {}
            for attr in vattrlist:
                vattr[item][attr] = getattr(ncFile.variables[item],attr)
    ncFile.close()
    # convert the dictionary to a Pandas data frame
    df = pd.DataFrame(d,index=dates_list)
    return df,gattr,vattr

def df_droprecords(df,qc_list=[0,10]):
    # replace configured error values with NaNs
    df.replace(c.missing_value,np.nan)
    # replace unacceptable QC flags with NaNs
    var_list = df.columns.values.tolist()
    data_list = [item for item in var_list if "QCFlag" not in item]
    flag_list = [item for item in var_list if "QCFlag" in item]
    if len(data_list)!=len(flag_list): raise Exception("df_droprecords: number of data and flag series differ")
    eval_string='|'.join(['(df[flag_list[i]]=='+str(i)+')' for i in qc_list])
    #for i in xrange(len(data_list)):
    for i in range(len(data_list)):
        df[data_list[i]]=np.where(eval(eval_string),df[data_list[i]],np.nan)
    # drop the all records with NaNs
    df=df[data_list]
    # return the data frame
    return df

def nc_read_var(ncFile,ThisOne):
    """ Reads a variable from a netCDF file and returns the data, the QC flag and the variable
        attribute dictionary.
    """
    # check the number of dimensions
    nDims = len(ncFile.variables[ThisOne].shape)
    if nDims not in [1,3]:
        msg = "nc_read_var: unrecognised number of dimensions ("+str(nDims)
        msg = msg+") for netCDF variable "+ ThisOne
        raise Exception(msg)
    if nDims==1:
        # single dimension
        data = ncFile.variables[ThisOne][:]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data,dummy = qcutils.MAtoSeries(data)
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    elif nDims==3:
        # 3 dimensions
        data = ncFile.variables[ThisOne][:,0,0]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data,dummy = qcutils.MAtoSeries(data)
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:,0,0]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    # force float32 to float64
    if data.dtype=="float32": data = data.astype(numpy.float64)
    # check for Year, Month etc as int64, force to int32 if required
    if ThisOne in ["Year","Month","Day","Hour","Minute","Second"]:
        if data.dtype=="int64": data = data.astype(numpy.int32)
    # get the variable attributes
    vattrlist = ncFile.variables[ThisOne].ncattrs()
    attr = {}
    if len(vattrlist)!=0:
        for vattr in vattrlist:
            attr[vattr] = getattr(ncFile.variables[ThisOne],vattr)
    return data,flag,attr

def nc_open_write(ncFullName,nctype='NETCDF4'):
    """
    Purpose:
     Opens a netCDF file object for writing.  The opened netCDF file
     object is then passed as an argument to qcio.nc_write_series, where
     the actual writing of data occurs.
    Usage:
     nc_name = '../Sites/Whroo/Data/Processed/all/Whroo_L4.nc'
     nc_file = qcio.nc_open_write(nc_name)
     where nc_name is the ful file name of the netCDF to be written
           nc_file is the returned netCDF file object
    Author: PRI
    Date: Back in the day
    """
    log.info(' Opening netCDF file '+ncFullName+' for writing')
    try:
        ncFile = netCDF4.Dataset(ncFullName,'w',format=nctype)
    except:
        log.error(' Unable to open netCDF file '+ncFullName+' for writing')
        ncFile = ''
    return ncFile

def nc_write_series(ncFile,ds,outputlist=None,ndims=3):
    """
    Purpose:
     Write the contents of a data structure to a netCDF file.
    Usage:
     nc_file = qcio.nc_open_write(nc_name)
     qcio.nc_write_series(nc_file,ds)
     where nc_file is a netCDF file object returned by qcio.nc_open_write
           ds is a data structure
    Author: PRI
    Date: Back in the day
    """
    ldt = ds.series["DateTime"]["Data"]
    ds.globalattributes['QC_version'] = str(cfg.version_name)+' '+str(cfg.version_number)
    for ThisOne in ds.globalattributes.keys():
        if "int" in str(type(ds.globalattributes[ThisOne])):
            ds.globalattributes[ThisOne] = numpy.int32(ds.globalattributes[ThisOne])
        setattr(ncFile,ThisOne,ds.globalattributes[ThisOne])
    t = time.localtime()
    rundatetime = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    setattr(ncFile,'nc_rundatetime',rundatetime)
    # we specify the size of the Time dimension because netCDF4 is slow to write files
    # when the Time dimension is unlimited
    nRecs = int(ds.globalattributes['nc_nrecs'])
    ncFile.createDimension("time",nRecs)
    if ndims==3:
        ncFile.createDimension("latitude",1)
        ncFile.createDimension("longitude",1)
        dims = ("time","latitude","longitude")
    else:
        dims = ("time",)
    if outputlist is None:
        outputlist = ds.series.keys()
    else:
        for ThisOne in outputlist:
            if ThisOne not in ds.series.keys():
                log.warning(" Requested series "+ThisOne+" not found in data structure")
                outputlist.remove(ThisOne)
        if len(outputlist)==0: outputlist = ds.series.keys()
    # can't write an array of Python datetime objects to a netCDF file
    # actually, this could be written as characters
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write the time variable
    nc_time = netCDF4.date2num(ldt,"days since 1800-01-01 00:00:00.0",calendar="gregorian")
    ncVar = ncFile.createVariable("time","d",("time",))
    ncVar[:] = nc_time
    setattr(ncVar,"long_name","time")
    setattr(ncVar,"standard_name","time")
    setattr(ncVar,"units","days since 1800-01-01 00:00:00.0")
    setattr(ncVar,"calendar","gregorian")
    if "time" in outputlist: outputlist.remove("time")
    # now write the latitude and longitude variables
    if ndims==3:
        if "latitude" not in outputlist:
            ncVar = ncFile.createVariable("latitude","d",("latitude",))
            ncVar[:] = qcutils.convert_anglestring(str(ds.globalattributes["latitude"]))
            setattr(ncVar,'long_name','latitude')
            setattr(ncVar,'standard_name','latitude')
            setattr(ncVar,'units','degrees north')
        if "longitude" not in outputlist:
            ncVar = ncFile.createVariable("longitude","d",("longitude",))
            ncVar[:] = qcutils.convert_anglestring(str(ds.globalattributes["longitude"]))
            setattr(ncVar,'long_name','longitude')
            setattr(ncVar,'standard_name','longitude')
            setattr(ncVar,'units','degrees east')
    # now make sure the date and time series are in outputlist
    datetimelist = ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh','Ddd']
    # and write them to the netCDF file
    for ThisOne in sorted(datetimelist):
        if ThisOne in ds.series.keys(): nc_write_var(ncFile,ds,ThisOne,dims)
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write everything else to the netCDF file
    for ThisOne in sorted(outputlist):
        nc_write_var(ncFile,ds,ThisOne,dims)
    # write the coordinate reference system (crs) variable
    if "crs" not in outputlist:
        ncVar = ncFile.createVariable("crs","i",())
        setattr(ncVar,"grid_mapping_name","latitude_longitude")
        setattr(ncVar,"long_name","WGS 1984 datum")
        setattr(ncVar,"longitude_of_prime_meridian","0.0")
        setattr(ncVar,"semi_major_axis","6378137.0")
        setattr(ncVar,"inverse_flattening","298.257223563")
    ncFile.close()

def nc_write_var(ncFile,ds,ThisOne,dim):
    """
    Purpose:
     Function to write data from a series in the data structure to a netCDF variable.
    Usage:
     nc_write_var(ncFile,ds,ThisOne,("time","latitude","longitude"))
      where ncFile is a netCDF file object
            ds is the data structure
            ThisOne is the label of a series in ds
            ("time","latitude","longitude") is the dimension tuple
    Author: PRI
    Date: August 2014
    """
    # get the data type of the series in ds
    dt = get_ncdtype(ds.series[ThisOne]['Data'])
    # force data type to float64 or int32
    if dt not in ["d","i"]:
        dt = "d"
        if ThisOne in ["Year","Month","Day","Hour","Minute","Second"]: dt = "i"
    # create the netCDF variable
    try:
        ncVar = ncFile.createVariable(ThisOne,dt,dim)
    except RuntimeError:
        print ThisOne
        raise Exception("Error writing variable to netCDF file")
    # different writes to the variable depending on whether it is 1D or 3D
    #print ds.globalattributes["nc_nrecs"],ThisOne,len(ds.series[ThisOne]["Data"])
    if len(dim)==1: ncVar[:] = ds.series[ThisOne]['Data'].tolist()
    if len(dim)==3: ncVar[:,0,0] = ds.series[ThisOne]['Data'].tolist()
    # write the attributes
    for attr in ds.series[ThisOne]['Attr']:
        if attr!="_FillValue":
            setattr(ncVar,attr,ds.series[ThisOne]['Attr'][attr])
    # make sure the missing_value attribute is written
    if "missing_value" not in ds.series[ThisOne]['Attr']: setattr(ncVar,"missing_value",c.missing_value)
    # get the data type of the QC flag
    dt = get_ncdtype(ds.series[ThisOne]['Flag'])
    # create the variable
    ncVar = ncFile.createVariable(ThisOne+'_QCFlag',dt,dim)
    # write 1D or 3D
    if len(dim)==1: ncVar[:] = ds.series[ThisOne]['Flag'].tolist()
    if len(dim)==3: ncVar[:,0,0] = ds.series[ThisOne]['Flag'].tolist()
    # set the attributes
    setattr(ncVar,'long_name',ThisOne+'QC flag')
    setattr(ncVar,'units','none')

def xl_open_write(xl_name):
    log.info(' Opening Excel file '+xl_name+' for writing')
    try:
        xl_file = xlwt.Workbook()
    except:
        log.error(' Unable to open Excel file '+xl_name+' for writing')
        xl_file = ''
    return xl_file

def xl_read_flags(cf,ds,level,VariablesInFile):
    # First data row in Excel worksheets.
    FirstDataRow = int(qcutils.get_keyvaluefromcf(cf,["Files",level],"first_data_row")) - 1
    HeaderRow = int(qcutils.get_keyvaluefromcf(cf,['Files','in'],'header_row')) - 1
    # Get the full name of the Excel file from the control file.
    xlFullName = get_filename_from_cf(cf,level)
    # Get the Excel workbook object.
    if os.path.isfile(xlFullName):
        xlBook = xlrd.open_workbook(xlFullName)
    else:
        log.error(' Excel file '+xlFullName+' not found, choose another')
        xlFullName = get_filename_dialog(path='.',title='Choose an Excel file')
        if len(xlFullName)==0:
            return
        xlBook = xlrd.open_workbook(xlFullName)
    ds.globalattributes['xlFullName'] = xlFullName

    for ThisOne in VariablesInFile:
        if 'xl' in cf['Variables'][ThisOne].keys():
            log.info(' Getting flags for '+ThisOne+' from spreadsheet')
            ActiveSheet = xlBook.sheet_by_name('Flag')
            LastDataRow = int(ActiveSheet.nrows)
            HeaderList = [x.lower() for x in ActiveSheet.row_values(HeaderRow)]
            if cf['Variables'][ThisOne]['xl']['name'] in HeaderList:
                xlCol = HeaderRow.index(cf['Variables'][ThisOne]['xl']['name'])
                Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                ds.series[ThisOne]['Flag'] = numpy.array([c.missing_value]*len(Values),numpy.int32)
                for i in range(len(Values)):
                    if Types[i]==2: #xlType=3 means a date/time value, xlType=2 means a number
                        ds.series[ThisOne]['Flag'][i] = numpy.int32(Values[i])
                    else:
                        log.error('  xl_read_flags: flags for '+ThisOne+' not found in xl file')
    return ds

def xl_read_series(cf):
    # Instance the data structure object.
    ds = DataStructure()
    # get the filename
    FileName = get_infilenamefromcf(cf)
    if len(FileName)==0:
        log.error(' in_filename not found in control file')
        return ds
    if not os.path.exists(FileName):
        log.error(' Input file '+FileName+' specified in control file not found')
        return ds
    # convert from Excel row number to xlrd row number
    FirstDataRow = int(qcutils.get_keyvaluefromcf(cf,["Files"],"in_firstdatarow")) - 1
    HeaderRow = int(qcutils.get_keyvaluefromcf(cf,["Files"],"in_headerrow")) - 1
    # get the Excel workbook object.
    log.info(" Opening and reading Excel file "+FileName)
    xlBook = xlrd.open_workbook(FileName)
    log.info(" Opened and read Excel file "+FileName)
    ds.globalattributes['featureType'] = 'timeseries'
    ds.globalattributes['xl_filename'] = FileName
    ds.globalattributes['xl_datemode'] = str(xlBook.datemode)
    xlsheet_names = [x.lower() for x in xlBook.sheet_names()]
    # Get the Excel file modification date and time, these will be
    # written to the netCDF file to uniquely identify the version
    # of the Excel file used to create this netCDF file.
    s = os.stat(FileName)
    t = time.localtime(s.st_mtime)
    ds.globalattributes['xl_moddatetime'] = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    # Loop over the variables defined in the 'Variables' section of the
    # configuration file.
    for ThisOne in cf['Variables'].keys():
        if "Function" in cf['Variables'][ThisOne].keys(): continue
        if 'xl' in cf['Variables'][ThisOne].keys():
            if 'sheet' in cf['Variables'][ThisOne]['xl'].keys():
                xlsheet_name = cf['Variables'][ThisOne]['xl']['sheet']
                if xlsheet_name.lower() in xlsheet_names:
                    log.info(' Getting data for '+ThisOne+' from spreadsheet')
                    xlsheet_index = xlsheet_names.index(xlsheet_name.lower())
                    ActiveSheet = xlBook.sheet_by_index(xlsheet_index)
                    HeaderList = [x.lower() for x in ActiveSheet.row_values(HeaderRow)]
                    if cf['Variables'][ThisOne]['xl']['name'].lower() in HeaderList:
                        LastDataRow = int(ActiveSheet.nrows)
                        ds.series[unicode(ThisOne)] = {}
                        xlCol = HeaderList.index(cf['Variables'][ThisOne]['xl']['name'].lower())
                        Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                        Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                        ds.series[ThisOne]['Data'] = numpy.ones(len(Values),dtype=numpy.float64)*float(c.missing_value)
                        ds.series[ThisOne]['Flag'] = numpy.ones(len(Values),dtype=numpy.int32)
                        # we could use "where" and get rid of this for loop
                        for i in range(len(Values)):
                            if (Types[i]==3) or (Types[i]==2): #xlType=3 means a date/time value, xlType=2 means a number
                                ds.series[ThisOne]['Data'][i] = numpy.float64(Values[i])
                                ds.series[ThisOne]['Flag'][i] = numpy.int32(0)
                    else:
                        log.error('  xl_read_series: series '+ThisOne+' not found in xl file')
                else:
                    log.error('  xl_read_series: sheet '+xlsheet_name+' not found in xl file')
            else:
                log.error('  xl_read_series: key "sheet" not found in control file entry for '+ThisOne)
        else:
            log.error('  xl_read_series: key "xl" not found in control file entry for '+ThisOne)
    ds.globalattributes['nc_nrecs'] = str(len(ds.series['xlDateTime']['Data']))
    return ds

def xl_write_AlternateStats(ds):
    if "alternate" not in dir(ds): return
    # open an Excel file for the fit statistics
    cfname = ds.globalattributes["controlfile_name"]
    cf = get_controlfilecontents(cfname,mode="quiet")
    out_filename = get_outfilenamefromcf(cf)
    # get the Excel file name
    xl_filename = out_filename.replace('.nc','_AlternateStats.xls')
    log.info(' Writing alternate fit statistics to Excel file '+xl_filename)
    # open the Excel file
    xlfile = xlwt.Workbook()
    # list of outputs to write to the Excel file
    date_list = ["startdate","enddate"]
    # loop over the series that have been gap filled using alternate data
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for label in ds.alternate.keys():
        # get the list of values to output with the start and end dates removed
        output_list = ds.alternate[label]["results"].keys()
        for item in date_list:
            if item in output_list: output_list.remove(item)
        # add a sheet with the series label
        xlResultsSheet = xlfile.add_sheet(label)
        xlRow = 9
        xlCol = 0
        for dt in date_list:
            xlResultsSheet.write(xlRow,xlCol,dt)
            for item in ds.alternate[label]["results"][dt]:
                xlRow = xlRow + 1
                xlResultsSheet.write(xlRow,xlCol,item,d_xf)
            xlRow = 9
            xlCol = xlCol + 1
        for output in output_list:
            xlResultsSheet.write(xlRow,xlCol,output)
            for item in ds.alternate[label]["results"][output]:
                xlRow = xlRow + 1
                # xlwt under Anaconda seems to only allow float64!
                xlResultsSheet.write(xlRow,xlCol,numpy.float64(item))
            xlRow = 9
            xlCol = xlCol + 1
    xlfile.save(xl_filename)

def xl_write_SOLOStats(ds):
    if "solo" not in dir(ds): return
    # open an Excel file for the fit statistics
    cfname = ds.globalattributes["controlfile_name"]
    cf = get_controlfilecontents(cfname)
    out_filename = get_outfilenamefromcf(cf)
    # get the Excel file name
    xl_filename = out_filename.replace('.nc','_SOLOStats.xls')
    log.info(' Writing SOLO fit statistics to Excel file '+xl_filename)
    # open the Excel file
    xlfile = xlwt.Workbook()
    # list of outputs to write to the Excel file
    date_list = ["startdate","enddate"]
    output_list = ["n","r_max","bias","rmse","var_obs","var_mod","m_ols","b_ols"]
    # loop over the series that have been gap filled using ACCESS data
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for label in ds.solo.keys():
        # get the list of values to output with the start and end dates removed
        output_list = ds.solo[label]["results"].keys()
        for item in date_list:
            if item in output_list: output_list.remove(item)
        # add a sheet with the series label
        xlResultsSheet = xlfile.add_sheet(label)
        xlRow = 10
        xlCol = 0
        for dt in date_list:
            xlResultsSheet.write(xlRow,xlCol,dt)
            for item in ds.solo[label]["results"][dt]:
                xlRow = xlRow + 1
                xlResultsSheet.write(xlRow,xlCol,item,d_xf)
            xlRow = 10
            xlCol = xlCol + 1
        for output in output_list:
            xlResultsSheet.write(xlRow,xlCol,output)
            for item in ds.solo[label]["results"][output]:
                xlRow = xlRow + 1
                # xlwt under Anaconda seems to only allow float64!
                xlResultsSheet.write(xlRow,xlCol,numpy.float64(item))
            xlRow = 10
            xlCol = xlCol + 1
    xlfile.save(xl_filename)

def xl_write_data(xl_sheet,data):
    """
    Purpose:
     Writes a dictionary to a worksheet in an Excel workbook.
     This routine has 2 arguments,an Excel worksheet instance and
     a dictionary of data to be written out.  The dictionary
     format needs to be:
      1) data["DateTime"]["data"]   - a list of Python datetimes, these will
                                      be written tp the first column of the
                                      worksheet
         data["DateTime"]["units"]  - units of the date time eg "Days", "Years"
         data["DateTime"]["format"] - a format string for xlwt.easyxf eg "dd/mm/yyy"
      2) data[variable]["data"]   - a numpy array of data values
         data[variable]["units"]  - units of the data
         data[variable]["format"] - an xlwt.easyxf format string eg "0.00" for 2 decimal places
         There can be multiple variables but each must follow the above template.
    Usage:
     qcio.xl_write_data(xl_sheet,data)
      where xl_sheet is an Excel worksheet instance
            data     is a dictionary as defined above
    Side effects:
     Writes to an Excel worksheet instance
    Called by:
    Calls:
    Author: PRI
    Date: June 2015
    """
    xlCol = 0
    # write the data to the xl file
    series_list = data.keys()
    xl_sheet.write(1,xlCol,data["DateTime"]["units"])
    nrows = len(data["DateTime"]["data"])
    ncols = len(series_list)
    d_xf = xlwt.easyxf(num_format_str=data["DateTime"]["format"])
    for j in range(nrows):
        xl_sheet.write(j+2,xlCol,data["DateTime"]["data"][j],d_xf)
    series_list.remove("DateTime")
    series_list.sort()
    for item in series_list:
        xlCol = xlCol + 1
        xl_sheet.write(0,xlCol,data[item]["units"])
        xl_sheet.write(1,xlCol,item)
        d_xf = xlwt.easyxf(num_format_str=data[item]["format"])
        for j in range(nrows):
            xl_sheet.write(j+2,xlCol,float(data[item]["data"][j]),d_xf)

def xl_write_series(ds, xlfullname, outputlist=None):
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = int(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlfullname)
    xlfile = xlwt.Workbook(encoding="latin-1")
    # set the datemode
    if "xl_datemode" not in ds.globalattributes:
        if platform.system()=="darwin":
            ds.globalattributes["xl_datemode"] = 0
        else:
            ds.globalattributes["xl_datemode"] = 1
    xlfile.dates_1904 = int(ds.globalattributes["xl_datemode"])
    # add sheets to the Excel file
    xlAttrSheet = xlfile.add_sheet('Attr')
    xlDataSheet = xlfile.add_sheet('Data')
    xlFlagSheet = xlfile.add_sheet('Flag')
    # write the global attributes
    log.info(' Writing the global attributes to Excel file '+xlfullname)
    xlcol = 0
    xlrow = 0
    xlAttrSheet.write(xlrow,xlcol,'Global attributes')
    xlrow = xlrow + 1
    globalattrlist = ds.globalattributes.keys()
    globalattrlist.sort()
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' not in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    # write the variable attributes
    log.info(' Writing the variable attributes to Excel file '+xlfullname)
    xlrow = xlrow + 1
    xlAttrSheet.write(xlrow,xlcol,'Variable attributes')
    xlrow = xlrow + 1
    xlcol_varname = 0
    xlcol_attrname = 1
    xlcol_attrvalue = 2
    variablelist = ds.series.keys()
    if outputlist is None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.warning(" Requested series "+ThisOne+" not found in data structure")
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    outputlist.sort()
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    for ThisOne in outputlist:
        xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
        attributelist = ds.series[ThisOne]['Attr'].keys()
        attributelist.sort()
        for Attr in attributelist:
            xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
            xlAttrSheet.write(xlrow,xlcol_attrvalue,str(ds.series[ThisOne]['Attr'][Attr]))
            xlrow = xlrow + 1
    # write the Excel date/time to the data and the QC flags as the first column
    xlDateTime,f,a = qcutils.GetSeries(ds,"xlDateTime")
    log.info(' Writing the datetime to Excel file '+xlfullname)
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    xlDataSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write(j+3,xlcol,xlDateTime[j],d_xf)
        xlFlagSheet.write(j+3,xlcol,xlDateTime[j],d_xf)
    # remove xlDateTime from the list of variables to be written to the Excel file
    if "xlDateTime" in outputlist: outputlist.remove("xlDateTime")
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
        # put up a progress message
        log.info(' Writing '+ThisOne+' into column '+str(xlcol)+' of the Excel file')
        # write the units and the variable name to the header rows in the xl file
        attrlist = ds.series[ThisOne]['Attr'].keys()
        if 'long_name' in attrlist:
            longname = ds.series[ThisOne]['Attr']['long_name']
        elif 'Description' in attrlist:
            longname = ds.series[ThisOne]['Attr']['Description']
        else:
            longname = None
        if 'units' in attrlist:
            units = ds.series[ThisOne]['Attr']['units']
        elif 'Units' in attrlist:
            units = ds.series[ThisOne]['Attr']['Units']
        else:
            units = None
        xlDataSheet.write(0,xlcol,longname)
        xlDataSheet.write(1,xlcol,units)
        xlDataSheet.write(2,xlcol,ThisOne)
        # loop over the values in the variable series (array writes don't seem to work)
        for j in range(nRecs):
            xlDataSheet.write(j+3,xlcol,float(ds.series[ThisOne]['Data'][j]))
        # check to see if this variable has a quality control flag
        if 'Flag' in ds.series[ThisOne].keys():
            # write the QC flag name to the xk file
            xlFlagSheet.write(2,xlcol,ThisOne)
            # specify the format of the QC flag (integer)
            d_xf = xlwt.easyxf(num_format_str='0')
            # loop over QV flag values and write to xl file
            for j in range(nRecs):
                xlFlagSheet.write(j+3,xlcol,int(ds.series[ThisOne]['Flag'][j]),d_xf)
        # increment the column pointer
        xlcol = xlcol + 1
    print xlfullname
    xlfile.save(xlfullname)

def xlsx_write_series(ds, xlsxfullname, outputlist=None):
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = int(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlsxfullname)
    if "xl_datemode" not in ds.globalattributes:
        if platform.system()=="darwin":
            ds.globalattributes["xl_datemode"] = 0
        else:
            ds.globalattributes["xl_datemode"] = 1
    if int(ds.globalattributes["xl_datemode"])==1:
        xlfile = xlsxwriter.Workbook(xlsxfullname, {'date_1904': True})
    else:
        xlfile = xlsxwriter.Workbook(xlsxfullname, {'date_1904': False})
    # add sheets to the Excel file
    xlAttrSheet = xlfile.add_worksheet('Attr')
    xlDataSheet = xlfile.add_worksheet('Data')
    xlFlagSheet = xlfile.add_worksheet('Flag')
    # write the global attributes
    log.info(' Writing the global attributes to Excel file '+xlsxfullname)
    xlcol = 0
    xlrow = 0
    xlAttrSheet.write(xlrow,xlcol,'Global attributes')
    xlrow = xlrow + 1
    globalattrlist = ds.globalattributes.keys()
    globalattrlist.sort()
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' not in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    # write the variable attributes
    log.info(' Writing the variable attributes to Excel file '+xlsxfullname)
    xlrow = xlrow + 1
    xlAttrSheet.write(xlrow,xlcol,'Variable attributes')
    xlrow = xlrow + 1
    xlcol_varname = 0
    xlcol_attrname = 1
    xlcol_attrvalue = 2
    variablelist = ds.series.keys()
    if outputlist is None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.warning(" Requested series "+ThisOne+" not found in data structure")
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    outputlist.sort()
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    for ThisOne in outputlist:
        xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
        attributelist = ds.series[ThisOne]['Attr'].keys()
        attributelist.sort()
        for Attr in attributelist:
            xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
            xlAttrSheet.write(xlrow,xlcol_attrvalue,str(ds.series[ThisOne]['Attr'][Attr]))
            xlrow = xlrow + 1
    # write the Excel date/time to the data and the QC flags as the first column
    ldt = ds.series["DateTime"]["Data"]
    log.info(' Writing the datetime to Excel file '+xlsxfullname)
    dt_format = xlfile.add_format({'num_format': 'dd/mm/yyyy hh:mm'})
    xlDataSheet.write(2,xlcol,'xlDateTime')
    xlFlagSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write_datetime(j+3,xlcol,ldt[j],dt_format)
        xlFlagSheet.write_datetime(j+3,xlcol,ldt[j],dt_format)
    # remove xlDateTime from the list of variables to be written to the Excel file
    if "xlDateTime" in outputlist: outputlist.remove("xlDateTime")
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
        # put up a progress message
        log.info(' Writing '+ThisOne+' into column '+str(xlcol)+' of the Excel file')
        # write the units and the variable name to the header rows in the xl file
        attrlist = ds.series[ThisOne]['Attr'].keys()
        if 'long_name' in attrlist:
            longname = ds.series[ThisOne]['Attr']['long_name']
        elif 'Description' in attrlist:
            longname = ds.series[ThisOne]['Attr']['Description']
        else:
            longname = None
        if 'units' in attrlist:
            units = ds.series[ThisOne]['Attr']['units']
        elif 'Units' in attrlist:
            units = ds.series[ThisOne]['Attr']['Units']
        else:
            units = None
        xlDataSheet.write(0,xlcol,longname)
        xlDataSheet.write(1,xlcol,units)
        xlDataSheet.write(2,xlcol,ThisOne)
        # loop over the values in the variable series (array writes don't seem to work)
        for j in range(nRecs):
            xlDataSheet.write(j+3,xlcol,float(ds.series[ThisOne]['Data'][j]))
        # check to see if this variable has a quality control flag
        if 'Flag' in ds.series[ThisOne].keys():
            # write the QC flag name to the Excel file
            xlFlagSheet.write(2,xlcol,ThisOne)
            # specify the format of the QC flag (integer)
            flag_format = xlfile.add_format({'num_format': '0'})
            # loop over QC flag values and write to xl file
            for j in range(nRecs):
                xlFlagSheet.write(j+3,xlcol,int(ds.series[ThisOne]['Flag'][j]),flag_format)
        # increment the column pointer
        xlcol = xlcol + 1
    
    xlfile.close()
