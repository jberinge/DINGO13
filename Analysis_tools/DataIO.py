# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 09:29:43 2014

@author: imchugh
"""

import Tkinter, tkFileDialog, tkMessageBox
from configobj import ConfigObj
import netCDF4
import numpy as np
import pandas as pd
import datetime as dt
import ast
import csv
import pdb

def ask_question_dialog(header, question):
    
    root = Tkinter.Tk(); root.withdraw()
    do = tkMessageBox.askyesno(header, question)
    root.destroy()    
    return do

def file_select_dialog():
    """ Open a file select dialog to get path for file retrieval"""
    
    root = Tkinter.Tk(); root.withdraw()
    file_in = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    return file_in

def get_nc_global_attr(file_in):
    
    nc_obj=netCDF4.Dataset(file_in)
    return nc_obj.__dict__
    

def OzFluxQCnc_to_pandasDF(file_in, var_list = None):
    
    nc_obj=netCDF4.Dataset(file_in)
    
    dates_list = netCDF4.num2date(nc_obj.variables['time'][:], 
                                  'days since 1800-01-01 00:00:00')
    
    d_data = {}
    d_attr = nc_obj.__dict__

    if var_list == None: 
        var_list = nc_obj.variables.keys()

    for var in var_list:
        ndims = len(nc_obj.variables[var].shape)
        if ndims == 3:
            d_data[var] = nc_obj.variables[var][:, 0, 0]
        elif ndims == 1:
            d_data[var] = nc_obj.variables[var][:]
    nc_obj.close()

    df = pd.DataFrame(d_data, index = dates_list)
    
    df.replace(-9999, np.nan, inplace = True)    

    return df, d_attr

#------------------------------------------------------------------------------
def OzFluxQCnc_add_variable(append_file, data_array, var_attr_dict):
    
    nc_obj = netCDF4.Dataset(append_file, 'w')
    return nc_obj

#------------------------------------------------------------------------------
def OzFluxQCnc_to_data_structure(file_in, 
                                 var_list = None,
                                 fill_missing_with_nan = True,
                                 output_structure = None,
                                 QC_accept_codes = [],
                                 return_global_attr = False,
                                 return_QC = False):

    """
    Pass the following args: 1) valid file path to .nc file \n
    Optional kwargs: 1) 'var_list' - list of variable name strings
                     2) 'fill_missing_with_nan' - fill any record containing
                         the missing data placeholder [-9999 for OzFluxQC] 
                         with nan
                     3) 'output_structure' - if None returns dict by default, 
                         if 'pandas' then returns pandas DataFrame
                     4) 'QC_accept_codes' - if list contains numbers, finds the
                         QC_flag corresponding to the desired variables and 
                         screens out all records where the flag doesn't equal 
                         the accept code
                     5) 'return_global_attr' - return global attributes as a dict
                     6) 'return_QC' - if variable list is unspecified, returns
                        QC variables as well as primary variables (if var_list
                        is specified, this option is inoperative)  
    
    Returns: 1) .nc file data as a dictionary (or pandas DataFrame - see above) 
                of numpy arrays
             2) global attributes of the .nc file (if specified)
    
    """

    def get_var_from_nc(nc_obj, var):
        
        ndims = len(nc_obj.variables[var].shape)
        if ndims == 3:
            return nc_obj.variables[var][:, 0, 0]
        elif ndims == 1:
            return nc_obj.variables[var][:]
        
    nc_obj=netCDF4.Dataset(file_in)
    dates_array = netCDF4.num2date(nc_obj.variables['time'][:], 
                                   'days since 1800-01-01 00:00:00')

    data_dict = {'date_time': dates_array}
    attr_dict = nc_obj.__dict__

    data_var_list = [var for var in nc_obj.variables.keys() if not 'QCFlag' in var]
    QC_var_list = [var for var in nc_obj.variables.keys() if 'QCFlag' in var]

    if var_list == None:
        if return_QC:
            var_list = data_var_list + QC_var_list
        else:
            var_list = data_var_list
    else:
        if not isinstance(var_list, list): var_list = [var_list]

    # Iterate through vars
    for var in var_list:

        if not var in nc_obj.variables.keys():
            print 'Variable ' + var + ' was not found in specified target file'
            continue

        # Check dimensions and get data from variable
        arr = get_var_from_nc(nc_obj, var)
        
        # Check if returned masked array; if so, fill and replace with nan 
        # (if specified by user)
        if isinstance(arr, np.ma.core.MaskedArray):
            arr = arr.filled()
        try:
            if arr.dtype == 'float64':
                if fill_missing_with_nan:
                    arr[arr == -9999] = np.nan
        except AttributeError:
            continue
        
        # If user has specified flag-dependent data screening, then do (unless
        # var is a QC var)
        if QC_accept_codes:
            if not 'QC' in var:
                QC_var = var + '_QCFlag'
                if QC_var in QC_var_list:
                    QC_arr = get_var_from_nc(nc_obj, QC_var)
                    new_arr = np.empty(len(arr))
                    new_arr[:] = np.nan
                    for code in QC_accept_codes:
                        new_arr[QC_arr == code] = arr[QC_arr == code]
                    arr = new_arr
                else:
                    print 'Missing QC variable for variable ' + var
        
        # Add to dict
        data_dict[var] = arr
        
    nc_obj.close()

    if output_structure == 'pandas':
        data_structure = pd.DataFrame(data_dict, index = dates_array)
        data_structure.drop('date_time', axis = 1, inplace = True)
    else:
        data_structure = data_dict

    if return_global_attr:
        return data_structure, attr_dict
    else:
        return data_structure

#------------------------------------------------------------------------------
def DINGO_df_to_data_structure(DataDF = None, file_in= None, 
                               var_list = None,
                               fill_missing_with_nan = True,
                               output_structure = None,
                               return_global_attr = False):

    if not DataDF.empty:
        df = DataDF
    else:
        df = pd.read_pickle(file_in)  
    
    time_step =  pd.infer_freq(df.index)
    attr_dict = {'time_step': int(time_step[: len(time_step) - 1])}

    all_var_list = df.columns
    
    if var_list == None: 
        var_list = all_var_list
    else:
        if not isinstance(var_list, list): var_list = [var_list]
    
    if output_structure == 'pandas':
        data_structure = df[var_list]
    else:
        data_dict = {'date_time': np.array([pd.Timestamp(rec).to_datetime() 
                                            for rec in df.index])}
        for var in var_list:
            data_dict[var] = np.array(df[var])
        data_structure = data_dict

    if return_global_attr:
        return data_structure, attr_dict
    else:
        return data_structure

#------------------------------------------------------------------------------
def config_to_dict(file_in):
    
    def this_funct(section, key):
        try:
            section[key] = ast.literal_eval(section[key])
        except:
            try:
                section[key] = float(section[key])
            except:
                pass

    cf = ConfigObj(file_in)
    
    cf.walk(this_funct)

    return cf  
    
def dict_to_csv(data_dict, outfile, keyorder = False):
    """
    Writes a csv file from a dictionary (or dictionary of dictionaries);
    Pass the following arguments: 1) dict or dict of dicts
                                  2) the order in which the keys should be 
                                     printed as headers (list of strings)
                                  3) the ouput file path and name (string)
    Returns None - writes to file
    """    
    if not keyorder: 
        keyorder = data_dict.keys()
        keyorder.sort()
    with open(outfile, 'wb') as f:
       writer = csv.DictWriter(f, delimiter = ',', fieldnames = keyorder)
       writer.writeheader()
       if np.all(np.array([isinstance(data_dict[key], dict) for key in data_dict.keys()])):
           for d in data_dict.keys():
               writer.writerow(data_dict[d])
       else:
           writer.writerow(data_dict)
    
    return

def array_dict_to_csv(data_dict, outfile, keyorder = False):
    
    keys = data_dict.keys() if not keyorder else keyorder
    try:
        arr_len_list = list(set([len(data_dict[key]) for key in keys]))
        if len(arr_len_list) > 1:
            print 'Dictionary must contain arrays of equal length!'
            return
    except KeyError:
        print 'Could not find key "' + key + '" in passed dictionary'
        return

    arr_len = arr_len_list[0]
    arr = np.empty([arr_len, len(keys)])
    arr[:] = np.nan
    arr = arr.astype(object)
    for i, key in enumerate(keys):
        arr[:, i] = data_dict[key]
       
    with open(outfile, 'wb') as f:
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(keys)
        for i in range(arr_len):
            writer.writerow(arr[i, :])
            
    return
        
def text_dict_to_text_file(write_dict, outfile):        
        
    with open(outfile, 'w') as f:
        for key in write_dict.keys():
            this_key = str(key)
            this_val = str(write_dict[key])
            output = this_key + ': ' + this_val + '\n'
            f.write(output)        
    
def xlsx_to_pandas(file_in,header=True,header_row=0,skiprows_after_header=0,date_col=True,regularise=True,worksheets=[]):

    xl_book=xlrd.open_workbook(file_in)
    
    d={}
    
    start_date='1900-01-01'
    end_date='2100-01-01'    
    
    if not worksheets:
        get_sheets=xl_book.sheet_names()
    else:
        get_sheets=worksheets
    
    for sheet_name in get_sheets:

        sheet=xl_book.sheet_by_name(sheet_name)
       
        rows=sheet.nrows
        cols=sheet.ncols
        
        if rows==0: print 'Could not find any valid rows'
        if cols==0: print 'Could not find any valid columns'

        if rows!=0 and cols!=0:
            if header==True:
                if date_col==True:
                    column_names=[str(sheet.cell_value(header_row,i)) for i in range(cols)]
                    index=[]
                    for i in xrange(header_row+skiprows_after_header+1,sheet.nrows):
                        try:
                            index.append(dt.datetime(*xlrd.xldate_as_tuple(sheet.cell_value(i,0), xl_book.datemode)))
                        except ValueError:
                            index.append('')
                            print 'Error in sheet '+sheet_name+' at row '+str(i)+'; missing or invalid datetime stamp! Skipping...'
                    df=pd.DataFrame(columns=column_names[1:],index=index)
                    for i in range(1,cols):
                        arr=np.array(sheet.col_values(i)[header_row+skiprows_after_header+1:])
                        arr[arr=='']='-9999'
                        df[column_names[i]]=arr.astype(np.float)
                    if regularise==True:
                        df_freq=pd.infer_freq(df.index)
                        df_ind=pd.date_range(start=df.index[0],end=df.index[-1],freq=df_freq)
                        df=df.reindex(df_ind)
                    d[sheet_name]=df
        else:
            d[sheet_name]=pd.DataFrame()
    return d          
    
    # Multiple sheets are returned as dictionary (pandas dataframe) objects
    # Note XLRD cell type codes:
    #    XL_CELL_EMPTY: 0
    #    XL_CELL_TEXT: 1 (STRING)
    #    XL_CELL_NUMBER: 2 (FLOAT)
    #    XL_CELL_DATE: 3 (FLOAT)
    #    XL_CELL_BOOLEAN: 4 (INT)
    #    XL_CELL_ERROR: 5 (INTERNAL EXCEL CODE)
    #    XL_CELL_BLANK: 6 (EMPTY STRING)