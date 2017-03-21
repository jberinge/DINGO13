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
#   The index name is expected to be found on the first row
#
# Programmed by Jason (Dec 1, 2012)
############################################################################
#import httplib2
#from BeautifulSoup import BeautifulSoup, SoupStrainer
import pandas as pd
import csv
import datetime as dt
import numpy as np
from pylab import *
import os
import netCDF4
import scipy
import matplotlib.pyplot as plt

from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from pylab import plot, title, show , legend, xlim, ylim, xlabel,ylabel, title, savefig, close,figtext
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    
    #r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    #The Savitzky-Golay filter removes high frequency noise from data.
    #It has the advantage of preserving the original shape and
    #features of the signal better than other types of filtering
    #approaches, such as moving averages techniques.
    #Parameters
    #----------
    #y : array_like, shape (N,)
        #the values of the time history of the signal.
    #window_size : int
        #the length of the window. Must be an odd integer number.
    #order : int
        #the order of the polynomial used in the filtering.
        #Must be less then `window_size` - 1.
    #deriv: int
        #the order of the derivative to compute (default = 0 means only smoothing)
    #Returns
    #-------
    #ys : ndarray, shape (N)
        #the smoothed signal (or it's n-th derivative).
    #Notes
    #-----
    #The Savitzky-Golay is a type of low-pass filter, particularly
    #suited for smoothing noisy data. The main idea behind this
    #approach is to make for each point a least-square fit with a
    #polynomial of high order over a odd-sized window centered at
    #the point.
    #Examples
    #--------
    #t = np.linspace(-4, 4, 500)
    #y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    #ysg = savitzky_golay(y, window_size=31, order=4)
    #import matplotlib.pyplot as plt
    #plt.plot(t, y, label='Noisy signal')
    #plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    #plt.plot(t, ysg, 'r', label='Filtered signal')
    #plt.legend()
    #plt.show()
    #References
    #----------
    #.. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       #Data by Simplified Least Squares Procedures. Analytical
       #Chemistry, 1964, 36 (8), pp 1627-1639.
    #.. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       #W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       #Cambridge University Press ISBN-13: 9780521880688
    #"""
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def excel_date(date1):
    temp = dt.datetime(1899, 12, 31)
    delta = date1 - temp
    return float(delta.days) + (float(delta.seconds) / 86400)

def MODIS_average_calc(x):
    columns=shape(x)[0]
    number_of_grids=(columns - 3)/2
    if MODIS_productID=='MOD09A1':               #MODIS reflectance
        if QC_tolerance=='tight':
            QC_values=[1073741824]
        else:
            QC_values=[1073741824,1075576832,1075576833,1075576834,1075576835,1610612736,1077411841]
    elif MODIS_productID=='MOD13Q1':        #MODIS EVI, NDVI
        if QC_tolerance=='tight':
            QC_values=[2116]
        else:
            QC_values=[2048, 2049, 2052, 2053, 2112, 2113, 2116, 2117, 2560, 2561, 2564, 2565, 2624, 2625, 2628, 2629]        
    elif MODIS_productID=='MOD15A2':        #MODIS LAI/fPAR
        if QC_tolerance=='tight':
            QC_values=[0]
        else:
            QC_values=[0, 2, 24, 26, 32, 34, 56, 58] 
    elif MODIS_productID=='MOD17A2':        #MODIS GPP and PsNet
        if QC_tolerance=='tight':
            QC_values=[0]
        else:
            QC_values=[0, 2, 24, 26, 32, 34, 56, 58] 
    elif MODIS_productID=='MOD16A2':        #MODIS ET, PET, LE
        if QC_tolerance=='tight':
            QC_values=[0]
        else:
            QC_values=[0, 2, 24, 26, 32, 34, 56, 58] 
    elif MODIS_productID=='MOD11A2':        #MODIS LST day and night
        if QC_tolerance=='tight':
            QC_values=[0]
        else:
            QC_values=[0, 17, 65] 
    else:
        if QC_tolerance=='tight':
            QC_values=[0]
        else:
            QC_values=[0,1]
            
    i=0           #reset counter for averaging
    j=0           #reset cummulative  for averaging
    for gridcell in range(number_of_grids):
        Data_column=(gridcell*2)+3
        QC_column=(gridcell*2)+4  
        if int(x.iloc[QC_column]) in QC_values: 
            j=x.iloc[Data_column]+j
            i=i+1
    if i>0:
        x=j/i
    else:
        x=np.nan
    return x

#################################
# Configuration
#################################


def MODIS_ingest_interp(FLUXDataframe,myBaseforResults,inputMODIS_base,Site_ID,MODIS_key,FluxFreq):
   
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/MODIS"):
        os.mkdir(myBaseforResults+"/MODIS")
    mypathforResults=myBaseforResults+"/MODIS/"
    
    files_to_process = [MODISfile for MODISfile in os.listdir(inputMODIS_base) if (MODIS_key in MODISfile) ]
    
    #the files have an average but we will calulcate our own and use different
    #QC criteria either 'tight' or 'normal'
    global QC_tolerance
    QC_tolerance='normal'    #either 'tight' or 'normal'   
    
    for new_file in files_to_process:   
        
        #Get some inout from the csv file
        filepath=inputMODIS_base+'/'+new_file
        with open(filepath, 'rb') as f:
            mycsv = csv.reader(f)
            mycsv = list(mycsv)
            global MODIS_label1
            MODIS_label1 = mycsv[1][2]
            global MODIS_label
            MODIS_label=MODIS_label1+"_new"
            global MODIS_productID
            MODIS_productID = mycsv[0][7]
            global MODIS_units
            MODIS_units = mycsv[2][2]
        print "Processing MODIS label", MODIS_label1, "MODIS productID", MODIS_productID, "MODIS units", MODIS_units
        # For MODIS data assume the following format forM
        #TOA5	Howard Springs	CR800	12345	CR800.Std.21	CPU:NONE	12345	MOD11A2														
        #TS	RN	Kelvin	AVE_POINTS	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx	Kelvin	0bx
        #		Avg	Sum																		
        #5/03/2000 0:00	1	301.4244444	9	300.82	0	300.88	0	301.42	0	300.94	0	301.54	0	301.56	0	301.88	0	301.9	0	301.88	0
        #13/03/2000 0:00	2	-9999	0	0	2	0	2	0	2	0	2	0	2	0	2	0	2	0	2	0	2
        #21/03/2000 0:00	3	298.6355556	9	299.32	0	299.22	0	298.14	65	299.2	0	299.16	0	298.04	65	299.24	0	297.7	65	297.7	65
        #29/03/2000 0:00	4	299.4355556	9	299.68	65	300.12	65	299.36	65	299.7	65	299.24	65	299.38	65	299.72	65	299	65	298.72	65
        
        
        #read in csv file using Pandas DF
        MODIS_product=pd.io.parsers.read_csv(filepath,index_col=0,prefix=MODIS_label1, skiprows=[0,2,3], na_values=['-9999.0'], parse_dates=True, keep_date_col=False)
        
        #The MODIS csv files have already done the averaging based on QC Filter Conditions 000 and 001. This is column 2
        #Condition 000 represents the highest QC Filter possible and 001 represents a reliable and usable QC Filter, though not to the standard of 000.
        
        # Now the hard bit.  Take data frame and get the MODIS average across columns depending on MODIS QC flag
        MODIS_product[MODIS_label]= MODIS_product.apply(MODIS_average_calc, axis=1)
    
        #Prepare the DF
        #Now make a copy of the index and asign it to float.  Then use that as input to teh spline function
        MODIS_product['MODISdatecopy']=MODIS_product.index
        #MODIS_converted['MODISdatecopy'] = MODIS_converted['MODISdatecopy'].astype('d') 
        MODIS_product['MODISdatecopy_flt']=np.nan
        #Convert the date column which is datetime to float so it can be passed to interp or spline function
        MODIS_product.MODISdatecopy_flt = MODIS_product.MODISdatecopy.apply(excel_date)
        
        #Get the function for interpolation now before expanding the series to 30 minute 
        # interpolate on the non nan drop all other rows
        #Assign index values from dataframe to temporary series
        df2 = MODIS_product.dropna(subset=[MODIS_label,'MODISdatecopy_flt']) 
        
         
        
        xvalues= df2['MODISdatecopy_flt']
        yvalues= df2[MODIS_label]
        
        #f1=scipy.interpolate.interp1d(xvalues, yvalues, bounds_error=False, kind='cubic')
        #f2=scipy.interpolate.BarycentricInterpolator(xvalues, yvalues)
        #f4=scipy.interpolate.interp1d(xvalues, yvalues, bounds_error=False, kind='quadratic')
        f=scipy.interpolate.interp1d(xvalues, yvalues, bounds_error=False, kind='slinear')
        
        #Now INdex the Dataframe using the datestamp MODISdate
        #MODIS_indexed=MODIS_product.set_index(['MODISdate'])
        #Now we have padded out to 30 minutes we need to create a new continuous ORDINAL time column
        # Change to 30 minute frequency and no fill
        MODIS_converted = MODIS_product.asfreq(FluxFreq)
        
        #Now make a copy of the index and asign it to float.  Then use that as input to teh spline function
        MODIS_converted['MODISdatecopy_temp']=MODIS_converted.index
        #Convert the date column which is datetime to float so it can be passed to interp or spline function
        MODIS_converted['MODISdatecopy_new'] = MODIS_converted.MODISdatecopy_temp.apply(excel_date)
        
        MODIS_interp_label=MODIS_label+'_interp'
        MODIS_smooth_label=MODIS_label+'_smooth'
        
        #apply the interpolation to new expanded dataframe
        MODIS_converted[MODIS_interp_label] = MODIS_converted['MODISdatecopy_new'].apply(f)
        MODIS_converted.name = ['MODIS']
        
        #Try the filter savitzky_golay(y, window_size, order, deriv=0, rate=1):
        SGinput=MODIS_converted[MODIS_interp_label]
        
        MODIS_avg_SG=savitzky_golay(SGinput, 10001, 4, 0, 1)
        #new6=savitzky_golay(SGinput, 65, 4, 0, 1)
        #new7=savitzky_golay(SGinput, 65, 5, 0, 1)
        #new8=savitzky_golay(SGinput, 65, 6, 0, 1)
        
        #new10=savitzky_golay(SGinput, 151, 4, 0, 1)
        #new11=savitzky_golay(SGinput, 201, 4, 0, 1)
        #new12=savitzky_golay(SGinput, 251, 4, 0, 1)
        #new13=savitzky_golay(SGinput, 301, 4, 0, 1)
        
        #Do some plots
        #################
        
        
        plot_date(MODIS_converted.index, MODIS_converted[MODIS_label],'yo',label='MODIS 8-day values')
        plot_date(MODIS_converted.index, MODIS_converted[MODIS_interp_label],'g:',label='MODIS spline interpolate')
        plot_date(MODIS_converted.index, MODIS_avg_SG,'r--',label='MODIS SG smoothing')
        title='MODIS_output_interpolated and smoothed_for_'+Site_ID+' '+MODIS_label+ ', '+MODIS_productID+'QC_'+QC_tolerance
        xlabel('Year')
        ylabel(MODIS_label+'  '+MODIS_units) 
        legend()
        savefig(mypathforResults+'MODIS_output_interp and smooth_for_'+Site_ID+' '+MODIS_label+'QC_'+QC_tolerance)
        close()
        
        #Now to get the LIST of data from the SG filter we need to make it into a Dataframe. 
        #Series wont work.  So get a temp date series from exisiting dataframe  
        dateseries=MODIS_converted.index
        to_merge=pd.DataFrame(MODIS_avg_SG, index=dateseries, columns=[MODIS_smooth_label])
        to_merge.name = [MODIS_smooth_label]
        
        #Then CONCAT these together.
        #Use the CONCAT function to join columns (Axis =1)
        NewMODIS=pd.concat([MODIS_converted, to_merge], axis=1)
        
        #Write out to CSV
        NewMODIS.to_csv(mypathforResults+'MODIS_output_interpolated and smoothed_for_'+Site_ID+' '+MODIS_label+'QC_'+QC_tolerance+'.csv', sep=',')
        
        #Remove any existing columns before joing or it doesnt work
        try:
            del FLUXDataframe[MODIS_label]
            del FLUXDataframe[MODIS_interp_label]
            del FLUXDataframe[MODIS_smooth_label]
            print "Deleted existing MODIS columns in dataframe"
        except :
            pass
           
        #Concatenate all the series together
        FLUXDataframe=FLUXDataframe.join(NewMODIS[[MODIS_label,MODIS_interp_label,MODIS_smooth_label]],how='left')
        
        print "DONE processing MODIS "+MODIS_label
    
    print "DONE saving files"
    #Write out final file with all MODIS products AND Fluxdrame to CSV
    FLUXDataframe.to_csv(mypathforResults+'MODIS merged with tower_for_'+Site_ID+'QC_'+QC_tolerance+'.csv', sep=',')
    # Save dataframe
    FLUXDataframe.to_pickle(mypathforResults+'MODIS merged with tower_for_'+Site_ID+'QC_'+QC_tolerance+'.df')
    
    print "DONE processing all MODIS"
    return (FLUXDataframe)