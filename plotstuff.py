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

import pandas as pd
from pandas.tools.plotting import scatter_matrix
import datetime as dt
import xlrd
import string
import numpy as np
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['legend.fancybox'] = True
from matplotlib.backends.backend_pdf import PdfPages
from pylab import * 

#from reportlab.lib.enums import TA_JUSTIFY
#from reportlab.lib import colors
#from reportlab.lib.pagesizes import A4
#from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image
#from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
#import reportlab.pdfgen 
#from reportlab.lib.units import inch


#import AWS_Flux_Correlate_plot_V1
from scipy import stats


def plotandPDF2timeseries(title,xdata1,xdata2,xlabel,ylabel,xunits,yunits,Site_ID,VarToCorrelate,mypathforResults):
	
    plt.plot(xdata1.values, 'g-', label='Gap Filled')
    plt.plot(xdata2.values, 'y--', label='Flux Tower original')	
    #Set the scale mins and maxs
    plt.title(title +' Variable  '+VarToCorrelate+ ' at ' +Site_ID)
    plt.xlabel(xlabel + ' ('+xunits+')')
    plt.ylabel(ylabel + ' ('+yunits+')')
    plt.legend(shadow=True, fancybox=True,loc='best')
    
    #Output to PDF using PdfPages a backend for MatPlotLib
    fname_graph=mypathforResults+'/'+title +'-Variable'+VarToCorrelate+ ' at ' +Site_ID+'.pdf'
    # Create the PdfPages object to which we will save the pages:
    pdf = PdfPages(fname_graph)
    savefig(pdf, format='pdf',facecolor='w', edgecolor='w') # note the format='pdf' argument!
    #pdf.close()
    #show()	
    close()		
    pdf.close()
    
	
    
	