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
import numpy as np
import os
import datetime as dt
import pylab as pl
import meteorologicalfunctions as metfuncs
import time

from pylab import figure, ioff, clf, contourf, ion, draw, show
from ffnet._tests import runtest
from ffnet import ffnet, mlgraph, readdata, tmlgraph, imlgraph
from numpy import array


def Doplots_diurnal(mypathforResults,PlottingDF, Site_ID,versionID):

    print "Doing diurnal plot for month "
    #Do Diurnal Plots for all 12 months
    #create an X axis series for all 24 hours
    t = np.arange(1, 25, 1)

    item_list=['Fre_Con','Fc_ustar','GPP_Con','Fc_Con','Fc_Lasslop']
    Plottemp = PlottingDF[item_list]

	    
    #figure(1)
    pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    pl.subplot(321)
    pl.title("Diurnal partitioning month = 1")
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==1)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==1)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==1)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==1)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux')    
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(322)
    pl.title('Diurnal partitioning month = 3')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==3)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==3)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==3)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==3)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(323)
    pl.title('Diurnal partitioning month = 5')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==5)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==5)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==5)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==5)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    
    pl.subplot(324)
    pl.title('Diurnal partitioning month = 7')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==7)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==7)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==7)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==7)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(325)
    pl.title('Diurnal partitioning month = 9')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==9)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
	try:
	    xdata1a=Plottemp[(PlottingDF.index.month==9)][item_list[0]].groupby([lambda x: x.hour]).mean()
	    plotxdata1a=True
	except:
	    plotxdata1a=False
	try:
	    xdata1b=Plottemp[(PlottingDF.index.month==9)][item_list[1]].groupby([lambda x: x.hour]).mean()
	    plotxdata1b=True
	except:
	    plotxdata1b=False 
	try:
	    xdata1c=Plottemp[(PlottingDF.index.month==9)][item_list[2]].groupby([lambda x: x.hour]).mean()
	    plotxdata1c=True
	except:
	    plotxdata1c=False 	
	try:
	    xdata1d=Plottemp[(PlottingDF.index.month==9)][item_list[3]].groupby([lambda x: x.hour]).mean()
	    plotxdata1d=True
	except:
	    plotxdata1d=False 	
	    
	if plotxdata1a==True:
	    pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
	if plotxdata1b==True:
	    pl.plot(t,xdata1b,'b',label=str(item_list[1]))
	if plotxdata1c==True:
	    pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
	if plotxdata1d==True:
	    pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
	pl.ylabel('Flux') 
	pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(326)
    pl.title('Diurnal partitioning month = 11')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==11)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==11)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==11)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==11)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    figure(1)
    pl.suptitle('Carbon partitioning ensemble diurnal average at '+Site_ID + '_' + versionID)
    pl.subplots_adjust(top=0.85)
    pl.tight_layout()  
    pl.savefig(mypathforResults+'/ANN ensemble diurnal average at '+Site_ID + '_' + versionID)
    #pl.show() 
    pl.close()
    time.sleep(1)

def Doplots_diurnal_Fre(mypathforResults,PlottingDF, Site_ID,versionID):

    print "Doing diurnal Fre plot for month "
    #Do Diurnal Plots for all 12 months
    #create an X axis series for all 24 hours
    t = np.arange(1, 25, 1)

    item_list=['Fre_Con','Fc_ustar','Fre_Lasslop']
    Plottemp = PlottingDF[item_list]

	    
    #figure(1)
    pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    pl.subplot(321)
    pl.title('Diurnal partitioning month = 1')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==1)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==1)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==1)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==1)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux')    
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(322)
    pl.title('Diurnal partitioning month = 3')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==3)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==3)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==3)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==3)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(323)
    pl.title('Diurnal partitioning month = 5')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==5)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==5)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==5)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==5)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    
    pl.subplot(324)
    pl.title('Diurnal partitioning month = 7')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==7)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==7)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==7)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==7)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(325)
    pl.title('Diurnal partitioning month = 9')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==9)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
	try:
	    xdata1a=Plottemp[(PlottingDF.index.month==9)][item_list[0]].groupby([lambda x: x.hour]).mean()
	    plotxdata1a=True
	except:
	    plotxdata1a=False
	try:
	    xdata1b=Plottemp[(PlottingDF.index.month==9)][item_list[1]].groupby([lambda x: x.hour]).mean()
	    plotxdata1b=True
	except:
	    plotxdata1b=False 
	try:
	    xdata1c=Plottemp[(PlottingDF.index.month==9)][item_list[2]].groupby([lambda x: x.hour]).mean()
	    plotxdata1c=True
	except:
	    plotxdata1c=False 	
	try:
	    xdata1d=Plottemp[(PlottingDF.index.month==9)][item_list[3]].groupby([lambda x: x.hour]).mean()
	    plotxdata1d=True
	except:
	    plotxdata1d=False 	
	    
	if plotxdata1a==True:
	    pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
	if plotxdata1b==True:
	    pl.plot(t,xdata1b,'b',label=str(item_list[1]))
	if plotxdata1c==True:
	    pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
	if plotxdata1d==True:
	    pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
	pl.ylabel('Flux') 
	pl.legend(shadow=True, fancybox=True,loc='best')
    
    pl.subplot(326)
    pl.title('Diurnal partitioning month = 11')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==11)][item_list[0]].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==11)][item_list[1]].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    try:
	xdata1c=Plottemp[(PlottingDF.index.month==11)][item_list[2]].groupby([lambda x: x.hour]).mean()
	plotxdata1c=True
    except:
	plotxdata1c=False 	
    try:
	xdata1d=Plottemp[(PlottingDF.index.month==11)][item_list[3]].groupby([lambda x: x.hour]).mean()
	plotxdata1d=True
    except:
	plotxdata1d=False 	
	
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=str(item_list[0])) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=str(item_list[1]))
    if plotxdata1c==True:
	pl.plot(t,xdata1c,'k',label=str(item_list[2]))	
    if plotxdata1d==True:
	pl.plot(t,xdata1d,'g',label=str(item_list[3]))    
    pl.ylabel('Flux') 
    pl.legend(shadow=True, fancybox=True,loc='best')
    
    figure(1)
    pl.suptitle('Fre ensemble diurnal average at '+Site_ID + '_' + versionID)
    pl.subplots_adjust(top=0.85)
    pl.tight_layout()  
    pl.savefig(mypathforResults+'/Fre ensemble diurnal average at '+Site_ID + '_' + versionID)
    #pl.show() 
    pl.close()
    time.sleep(1)

    
def Do_sensitivity_plot1(mypathforResults,New_combined, Site_ID,versionID):
    x_data=[130,100,80,60,40,20,10]
    ydata_NEE=[New_combined['Fc_ustar_130pc'].mean(),New_combined['Fc_ustar_100pc'].mean(),New_combined['Fc_ustar_80pc'].mean(),New_combined['Fc_ustar_60pc'].mean(), New_combined['Fc_ustar_40pc'].mean(),New_combined['Fc_ustar_20pc'].mean(),New_combined['Fc_ustar_10pc'].mean()]
    ydata_Re=[New_combined['Fre_Con_130pc'].mean(),New_combined['Fre_Con_100pc'].mean(),New_combined['Fre_Con_80pc'].mean(),New_combined['Fre_Con_60pc'].mean(), New_combined['Fre_Con_40pc'].mean(),New_combined['Fre_Con_20pc'].mean(),New_combined['Fre_Con_10pc'].mean()]

    #ydata_GPP=[New_combined['GPP_Con_130pc'].mean(),New_combined['GPP_Con_100pc'].mean(),New_combined['GPP_Con_80pc'].mean(),New_combined['GPP_Con_60pc'].mean(), New_combined['GPP_Con_40pc'].mean(),New_combined['GPP_Con_20pc'].mean(),New_combined['GPP_Con_10pc'].mean(),New_combined['GPP_Con_none'].mean()]
    #originally used the GPP variable but this is predicted by ANN so doesnt vary with ustar.  To do this we would need to re run the entire data
    #set back throiugh the annalysis.  So simply here we take GPP as (Nee-Re).
    from operator import sub
    ydata_GPP = map(sub, ydata_NEE, ydata_Re)    
   
    print "Doing ustar sensitivity plot "

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')

    
    ##Convert all units to t C ha-2 yr-1
    convert= 60*60*24*365 *12 / 1000000 * 10000/1000000
    ydata_GPP[:]   =[x * convert for x in ydata_GPP]
    ydata_Re[:]    =[x * convert for x in ydata_Re]
    ydata_NEE[:]   =[x * convert for x in ydata_NEE]

    ax1=pl.subplot(1,1,1)
    
    ax1.plot(x_data,ydata_GPP, 'o', color='#66CC00',markersize=14) 
    ax1.plot(x_data,ydata_Re, 'd', color='#990033',markersize=14) 
    ax1.plot(x_data,ydata_NEE, 's', color='#CC9900',markersize=14) 

    pl.figtext(0.6,0.20,'Average U star used '+ str("{0:.2f}".format(New_combined['ustar_used'].mean())),
               fontsize = 20, bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    #ax1.plot(x_data,ydata_GPP, '-', color='#66CC00') 
    #ax1.plot(x_data,ydata_Re, '-', color='#990033') 
    #ax1.plot(x_data,ydata_NEE, '-', color='#CC9900') 
    
    ax1.legend(loc='upper left')
    pl.xlabel('Ustar threshold sensitivity (% of calculated value)',size=16)
    pl.ylabel('Carbon Flux components annual average  t C ha-2 yr-1',size=16)
    
    pl.suptitle('Plot of Ustar threshold sensitivity for '+Site_ID +'_'+versionID,size=20)  
    
    pl.savefig(mypathforResults+'/Plot of Ustar threshold sensitivity for '+Site_ID +'_'+versionID+'.png') 
    #pl.show()
    pl.close()
   
    print 'Closed Plot of Ustar threshold sensitivity for '+Site_ID +'_'+versionID

def Do_sensitivity_plot2(mypathforResults,New_combined, Site_ID,versionID):
    x_data=[130,100,80,60,40,20,10,0]    
    #ydata_NEE=[New_combined['Fre_Con_130pc'].mean(),New_combined['Fre_Con_100pc'].mean(),New_combined['Fre_Con_80pc'].mean(),New_combined['Fre_Con_60pc'].mean(), New_combined['Fre_Con_40pc'].mean(),New_combined['Fre_Con_20pc'].mean(),New_combined['Fre_Con_10pc'].mean(),New_combined['Fre_Con'].mean()]
    ydata_NEE=[New_combined['Fc_ustar_130pc'].mean(),New_combined['Fc_ustar_100pc'].mean(),New_combined['Fc_ustar_80pc'].mean(),New_combined['Fc_ustar_60pc'].mean(), New_combined['Fc_ustar_40pc'].mean(),New_combined['Fc_ustar_20pc'].mean(),New_combined['Fc_ustar_10pc'].mean(),New_combined['Fc_Con'].mean()]
    
    print "Doing ustar sensitivity plot Fc only"

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')

    ##Convert all units to t C ha-2 yr-1
    convert= 60*60*24*365 *12 / 1000000 * 10000/1000000
    #ydata_GPP[:]   =[x * convert for x in ydata_GPP]
    #ydata_Re[:]    =[x * convert for x in ydata_Re]
    ydata_NEE[:]   =[x * convert for x in ydata_NEE]

    ax1=pl.subplot(1,1,1)
    ax1.plot(x_data,ydata_NEE, 's', color='#CC9900',markersize=14) 

    pl.figtext(0.6,0.20,'Average U star used '+ str("{0:.2f}".format(New_combined['ustar_used'].mean())),
               fontsize = 20, bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    ax1.legend(loc='upper left')
    pl.xlabel('Ustar threshold sensitivity (% of calculated value)',size=16)
    pl.ylabel('Carbon Flux components annual average  t C ha-2 yr-1',size=16)
    
    pl.suptitle('Plot of Ustar threshold sensitivity NEE for '+Site_ID +'_'+versionID,size=20)  
    
    pl.savefig(mypathforResults+'/Plot of Ustar threshold sensitivity NEE for '+Site_ID +'_'+versionID+'.png') 
    #pl.show()
    pl.close()
    
    print 'Closed Plot of Ustar threshold sensitivity for '+Site_ID +'_'+versionID

def GPP_calculate(myBaseforResults,New_combined,Site_ID,versionID, ustarrange):     
    
    ###########################################################################################################
    ##                 START MAIN CODE
    ###########################################################################################################
    
    print "Starting GPP Calculations"
    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/GPP"):
        os.mkdir(myBaseforResults+"/GPP")
    mypathforResults=myBaseforResults+"/GPP"  
    
    #Set columns to missing first
    New_combined['Fre_Con'] = -9999
    New_combined['Fc_ustar'] = -9999
    New_combined['GPP_Con'] = -9999
    
    New_combined['Fre_Lasslop'] = -9999
    New_combined['Fc_Lasslop'] = -9999
    
    
    #New_combined['GPP_Lasslop'] = -9999
    
    #Remember   #Can select period of 'night' (2) and or 'evening' (3) and 'day' (1)
    
    # Use the foloowing threshold for ustar
    # Use the user chosen ustar thershols ' ustar_used'
    # This is chosen by the user in the configuration and can be a manual value, Barr or Reichstein or auto.
    # These are calculated in the FFNET_Fre code and saved as 'ustar_used'
    # Create a new column of constructed GPP based on nighttime and ustar greater than ustar_used
    New_combined['Fre_Con']=New_combined['Fc'][New_combined['ustar']>New_combined['ustar_used']][New_combined['day_night']!=1]

    New_combined[['Fre_Con', 'Fc_ustar', 'GPP_Con', 'Fre_Lasslop', 'Fc_Lasslop','ustar','Fc' ]] = New_combined[['Fre_Con', 'Fc_ustar', 'GPP_Con', 'Fre_Lasslop', 'Fc_Lasslop','ustar','Fc']].astype(np.float64)

    
    #Then fill any remaining with ANN Fre
    New_combined['Fre_Con'][New_combined['Fre_Con'].isnull()]=New_combined['Fre_NN']
    New_combined['Fre_Lasslop']=New_combined['Fre_day']

        
    #Then daytime make Fc_ustar equal to the obs if valid
    New_combined['Fc_ustar']=New_combined['Fc'][New_combined['day_night']==1]
    New_combined['Fc_Lasslop']=New_combined['Fc'][New_combined['day_night']==1]
    
    #Then any remaining values fill with ANN Fc
    New_combined['Fc_ustar'][New_combined['Fc_ustar'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
    New_combined['Fc_Lasslop'][New_combined['Fc_Lasslop'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
    
    #Simply make Fc_ustar = Fre at night
    New_combined['Fc_ustar'][New_combined['day_night']!=1]=New_combined['Fre_Con']
    #LAsslop nighttime use daytime derived Fre using Lasslop
    New_combined['Fc_Lasslop'][New_combined['day_night']!=1]=New_combined['Fre_Lasslop']
    
    #Create a GPP column based on Fc-Re
    New_combined['GPP_Con']=New_combined['Fc_ustar']-New_combined['Fre_Con']
    #New_combined['GPP_Lasslop']=New_combined['Fc_Lasslop']-New_combined['Fre_Lasslop']
	
    #Force to zero during the night
    New_combined['GPP_Con'][New_combined['day_night']!=1]=0
    #New_combined['GPP_Lasslop'][New_combined['day_night']!=1]=0
    

    
    if ustarrange=='Yes':

	#Here output fluxes for each method Reich max, Reich 3 month window and Barr
	############ Ustar threshold to Barr ustar over entire period #######################
	#Set columns to missing first
	New_combined['Fre_Con_Barr'] = -9999
	New_combined['Fc_ustar_Barr'] = -9999
	New_combined['GPP_Con_Barr'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_Barr']=New_combined['Fc'][New_combined['ustar']>New_combined['ustar_Barr']][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_Barr'][New_combined['Fre_Con_Barr'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_Barr']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_Barr'][New_combined['Fc_ustar_Barr'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_Barr'][New_combined['day_night']!=1]=New_combined['Fre_Con_Barr']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_Barr']=New_combined['Fc_ustar_Barr']-New_combined['Fre_Con_Barr']
	#Force to zero during the night
	New_combined['GPP_Con_Barr'][New_combined['day_night']!=1]=0	
	
	#Here we output a whole range of GPP, Fre and Fc for different methods and sensitivities for ustar_used (i.e. 0.25, 0.5, 1.5 etc.)
	############ Ustar threshold to Reich Max ustar over entire period #######################
	#Set columns to missing first
	New_combined['Fre_Con_Reich_max'] = -9999
	New_combined['Fc_ustar_Reich_max'] = -9999
	New_combined['GPP_Con_Reich_max'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_Reich_max']=New_combined['Fc'][New_combined['ustar']>New_combined['ustar_Reich_max']][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_Reich_max'][New_combined['Fre_Con_Reich_max'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_Reich_max']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_Reich_max'][New_combined['Fc_ustar_Reich_max'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_Reich_max'][New_combined['day_night']!=1]=New_combined['Fre_Con_Reich_max']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_Reich_max']=New_combined['Fc_ustar_Reich_max']-New_combined['Fre_Con_Reich_max']
	#Force to zero during the night
	New_combined['GPP_Con_Reich_max'][New_combined['day_night']!=1]=0	

	############ Ustar threshold to 3 month variable ustar #######################
	#Set columns to missing first
	New_combined['Fre_Con_Reich_var'] = -9999
	New_combined['Fc_ustar_Reich_var'] = -9999
	New_combined['GPP_Con_Reich_var'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_Reich_var']=New_combined['Fc'][New_combined['ustar']>New_combined['ustar_Reich_var']][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_Reich_var'][New_combined['Fre_Con_Reich_var'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_Reich_var']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_Reich_var'][New_combined['Fc_ustar_Reich_var'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_Reich_var'][New_combined['day_night']!=1]=New_combined['Fre_Con_Reich_var']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_Reich_var']=New_combined['Fc_ustar_Reich_var']-New_combined['Fre_Con_Reich_var']
	#Force to zero during the night
	New_combined['GPP_Con_Reich_var'][New_combined['day_night']!=1]=0

	#Here we output a whole range of GPP, Fre and Fc for different sensitivities for ustar_used (i.e. all nighttime removed, ustar 130%, 100%, 80%, 60%, 40% and none)
	
	############ All nightime data removed #######################
	#Set columns to missing first
	New_combined['Fre_Con_all'] = -9999
	New_combined['Fc_ustar_all'] = -9999
	New_combined['GPP_Con_all'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	#New_combined['Fre_Con_all']=New_combined['Fc'][New_combined['ustar']>New_combined['ustar_all']][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_all']=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_all']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_all'][New_combined['Fc_ustar_all'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_all'][New_combined['day_night']!=1]=New_combined['Fre_Con_all']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_all']=New_combined['Fc_ustar_all']-New_combined['Fre_Con_all']
	#Force to zero during the night
	New_combined['GPP_Con_all'][New_combined['day_night']!=1]=0

	############ Ustar threshold to 130% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_130pc'] = -9999
	New_combined['Fc_ustar_130pc'] = -9999
	New_combined['GPP_Con_130pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_130pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*1.3)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_130pc'][New_combined['Fre_Con_130pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_130pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_130pc'][New_combined['Fc_ustar_130pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_130pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_130pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_130pc']=New_combined['Fc_ustar_130pc']-New_combined['Fre_Con_130pc']
	#Force to zero during the night
	New_combined['GPP_Con_130pc'][New_combined['day_night']!=1]=0
	
	############ Ustar threshold to 100% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_100pc'] = -9999
	New_combined['Fc_ustar_100pc'] = -9999
	New_combined['GPP_Con_100pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_100pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*1.0)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_100pc'][New_combined['Fre_Con_100pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_100pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_100pc'][New_combined['Fc_ustar_100pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_100pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_100pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_100pc']=New_combined['Fc_ustar_100pc']-New_combined['Fre_Con_100pc']
	#Force to zero during the night
	New_combined['GPP_Con_100pc'][New_combined['day_night']!=1]=0	

	############ Ustar threshold to 80% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_80pc'] = -9999
	New_combined['Fc_ustar_80pc'] = -9999
	New_combined['GPP_Con_80pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_80pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*0.8)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_80pc'][New_combined['Fre_Con_80pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_80pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_80pc'][New_combined['Fc_ustar_80pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_80pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_80pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_80pc']=New_combined['Fc_ustar_80pc']-New_combined['Fre_Con_80pc']
	#Force to zero during the night
	New_combined['GPP_Con_80pc'][New_combined['day_night']!=1]=0
	
	############ Ustar threshold to 60% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_60pc'] = -9999
	New_combined['Fc_ustar_60pc'] = -9999
	New_combined['GPP_Con_60pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_60pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*0.6)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_60pc'][New_combined['Fre_Con_60pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_60pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_60pc'][New_combined['Fc_ustar_60pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_60pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_60pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_60pc']=New_combined['Fc_ustar_60pc']-New_combined['Fre_Con_60pc']
	#Force to zero during the night
	New_combined['GPP_Con_60pc'][New_combined['day_night']!=1]=0	
	
	############ Ustar threshold to 40% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_40pc'] = -9999
	New_combined['Fc_ustar_40pc'] = -9999
	New_combined['GPP_Con_40pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_40pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*0.4)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_40pc'][New_combined['Fre_Con_40pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_40pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_40pc'][New_combined['Fc_ustar_40pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_40pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_40pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_40pc']=New_combined['Fc_ustar_40pc']-New_combined['Fre_Con_40pc']
	#Force to zero during the night
	New_combined['GPP_Con_40pc'][New_combined['day_night']!=1]=0	
	
	
	############ Ustar threshold to 20% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_20pc'] = -9999
	New_combined['Fc_ustar_20pc'] = -9999
	New_combined['GPP_Con_20pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_20pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*0.2)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_20pc'][New_combined['Fre_Con_20pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_20pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_20pc'][New_combined['Fc_ustar_20pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_20pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_20pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_20pc']=New_combined['Fc_ustar_20pc']-New_combined['Fre_Con_20pc']
	#Force to zero during the night
	New_combined['GPP_Con_20pc'][New_combined['day_night']!=1]=0	
	
	############ Ustar threshold to 10% of max #######################
	#Set columns to missing first
	New_combined['Fre_Con_10pc'] = -9999
	New_combined['Fc_ustar_10pc'] = -9999
	New_combined['GPP_Con_10pc'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_10pc']=New_combined['Fc'][New_combined['ustar']>(New_combined['ustar_used']*0.1)][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_10pc'][New_combined['Fre_Con_10pc'].isnull()]=New_combined['Fre_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_10pc']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_10pc'][New_combined['Fc_ustar_10pc'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_10pc'][New_combined['day_night']!=1]=New_combined['Fre_Con_10pc']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_10pc']=New_combined['Fc_ustar_10pc']-New_combined['Fre_Con_10pc']
	#Force to zero during the night
	New_combined['GPP_Con_10pc'][New_combined['day_night']!=1]=0
	
	############ Ustar threshold to none #######################
	#Set columns to missing first
	New_combined['Fre_Con_none'] = -9999
	New_combined['Fc_ustar_none'] = -9999
	New_combined['GPP_Con_none'] = -9999
	#Create a new column of constructed GPP based on nighttime and ustar greater than ustar max
	New_combined['Fre_Con_none']=New_combined['Fc'][New_combined['day_night']!=1]
	#Then fill any remaining with ANN Fre
	New_combined['Fre_Con_none'][New_combined['Fre_Con_none'].isnull()]=New_combined['Fc_NN']
	#Then daytime make Fc_ustar equal to the obs if valid
	New_combined['Fc_ustar_none']=New_combined['Fc'][New_combined['day_night']==1]
	#Then any remaining values fill with ANN Fc
	New_combined['Fc_ustar_none'][New_combined['Fc_ustar_none'].isnull()]=New_combined['Fc_NN'][New_combined['day_night']==1]
	#Simply make Fc_ustar = Fre at night
	New_combined['Fc_ustar_none'][New_combined['day_night']!=1]=New_combined['Fre_Con_none']
	#Create a GPP column based on Fc-Re
	New_combined['GPP_Con_none']=New_combined['Fc_ustar_none']-New_combined['Fre_Con_none']
	#Force to zero during the night
	New_combined['GPP_Con_none'][New_combined['day_night']!=1]=0
	
	#Save seperate file for sensitivity
	output_list=['Fc','Fc_Con','Fc_ustar','Fc_Lasslop','Fre_NN','Fre_Con','Fre_Lasslop','GPP_Con','GPP_Lasslop','ustar','ustar_Barr','ustar_Reich_max','ustar_Reich_var','ustar_used','day_night',
	    'Fre_Con_Reich_max', 'Fc_ustar_Reich_max', 'GPP_Con_Reich_max', 'Fre_Con_Reich_var', 'Fc_ustar_Reich_var', 'GPP_Con_Reich_var', 
	    'Fre_Con_all', 'Fc_ustar_all', 'GPP_Con_all', 'Fre_Con_130pc', 'Fc_ustar_130pc', 'GPP_Con_130pc', 'Fre_Con_100pc', 'Fc_ustar_100pc', 'GPP_Con_100pc', 
	    'Fre_Con_80pc', 'Fc_ustar_80pc', 'GPP_Con_80pc','Fre_Con_60pc', 'Fc_ustar_60pc', 'GPP_Con_60pc',
	    'Fre_Con_40pc', 'Fc_ustar_40pc', 'GPP_Con_40pc','Fre_Con_20pc', 'Fc_ustar_20pc', 'GPP_Con_20pc',
	    'Fre_Con_10pc', 'Fc_ustar_10pc', 'GPP_Con_10pc', 'Fre_Con_none','Fc_ustar_none', 'GPP_Con_none']
	#Write out to CSV 
	New_combined[output_list].to_csv(mypathforResults+'Carbon flux components and ustar method and sensitivity_for_'+Site_ID + '_' + versionID + '.csv', sep=',')
	# Save dataframe.         
	New_combined[output_list].to_pickle(mypathforResults  +'Carbon flux components and ustar method and sensitivity_for_'+Site_ID+ '_' + versionID +'.df')      	
	
	################################################
	# Plots
	###############################################
	#Plot Fc, GPP and Re
	Do_sensitivity_plot1(mypathforResults,New_combined, Site_ID,versionID)
	#plot Fc only
	Do_sensitivity_plot2(mypathforResults,New_combined, Site_ID,versionID)
	
    ################################################
    # Plots
    ################################################
    Doplots_diurnal(mypathforResults,New_combined, Site_ID,versionID)
    Doplots_diurnal_Fre(mypathforResults,New_combined, Site_ID,versionID)

    
    #Remove excess variables from DF before return
    #Here is the list to remove
    listtodrop=['Fre_Con_Reich_max', 'Fc_ustar_Reich_max', 'GPP_Con_Reich_max', 'Fre_Con_Reich_var', 'Fc_ustar_Reich_var', 'GPP_Con_Reich_var', 
	    'Fre_Con_all', 'Fc_ustar_all', 'GPP_Con_all', 'Fre_Con_130pc', 'Fc_ustar_130pc', 'GPP_Con_130pc', 'Fre_Con_100pc', 'Fc_ustar_100pc', 'GPP_Con_100pc', 
	    'Fre_Con_80pc', 'Fc_ustar_80pc', 'GPP_Con_80pc','Fre_Con_60pc', 'Fc_ustar_60pc', 'GPP_Con_60pc',
            'Fre_Con_40pc', 'Fc_ustar_40pc', 'GPP_Con_40pc','Fre_Con_20pc', 'Fc_ustar_20pc', 'GPP_Con_20pc',
            'Fre_Con_10pc', 'Fc_ustar_10pc', 'GPP_Con_10pc', 'Fre_Con_none','Fc_ustar_none', 'GPP_Con_none'] 
    
    
    for item in listtodrop:
	try:
	    New_combined.drop(item,axis=1, inplace=True)
	except:
	    pass
    
    #####################################
    #  File stuff
    ###################################################
    
    
    print"Finished GPP outputs at "+ Site_ID    
    return (New_combined)




