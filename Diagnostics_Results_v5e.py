############################################################################
# This script produces Diagnostics and results
# DINGO produces a variety of diagnostic, summary and results plots that assist the user in immediate visualisation of the data.  These 
# plots enableenables rapid identification and correction of instrument or processing errors.  The suite of outputs includes:
#1.	plots for the identification of energy balance non-closure (Franssen et al., 2010; Leuning et al., 2012; Oken, 2008)  
#showing scatter plots difference between turbulent fluxes and available energy for; a) all hours of 30 minute values, b) daytime 
#hours of 30 minute values, c) nighttime hours of 30 minute values and d) daily means.  An example is shown in Fig. 8. 
#2. plots that show the difference between actual tower observations and best alternate 'guess' time series for each variable
#such as those constructed using BoM AWS station data or the ANN model output.  This essentially compares what we think the variable
#'should' be with what it really is.  If there is a big difference then it could indicate instrument or processing differences that
#may need require further examination. 
#3. calculation and plots of weekly timeseries of ecosystem scale water use efficiency (WUE) following Beer et al. (2009), 
#radiation use efficiency (RUE) following Garbulsky et al. (2010), energy balance closure following (Twine and Kucharik, 2008), 
#evaporative fraction (EF) following (Zhou and Wang, 2016) and the Bowen ratio (BR) (Bowen, 1926).  These plots can be used to 
#identify physically and physiologically inconsistent data periods.
#4.	graphs indicating the missing data for all variables including percentage of each month that is gap-filled and the percentage
#of any data not gap-filled (Fig. 9a) and monthly time series plots where data with more than 30% of data gap-filled is shown in 
#light blue (Fig. 9b).
#5.	summary figures of variables in fingerprint style for non-gap-filled and gap-filled variables using code from OzFluxQC 
#(Isaac et al., 2016) (Fig 10). 
#6.	summary timeseries plots of daily means and a 30 day running mean of net ecosystem exchange (Fc), ecosystem respiration 
#(Fre) and gross primary production (GPP) (Fig. 11).
#7.	results showing the cumulative carbon (GPP, Fc and Fre) and water (Fe and precipitation) by year (Fig. 12).
#
# Last updated by Jason Beringer (5/5/2016)
############################################################################

import pandas as pd
import numpy as np
import os
import datetime as dt
import pylab as pl
import meteorologicalfunctions as metfuncs
from dateutil.relativedelta import relativedelta
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

from pylab import *
from numpy import array
import matplotlib.dates as mdates
import matplotlib.cbook as cbook

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def Doplots_diurnal(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,item,versionID):
    ANN_label=str(item+"_NN")     #Do Monthly Plot
    print "Doing diurnal plot for month "
    #Do Diurnal Plots for all 12 months
    #create an X axis series for all 24 hours
    t = np.arange(1, 25, 1)
    NN_label='Fc'
    Plottemp = PlottingDF[[NN_label,item]][PlottingDF['day_night']!=3]
    #Plottemp = PlottingDF[[NN_label,item]].dropna(how='any')
	    
    figure(1)
    pl.subplot(321)
    pl.title('Diurnal '+item+' month = 1')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==1)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==1)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')    
 
    pl.subplot(322)
    pl.title('Diurnal '+item+' month = 3')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==3)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==3)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')      

    pl.subplot(323)
    pl.title('Diurnal '+item+' month = 5')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==5)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==5)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')      
    
    pl.subplot(324)
    pl.title('Diurnal '+item+' month = 7')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==7)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==7)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.subplot(325)
    pl.title('Diurnal '+item+' month = 9')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==9)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==9)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    pl.subplot(326)
    pl.title('Diurnal '+item+' month = 11')
    try:
	xdata1a=Plottemp[(PlottingDF.index.month==11)][item].groupby([lambda x: x.hour]).mean()
	plotxdata1a=True
    except:
	plotxdata1a=False
    try:
	xdata1b=Plottemp[(PlottingDF.index.month==11)][NN_label].groupby([lambda x: x.hour]).mean()
	plotxdata1b=True
    except:
	plotxdata1b=False 
    if plotxdata1a==True:
	pl.plot(t,xdata1a,'r',label=item) 
    if plotxdata1b==True:
	pl.plot(t,xdata1b,'b',label=NN_label)
    pl.ylabel('Flux')  
    
    
    figure(1)
    pl.suptitle('ANN ensemble diurnal average for variable '+item+' at '+Site_ID+ 'for' +versionID)
    pl.subplots_adjust(top=0.85)
    pl.tight_layout()  
    pl.savefig(mypathforResults+'/ANN ensemble diurnal average for variable '+item+' at '+Site_ID+'_'+versionID)
    #pl.show() 
    pl.close()

def Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID):   
   
    print "Doing Mean plots for "+Site_ID
    #create a string of the items to be plotted to be used later to save filename
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'.pdf')

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(1, number_of_subplots+1, 1)):
	item=list_in[i]
	item_Con=item+'_Con'
	item_Con_flag=item+'_Con_QCFlag'
	
	if item_Con=="Fc_ustar_Con" : item_Con="Fc_ustar"
	if item_Con=="Fc_ustar_Con" : item_Con_flag="Fc_Con_QCFlag"
	
	by = lambda x: lambda y: getattr(y, x)
	xdata1a=New_combined[item_Con].groupby([by('year'),by(plot_freq)]).mean()	
	xdata1b=New_combined[[item_Con,item_Con_flag]].groupby([by('year'),by(plot_freq)]).mean()
	#Here get a data series that is only periods where the percentage gap filled more than 30%
	xdata1b[xdata1b[item_Con_flag]>30.0]=np.nan
	xdata1b=xdata1b[item_Con]
	
	#ydata=np.arange(len(xdata1a))	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
	if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
	if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
	if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])
	
	ax=pl.subplot(number_of_subplots,1,v)
	#ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	ax.plot(ydata,xdata1a, color='#00BFFF',label=item) #DeepSkyBlue 
	ax.plot(ydata,xdata1b, color='#00008B',linewidth=2,label=item)	 #DarkBlue 
	
	ax.legend(loc='upper right')
	pl.ylabel(item+'\n'+ '('+get_units(item,Ws_label)+')')
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	#if v != number_of_subplots:
	#    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean plots for for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'_'+versionID)  
    
    #pl.savefig(pdf, format='pdf')
    pl.savefig(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
    pl.close()    
    
def Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID):   
   
    print "Doing Mean plots for "+Site_ID
    #create a string of the items to be plotted to be used later to save filename
    running_freq=30    
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'.pdf')

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(1, number_of_subplots+1, 1)):
	item=list_in[i]
	item_Con=item+'_Con'
	
	if item_Con=="Fc_ustar_Con" : item_Con="Fc_ustar"
	
	by = lambda x: lambda y: getattr(y, x)
	xdata1a=New_combined[item_Con].groupby([by('year'),by(plot_freq)]).mean()	
	
	#ydata=np.arange(len(xdata1a))	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
	if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
	if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
	if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])
	
	ax=pl.subplot(number_of_subplots,1,v)
	#ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	ax.plot(ydata,xdata1a, color='#90EE90',label=item)	#'LightGreen '
	#Plot running mean
	xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
	ax.plot(ydata,xdata1_smooth, color='#008000',linewidth=2, label=item +' '+str(running_freq)+' day run mean')  #green
	
	#ax.legend(loc='upper right')
	pl.ylabel(item+'\n'+ '('+get_units(item,Ws_label)+')')
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	#if v != number_of_subplots:
	#    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean plots for for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'_'+versionID + 'with 30 day running mean')  
    
    #pl.savefig(pdf, format='pdf')
    pl.savefig(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID + 'with 30 day running mean')
    #pl.show()
    pl.close()        


def Doplots_at_daily_other(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID):   
   
    print "Doing Other plots for "+Site_ID
    running_freq=30
    #Create new data and variables
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	    
    #Calculate Bowen ratio daily averages used
    tempdata_other['BR']=tempdata_other['Fh_Con']/tempdata_other['Fe_Con']
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['EF']=tempdata_other['Fe_Con']/tempdata_other['Fn_Con']    
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['WUE']=tempdata_other['Fc_Con']/tempdata_other['Fe_Con'] 
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['RUE']=tempdata_other['Fc_Con']/tempdata_other['Fsd_Con'] 
    #Calculate Evaporative fraction daily averages used    
    tempdata_other['EBC']=(tempdata_other['Fh_Con']+tempdata_other['Fe_Con'])/(tempdata_other['Fn_Con'] -tempdata_other['Fg_Con'] )
    
    
    #create a string of the items to be plotted to be used later to save filename
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'.pdf')

    fig=pl.figure(1, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(1, number_of_subplots+1, 1)):
	item=list_in[i]
	xdata1a=tempdata_other[item]	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
	if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
	if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
	if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])
	
	ax=pl.subplot(number_of_subplots,1,v)
	#ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
	ax.plot(ydata,xdata1a, color='#9370DB',label=(item + ' at '+ plot_freq))	#'LightGreen '
	#Plot running mean	
	xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
	
	ax.plot(ydata,xdata1_smooth, color='#800080',linewidth=2, label=(item +' '+str(running_freq)+' day running mean'))
	#Calculate limits for plot.  Take the running mean and add some margins
	ylimit_max=max(xdata1_smooth)
	if ylimit_max>0:
	    ylimit_max=ylimit_max*1.2
	else:
	    ylimit_max=ylimit_max*0.8	
	    
	ylimit_min=min(xdata1_smooth)
	if ylimit_min<0:
	    ylimit_min=ylimit_min*1.2
	else:
	    ylimit_min=ylimit_min*0.8
	pl.ylim((ylimit_min,ylimit_max))
	#ax.legend(loc='upper right')
	pl.ylabel(item,size=14)
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	#if v != number_of_subplots:
	#    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean plots for for '+list_string+ ' at '+Site_ID+ ' freq ' + str(plot_freq)+'_'+versionID +' '+str(running_freq)+' day running mean',size=20)  

    pl.savefig(mypathforResults+'/Mean plots for '+list_string+ ' at '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
    pl.close()

def Doplots_at_daily_carbon_g(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID):   
   
    print "Doing Carbon plots in units mg for "+Site_ID
    #Create new data and variables
    running_freq=30
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=pl.figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    xdata1a=tempdata_other['GPP_Con']*12/1000000*60*60*24
    #xdata2a=tempdata_other['Fc_Con']*12/1000000*60*60*24
    xdata2a=tempdata_other['Fc_ustar']*12/1000000*60*60*24
    xdata3a=tempdata_other['Fre_Con']*12/1000000*60*60*24
    
    startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
    enddate1=New_combined.index[len(New_combined)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(New_combined)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
 
    xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
    xdata2_smooth=smooth(xdata2a,running_freq,window='flat')[0:len(ydata)]
    xdata3_smooth=smooth(xdata3a,running_freq,window='flat')[0:len(ydata)]
	
    pl.plot(ydata,xdata1a, color='#90EE90',label=('GPP'))	
    pl.plot(ydata,xdata1_smooth, color='#006400',linewidth=2.5)

    pl.plot(ydata,xdata2a, color='#FFD700',label=('Fc'))	
    pl.plot(ydata,xdata2_smooth, color='#DAA520',linewidth=2.5)
    
    pl.plot(ydata,xdata3a, color='#FF6347',label=('Fre'))	
    pl.plot(ydata,xdata3_smooth, color='#B22222',linewidth=2.5)
 
    #Calculate limits for plot.  Take the running mean and add some margins
    ylimit_max=max(max(xdata1_smooth),max(xdata2_smooth),max(xdata3_smooth))
    if ylimit_max>0:
	ylimit_max=ylimit_max*1.2
    else:
	ylimit_max=ylimit_max*0.8	
	
    ylimit_min=min(min(xdata1_smooth),min(xdata2_smooth),min(xdata3_smooth))
    if ylimit_min<0:
	ylimit_min=ylimit_min*1.2
    else:
	ylimit_min=ylimit_min*0.8
    pl.ylim((ylimit_min,ylimit_max))

    pl.legend(loc='lower right')
    pl.ylabel('CO2 flux (g C m-2 d-1)',size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    pl.suptitle('Timeseries Carbon plot for '+Site_ID + ' freq ' + str(plot_freq)+' with 30 day running mean '+'_'+versionID,size=24)
    pl.savefig(mypathforResults+'/Results carbon plots in g C m-2 d-1_'+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
    pl.close()

def Doplots_at_daily_EB(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID):   
   
    print "Doing Energy Balance plots for "+Site_ID
    #Create new data and variables
    running_freq=30
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=pl.figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    xdata1a=tempdata_other['Fn_Con']
    xdata2a=tempdata_other['Fe_Con']	
    xdata3a=tempdata_other['Fh_Con'] 
    xdata4a=tempdata_other['Fg_Con']
    
    startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
    enddate1=New_combined.index[len(New_combined)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(New_combined)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
 
    xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
    xdata2_smooth=smooth(xdata2a,running_freq,window='flat')[0:len(ydata)]
    xdata3_smooth=smooth(xdata3a,running_freq,window='flat')[0:len(ydata)]
    xdata4_smooth=smooth(xdata4a,running_freq,window='flat')[0:len(ydata)]

    pl.plot(ydata,xdata1a, color='#848484')	
    pl.plot(ydata,xdata1_smooth, color='#000000',linewidth=2.5)

    pl.plot(ydata,xdata2a, color='#81BEF7')	
    pl.plot(ydata,xdata2_smooth, color='#0101DF',linewidth=2.5)
    
    pl.plot(ydata,xdata3a, color='#F5A9A9')	
    pl.plot(ydata,xdata3_smooth, color='#B40404',linewidth=2.5)
    
    pl.plot(ydata,xdata4a, color='#58FA58')	
    pl.plot(ydata,xdata4_smooth, color='#31B404',linewidth=2.5)    
 
    #Calculate limits for plot.  Take the running mean and add some margins
    ylimit_max=max(max(xdata1_smooth),max(xdata2_smooth),max(xdata3_smooth),max(xdata4_smooth))
    if ylimit_max>0:
	ylimit_max=ylimit_max*1.2
    else:
	ylimit_max=ylimit_max*0.8	
	
    ylimit_min=min(min(xdata1_smooth),min(xdata2_smooth),min(xdata3_smooth),min(xdata4_smooth))
    if ylimit_min<0:
	ylimit_min=ylimit_min*1.2
    else:
	ylimit_min=ylimit_min*0.8
    pl.ylim((ylimit_min,ylimit_max))

    pl.legend(loc='upper right')
    pl.ylabel('Energy Balance flux Daily average (W m-2)',size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    pl.suptitle('Timeseries EB plot for '+Site_ID + ' freq ' + str(plot_freq)+' with 30 day running mean '+'_'+versionID,size=24)
    pl.savefig(mypathforResults+'/Results EB plots '+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
    pl.close()
    
def Doplots_at_daily_carbon_umol(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID):   
   
    print "Doing Carbon plots for "+Site_ID
    #Create new data and variables
    running_freq=30
    by = lambda x: lambda y: getattr(y, x)
    tempdata_other=New_combined.groupby([by('year'),by(plot_freq)]).mean()	     

    fig=pl.figure(1, figsize=(17, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    monthsFmt = mdates.DateFormatter('%M')
    
    xdata1a=tempdata_other['GPP_Con']	
    #xdata2a=tempdata_other['Fc_Con']	
    xdata2a=tempdata_other['Fc_ustar']    
    xdata3a=tempdata_other['Fre_Con']  
    
    startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
    enddate1=New_combined.index[len(New_combined)-1]
    endplot=dt.datetime(enddate1.year,enddate1.month,1)
    totalpoints=len(New_combined)
    if plot_freq=='dayofyear': ydata=np.array([startplot + relativedelta(days=i) for i in xrange(len(xdata1a))])
    if plot_freq=='week': ydata=np.array([startplot + relativedelta(weeks=i) for i in xrange(len(xdata1a))])
    if plot_freq=='month': ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    if plot_freq=='year': ydata=np.array([startplot + relativedelta(years=i) for i in xrange(len(xdata1a))])    
 
    xdata1_smooth=smooth(xdata1a,running_freq,window='flat')[0:len(ydata)]
    xdata2_smooth=smooth(xdata2a,running_freq,window='flat')[0:len(ydata)]
    xdata3_smooth=smooth(xdata3a,running_freq,window='flat')[0:len(ydata)]
	
    pl.plot(ydata,xdata1a, color='#90EE90',label=('GPP'))	
    pl.plot(ydata,xdata1_smooth, color='#006400',linewidth=2.5)

    pl.plot(ydata,xdata2a, color='#FFD700',label=('Fc')	)
    pl.plot(ydata,xdata2_smooth, color='#DAA520',linewidth=2.5)
    
    pl.plot(ydata,xdata3a, color='#FF6347',label=('Fre'))	
    pl.plot(ydata,xdata3_smooth, color='#B22222',linewidth=2.5)
 
    #Calculate limits for plot.  Take the running mean and add some margins
    ylimit_max=max(max(xdata1_smooth),max(xdata2_smooth),max(xdata3_smooth))
    if ylimit_max>0:
	ylimit_max=ylimit_max*1.2
    else:
	ylimit_max=ylimit_max*0.8	
	
    ylimit_min=min(min(xdata1_smooth),min(xdata2_smooth),min(xdata3_smooth))
    if ylimit_min<0:
	ylimit_min=ylimit_min*1.2
    else:
	ylimit_min=ylimit_min*0.8
    pl.ylim((ylimit_min,ylimit_max))

    pl.legend(loc='lower right')
    pl.ylabel('CO2 flux (umol CO2 m-2 s-1)',size=14)
    ## format the ticks
    #ax.xaxis.set_major_locator(years)
    #ax.xaxis.set_major_formatter(yearsFmt)
    #ax.xaxis.set_minor_locator(months)
    #ax.xaxis.set_minor_formatter(monthsFmt)
    #pl.grid(True, which='both') 
    ax.axhline(y=0,color='k')
    #datemin = dt.date(ydata.min().year, 1, 1)
    #datemax = dt.date(ydata.max().year+1, 1, 1)
    #ax.set_xlim(datemin, datemax)
    
    #pl.savefig(pdf, format='pdf')
    pl.suptitle('Timeseries Carbon plot for '+Site_ID + ' freq ' + str(plot_freq)+' with 30 day running mean '+'_'+versionID,size=24)
    pl.savefig(mypathforResults+'/Results carbon plots in umol CO2 m-2 s-1_'+Site_ID + ' freq ' + str(plot_freq)+'_'+versionID)
    #pl.show()
    pl.close()


def get_units(x,Ws_label):
    if x == 'Ta':
	units = 'oC'
    elif x == 'ps':
	units = 'hPa'
    elif x == 'Ah':
	units = 'g m-3'
    elif x == Ws_label:
	units = 'm s-1'	
    elif x == 'Fsd':
	units = 'W m-2'
    elif x == 'Fsu':
	units = 'W m-2'
    elif x == 'Precip':
	units = 'mm'
    elif x == 'Fld':
	units = 'W m-2'
    elif x == 'Flu':
	units = 'W m-2'	
    elif x == 'Fc':
	units = 'umol m-2 s-1'
    elif x == 'Fc_ustar':
	units = 'umol m-2 s-1'
    elif x == 'Fe':
	units = 'W m-2'
    elif x == 'Fh':
	units = 'W m-2'
    elif x == 'Fg':
	units = 'W m-2'	
    else:
	units = ''
    return units

def Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,versionID):   

    print "Doing Mean Monthly DIFF (Best Alt Met - Tower) plot for "+Site_ID
    #create a string of the items to be plotted to be used later to save filename
    list_string=''
    for z in list_in:
	#cretae string list
	list_string=list_string+' '+z

    #pdf = PdfPages(mypathforResults+'/Mean Monthly DIFF plot for '+Site_ID+'.pdf')

    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    number_of_subplots=len(list_in)
    years    = mdates.YearLocator()   # every year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
      
    for i,v in enumerate(xrange(1, number_of_subplots+1, 1)):
	item=list_in[i]
	item_Con=item+'_Con'
	if item in ['Fc','Fe','Fh','Fg']:
	    item_Corr=item+'_NN'
	else:
	    item_Corr=item+'_Corr'	
	
	xdata1a=New_combined[item_Corr].groupby([lambda x: x.year,lambda x: x.month]).mean()-New_combined[item].groupby([lambda x: x.year,lambda x: x.month]).mean()
	#ydata=np.arange(len(xdata1a))	
	
	startplot=dt.datetime(New_combined.index[0].year,New_combined.index[0].month,1)
	enddate1=New_combined.index[len(New_combined)-1]
	endplot=dt.datetime(enddate1.year,enddate1.month,1)
	totalpoints=len(New_combined)
	ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])

	ax=pl.subplot(number_of_subplots,1,v)
	ax.bar(ydata,xdata1a, width=20, color='DarkMagenta',label=item)	
	ax.legend(loc='upper right')
	pl.ylabel(item+'\n'+ '('+get_units(item,Ws_label)+')')
	# format the ticks
	ax.xaxis.set_major_locator(years)
	ax.xaxis.set_major_formatter(yearsFmt)
	ax.xaxis.set_minor_locator(months)
	ax.grid(True, which='both') 
	#horiz line
	ax.axhline(y=0,color='k')	
	
	if v != number_of_subplots:
	    pl.setp(ax.get_xticklabels(), visible=False)
	
	datemin = dt.date(ydata.min().year, 1, 1)
	datemax = dt.date(ydata.max().year+1, 1, 1)
	ax.set_xlim(datemin, datemax)
    
    pl.suptitle('Mean Monthly DIFF (Best Alt Met - Tower) plot for for '+list_string+ ' at '+Site_ID+'_'+versionID)  
    
    #pl.savefig(pdf, format='pdf')
    pl.savefig(mypathforResults+'/Mean Monthly DIFF plot for '+list_string+ ' at '+Site_ID+'_'+versionID)
    #pl.show()
    pl.close()
        
def Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,versionID):   
        list_string=''
	for z in list_in:
	    #cretae string list
	    list_string=list_string+' '+z
	    
	print 'Doing Plots of missing data for '+list_string+ ' at '+Site_ID
	
	#pdf = PdfPages(mypathforResults+'/Plots of missing data for '+list_string+ ' at '+Site_ID+'.pdf')
    
	fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	number_of_subplots=len(list_in)
	years    = mdates.YearLocator()   # every year
	months   = mdates.MonthLocator()  # every month
	yearsFmt = mdates.DateFormatter('%Y') 
	  
	for i,v in enumerate(xrange(1, number_of_subplots+1, 1)):
	    item=list_in[i]
	    #v = v+1
	    
	    if item == 'Fc_ustar': item = 'Fc'
	    
	    Pct_nan_DF=pd.read_csv(mypathforResults+'/'+'Nan counts and Pct filled for '+item+' at ' +Site_ID+'_'+versionID+'.csv')
	    xdata1a=Pct_nan_DF['Pct_nan']
	    xdata1b=Pct_nan_DF['Pct_notfilled']
	    totalpoints=len(Pct_nan_DF)
	    startplot=dt.datetime(int(Pct_nan_DF.ix[0,0]),int(Pct_nan_DF.ix[0,1]),1)
	    endplot=dt.datetime(int(Pct_nan_DF.ix[(totalpoints-1),0]),int(Pct_nan_DF.ix[(totalpoints-1),1]),1)
	    ydata=np.array([startplot + relativedelta(months=i) for i in xrange(len(xdata1a))])
    
	    ax=pl.subplot(number_of_subplots,1,v)
	    ax.bar(ydata,xdata1a, width=20, color='orange',label='Percent Nan')	
	    ax.bar(ydata,xdata1b, width=20, color='g',label='Percent not filled')	
	    ax.legend(loc='upper right')
	    pl.ylabel(item+ ' % missing')
	    # format the ticks
	    ax.xaxis.set_major_locator(years)
	    ax.xaxis.set_major_formatter(yearsFmt)
	    ax.xaxis.set_minor_locator(months)
	    ax.grid(True, which='both') 
	    #horiz line
	    ax.axhline(y=0,color='k')	
	    
	    if v != number_of_subplots:
		pl.setp(ax.get_xticklabels(), visible=False)
	    
	    datemin = dt.date(ydata.min().year, 1, 1)
	    datemax = dt.date(ydata.max().year+1, 1, 1)
	    ax.set_xlim(datemin, datemax)
	    #except:
		#pass
	
	pl.suptitle('Plots of missing data for '+list_string+ ' at '+Site_ID+'_'+versionID)  

	#pl.savefig(pdf, format='pdf')
	pl.savefig(mypathforResults+'/Plots of missing data for '+list_string+ ' at '+Site_ID+'_'+versionID)
	#pl.show()
	pl.close()    
	print 'Closed'

def regressionEBclosure(mypathforResults,New_combined,Site_ID,startdate,enddate,versionID):    
    print "Doing Energy Balance closure plots for "+Site_ID

    #First plot use ALL data and years
    ###################################
    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    #Do first plot - All hours
    #============================
    #Get data for all hours - so dont select any
    #Also select only data that is QCFlag 1 for valid observations
    tempdata=New_combined[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][(New_combined['Fh_Con_QCFlag'] == 1) &  (New_combined['Fe_Con_QCFlag'] == 1) & (New_combined['Fn_Con_QCFlag'] == 1) & (New_combined['Fg_Con_QCFlag'] == 1)].dropna(axis=0,how='any')
    xdata=tempdata['Fn_Con'] - tempdata['Fg_Con']
    ydata=tempdata['Fe_Con'] + tempdata['Fh_Con']
    ax1=pl.subplot(2,2,1)
    pl.title=('Energy Balance closure ALL hours for '+Site_ID )
    ax1.plot(xdata,ydata, 'o', color='#00BFFF') #DeepSkyBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.36,0.57,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax1.plot(x, y, linewidth = 2,label='ALL hours')    
    
    ax1.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')

    
    #Do second plot - Daytime hours
    #=================================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][(New_combined['Fh_Con_QCFlag'] == 1) &  (New_combined['Fe_Con_QCFlag'] == 1) & (New_combined['Fn_Con_QCFlag'] == 1) & (New_combined['Fg_Con_QCFlag'] == 1)].dropna(axis=0,how='any')
    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
    ax2=pl.subplot(2,2,2)
    pl.title=('Energy Balance closure DAYTIME hours for '+Site_ID )  
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax2.plot(xdata,ydata,  'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.8,0.57,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax2.plot(x, y, linewidth = 2,label='DAYTIME hours')    
    
    ax2.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')
    
    #Do third plot - Nighttime hours
    #===============================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][(New_combined['Fh_Con_QCFlag'] == 1) &  (New_combined['Fe_Con_QCFlag'] == 1) & (New_combined['Fn_Con_QCFlag'] == 1) & (New_combined['Fg_Con_QCFlag'] == 1)].dropna(axis=0,how='any')
    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
    ax3=pl.subplot(2,2,3)
    pl.title=('Energy Balance closure NIGHTTIME hours for '+Site_ID )    
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax3.plot(xdata,ydata,  'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression independantly rather than using the NN stats.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.36,0.12,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax3.plot(x, y, linewidth = 2,label='NIGHTTIME hours')    
    
    ax3.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')
    
    #Do forth plot - DAILY AVG hours
    #============================
    #Get data for all hours - so dont select any
    tempdata=New_combined[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][(New_combined['Fh_Con_QCFlag'] == 1) &  (New_combined['Fe_Con_QCFlag'] == 1) & (New_combined['Fn_Con_QCFlag'] == 1) & (New_combined['Fg_Con_QCFlag'] == 1)].dropna(axis=0,how='any')
    
    by = lambda x: lambda y: getattr(y, x)
    tempdata=tempdata.groupby([by('year'),by('month'),by('day')]).mean()    


    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
    ax4=pl.subplot(2,2,4)
    pl.title=('Energy Balance closure DAILY average for '+Site_ID )  
    #ax.bar(ydata,xdata1a, width=20, color='g',label=item)	
    ax4.plot(xdata,ydata,  'o', color='#00BFFF') #DeepSkyBlue 
    #ax.plot(ydata,xdata, color='#00008B',linewidth=2,label=item)	 #DarkBlue 

    #Get regression stats
    #In this function calculate the linear regression.
    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)

    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
    pl.figtext(0.8,0.12,graphtext1, bbox=dict())

    x = np.linspace(min(xdata),max(xdata))
    y = slope * x + intercept
    ax4.plot(x, y, linewidth = 2,label='DAILY averages')    
    
    ax4.legend(loc='upper left')
    pl.ylabel('Fh + Fe (W m-2)')
    pl.xlabel('Fn - Fg (W m-2)')    
    pl.suptitle('Energy balance closure at '+Site_ID+ ' '+str(startdate.year)+ ' to '+str(enddate.year)+'_'+versionID,size=20)  
    
    pl.savefig(mypathforResults+'/'+'Energy balance closure at ' +Site_ID + ' all years'+'_'+versionID) 
    #pl.show()
    pl.close()    
    print 'Closed'

    #Now do PLOTS BY YEAR
    #===========================
    tempdata_year=New_combined[['Fe_Con','Fh_Con','Fg_Con','Fn_Con','day_night']][(New_combined['Fh_Con_QCFlag'] == 1) &  (New_combined['Fe_Con_QCFlag'] == 1) & (New_combined['Fn_Con_QCFlag'] == 1) & (New_combined['Fg_Con_QCFlag'] == 1)]
    tempdata_grouped=tempdata_year.groupby([lambda x: x.year])
    
    for name, group in tempdata_grouped:
	try:
	    #Get year for plot lables
	    plotyear=group.index[0].year
	    
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    tempdata=group.dropna(axis=0,how='any')
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax1=pl.subplot(2,2,1)
	    pl.title=('Energy Balance closure ALL hours for '+Site_ID )
	    ax1.plot(xdata,ydata, 'o', color='#00CED1') #DeepSkyBlue 

	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
	                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
	                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.36,0.57,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2,label='ALL hours')    
	    
	    ax1.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	
	    
	    #Do second plot - Daytime hours
	    #=================================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][group['day_night']==1].dropna(axis=0,how='any')
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax2=pl.subplot(2,2,2)
	    pl.title=('Energy Balance closure DAYTIME hours for '+Site_ID )  	
	    ax2.plot(xdata,ydata,  'o', color='#00CED1') #DeepSkyBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
	                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
	                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.8,0.57,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax2.plot(x, y, linewidth = 2,label='DAYTIME hours')    
	    
	    ax2.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	    
	    #Do third plot - Nighttime hours
	    #===============================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][group['day_night']!=1].dropna(axis=0,how='any')
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax3=pl.subplot(2,2,3)
	    pl.title=('Energy Balance closure NIGHTTIME hours for '+Site_ID )    	
	    ax3.plot(xdata,ydata,  'o', color='#00CED1') #DeepSkyBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
	                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
	                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.36,0.12,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax3.plot(x, y, linewidth = 2,label='NIGHTTIME hours')    
	    
	    ax3.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	    
	    #Do forth plot - DAILY AVG hours
	    #============================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']].dropna(axis=0,how='any')
	    
	    by = lambda x: lambda y: getattr(y, x)
	    tempdata=tempdata.groupby([by('month'),by('day')]).mean()    
	
	
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax4=pl.subplot(2,2,4)
	    pl.title=('Energy Balance closure DAILY average for '+Site_ID )  	
	    ax4.plot(xdata,ydata,  'o', color='#00CED1') #DeepSkyBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression 
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
	                   'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
	                   'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.8,0.12,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax4.plot(x, y, linewidth = 2,label='DAILY averages')    
	    
	    ax4.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')    
	    pl.suptitle('Energy balance closure at '+Site_ID+ ' for year '+str(plotyear)+'_'+versionID,size=20)  
	    
	    pl.savefig(mypathforResults+'/'+'Energy balance closure at '+Site_ID+ ' for year '+str(plotyear)+'_'+versionID) 
	    #pl.show()
	    pl.close()    
	except:
	    print "WARNING plot not completed for ebergy balance at " + Site_ID + " for "+str(plotyear)
	    pass

    #Now do PLOTS BY CATEGORICAL VARIABLE
    #This can be set here.  Currently this is just for RDMF site 
    #===========================
    
    if Site_ID=="RDMF":
	#define dates of phases
	startdate=dt.date(2011,9,11)
	phase1_end=dt.date(2012,3,2)
	phase2_end=dt.date(2012,3,6)
	phase3_end=dt.date(2012,8,6)
	phase4_end=dt.date(2012,8,28)
	phase5_end=dt.date(2013,1,22)
	phase6_end=dt.date(2013,5,23)
	phase7_end=dt.date(2013,7,23)
	
	New_combined['RDMF_Phase']=''
	New_combined['RDMF_Phase'][startdate:phase1_end]="Phase1"
	New_combined['RDMF_Phase'][phase1_end:phase2_end]="Phase2"
	New_combined['RDMF_Phase'][phase2_end:phase3_end]="Phase3"
	New_combined['RDMF_Phase'][phase3_end:phase4_end]="Phase4"
	New_combined['RDMF_Phase'][phase4_end:phase5_end]="Phase5"
	New_combined['RDMF_Phase'][phase5_end:phase6_end]="Phase6"
	New_combined['RDMF_Phase'][phase6_end:phase7_end]="Phase7"	
	
	tempdata_grouped=New_combined[['Fe_Con','Fh_Con','Fg_Con','Fn_Con','day_night']][(New_combined['Fh_Con_QCFlag'] == 1) &  (New_combined['Fe_Con_QCFlag'] == 1) & (New_combined['Fn_Con_QCFlag'] == 1) & (New_combined['Fg_Con_QCFlag'] == 1)].groupby(New_combined['RDMF_Phase'])
    
	for name, group in tempdata_grouped:
	    #try:
	    
	    #Get year for plot lables
	    plot_cat=name
	    
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    tempdata=group.dropna(axis=0,how='any')
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax1=pl.subplot(2,2,1)
	    pl.title=('Energy Balance closure ALL hours for '+Site_ID )
	    ax1.plot(xdata,ydata, 'o', color='#00CED1') #DeepSkyBlue 
    
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.36,0.57,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2,label='ALL hours')    
	    
	    ax1.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	
	    
	    #Do second plot - Daytime hours
	    #=================================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][group['day_night']==1].dropna(axis=0,how='any')
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax2=pl.subplot(2,2,2)
	    pl.title=('Energy Balance closure DAYTIME hours for '+Site_ID )  	
	    ax2.plot(xdata,ydata,  'o', color='#00CED1') #DeepSkyBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.8,0.57,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax2.plot(x, y, linewidth = 2,label='DAYTIME hours')    
	    
	    ax2.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	    
	    #Do third plot - Nighttime hours
	    #===============================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']][group['day_night']!=1].dropna(axis=0,how='any')
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax3=pl.subplot(2,2,3)
	    pl.title=('Energy Balance closure NIGHTTIME hours for '+Site_ID )    	
	    ax3.plot(xdata,ydata,  'o', color='#00CED1') #DeepSkyBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.36,0.12,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax3.plot(x, y, linewidth = 2,label='NIGHTTIME hours')    
	    
	    ax3.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')
	    
	    #Do forth plot - DAILY AVG hours
	    #============================
	    #Get data for all hours - so dont select any
	    tempdata=group[['Fe_Con','Fh_Con','Fg_Con','Fn_Con']].dropna(axis=0,how='any')
	    
	    by = lambda x: lambda y: getattr(y, x)
	    tempdata=tempdata.groupby([by('month'),by('day')]).mean()    
	
	
	    xdata=tempdata['Fn_Con']-tempdata['Fg_Con']
	    ydata=tempdata['Fe_Con']+tempdata['Fh_Con']
	    ax4=pl.subplot(2,2,4)
	    pl.title=('Energy Balance closure DAILY average for '+Site_ID )  	
	    ax4.plot(xdata,ydata,  'o', color='#00CED1') #DeepSkyBlue 
	
	    #Get regression stats
	    #In this function calculate the linear regression 
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'n          ' + str("{0:.2f}".format(len(tempdata))) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err))       )  
	    pl.figtext(0.8,0.12,graphtext1, bbox=dict())
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax4.plot(x, y, linewidth = 2,label='DAILY averages')    
	    
	    ax4.legend(loc='upper left')
	    pl.ylabel('Fh + Fe (W m-2)')
	    pl.xlabel('Fn - Fg (W m-2)')    
	    pl.suptitle('Energy balance closure at '+Site_ID+ ' for year '+str(plot_cat)+'_'+versionID,size=20)  
	    
	    pl.savefig(mypathforResults+'/'+'Energy balance closure at '+Site_ID+ ' for year '+str(plot_cat)+'_'+versionID) 
	    #pl.show()
	    pl.close()    
	    #except:
		#print "WARNING plot not completed for ebergy balance at " + Site_ID + " for "+str(plot_cat)
		#pass	    

def regressionLasslop(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in,versionID):    
    
    for plot_freq in freq_list: 
	for var_to_plot in list_in:
	    print "Doing Fre plots for "+Site_ID + " at freq " + plot_freq
		
	    std_variable=var_to_plot+"_Con"
	    Lasslop_variable=var_to_plot+"_Lasslop"
	    if var_to_plot=='Fc': std_variable='Fc_ustar'
	    by = lambda x: lambda y: getattr(y, x)	
	    
	    # Added below a statement to screen out 'nonsense' values predicted by lasslop
	    #First define the range in the observations then screen outside of those bounds.	
	    #Calculate by group and then reapply
	    tempdata1=New_combined[[std_variable,Lasslop_variable ]]
	    tempdata1=tempdata1.dropna(how='any')
	    tempdata=tempdata1.groupby([by('year'),by(plot_freq)]).mean()	
	    #Min and max of the grouped series
	    obs_min = min(tempdata[std_variable])
	    obs_max = max(tempdata[std_variable])
	    #now use Min and Max to exclude orginal data from Lasslop variable and recalculate
	    tempdata2=New_combined[(New_combined[Lasslop_variable]<obs_max) & (New_combined[Lasslop_variable]>obs_min)][[std_variable,Lasslop_variable ]]
	    tempdata2=tempdata2.dropna(how='any')
	    tempdata3=tempdata2.groupby([by('year'),by(plot_freq)]).mean()
	    #print the number of values excluded from this 	
	    print 'Plotting ' + var_to_plot + ' and excluded ' + str((len(tempdata)-len(tempdata3))) + 'values at frequency ' + plot_freq
	    
	   
	    	    
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    xdata=tempdata3[std_variable]
	    ydata=tempdata3[Lasslop_variable]
	    
	    #Convert all units to g C m-2 period-1
	    if plot_freq=='dayofyear': ydata=tempdata3[Lasslop_variable]*60*60*24/1000*12/44
	    if plot_freq=='dayofyear': xdata=tempdata3[std_variable]*60*60*24/1000*12/44
	    if plot_freq=='week': ydata=tempdata3[Lasslop_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='week': xdata= tempdata3[std_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='month': ydata=tempdata3[Lasslop_variable]*60*60*24/1000*12/44  *30
	    if plot_freq=='month': xdata=tempdata3[std_variable]*60*60*24/1000*12/44     *30
	    
	    ax1=pl.subplot(1,1,1)
	    pl.title=('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID )
	    
	    #Define colours
	    if var_to_plot=='Fre': plotcolour='#990033'
	    if var_to_plot=='GPP': plotcolour='#66CC00'
	    if var_to_plot=='Fc': plotcolour='#CC9900'
	    
	    ax1.plot(xdata,ydata, 'o', color=plotcolour) 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	    num_points=len(ydata)
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err)) +'\n' +
	                   'n points   ' + str("{0:.2f}".format(num_points))       )  
	    pl.figtext(0.7,0.20,graphtext1, bbox=dict(),size=16)
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2)    
	    
	    #plot linear regression 1:1
	    x2 = np.linspace(min(xdata),max(xdata))
	    y2 = 1 * x2 + 0
	    ax1.plot(x2, y2, linewidth = 1,label='1:1',color='k',linestyle='--') 
	    
	    print "Min X  " + str(min(xdata)) + "            Max X  " + str(max(xdata))
	    
	    
	    ax1.legend(loc='upper left')
	    pl.xlabel(std_variable + ' (g C m-2 '+plot_freq+'-1)',size=16)
	    pl.ylabel(Lasslop_variable +' (g C m-2 '+plot_freq+'-1)',size=16)
	    
	    pl.suptitle('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID +'_'+versionID,size=20)  
	    
	    pl.savefig(mypathforResults+'/Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID+'_'+versionID+'.png') 
	    #pl.show()
	    pl.close()    
	    print 'Closed Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID+'_'+versionID

def regressionFre(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in,versionID):    
    
    for plot_freq in freq_list: 
	for var_to_plot in list_in:
	    print "Doing Fre plots for "+Site_ID + " at freq " + plot_freq
		
	    std_variable=var_to_plot
	    Lasslop_variable="Fre_Lasslop"
	    
	    by = lambda x: lambda y: getattr(y, x)	     

	    # Added below a statement to screen out 'nonsense' values predicted by lasslop
	    #First define the range in the observations then screen outside of those bounds.	
	    #Calculate by group and then reapply
	    tempdata1=New_combined[[std_variable,Lasslop_variable ]]
	    tempdata1=tempdata1.dropna(how='any')
	    tempdata=tempdata1.groupby([by('year'),by(plot_freq)]).mean()	
	    #Min and max of the grouped series
	    obs_min = min(tempdata[std_variable])
	    obs_max = max(tempdata[std_variable])
	    #now use Min and Max to exclude orginal data from Lasslop variable and recalculate
	    tempdata2=New_combined[(New_combined[Lasslop_variable]<obs_max) & (New_combined[Lasslop_variable]>obs_min)][[std_variable,Lasslop_variable ]]
	    tempdata2=tempdata2.dropna(how='any')
	    tempdata3=tempdata2.groupby([by('year'),by(plot_freq)]).mean()
	    #print the number of values excluded from this 	
	    print 'Plotting ' + var_to_plot + ' and excluded ' + str((len(tempdata)-len(tempdata3))) + 'values at frequency ' + plot_freq

	
	    #First plot use ALL data and years
	    ###################################
	    fig=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	    #Do first plot - All hours
	    #============================
	    #Get data for all hours - so dont select any
	    xdata=tempdata3[std_variable]
	    ydata=tempdata3[Lasslop_variable]
	    
	    #Convert all units to g C m-2 period-1
	    if plot_freq=='dayofyear': ydata=tempdata3[Lasslop_variable]*60*60*24/1000*12/44
	    if plot_freq=='dayofyear': xdata=tempdata3[std_variable]*60*60*24/1000*12/44
	    if plot_freq=='week': ydata=tempdata3[Lasslop_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='week': xdata= tempdata3[std_variable]*60*60*24/1000*12/44   *7
	    if plot_freq=='month': ydata=tempdata3[Lasslop_variable]*60*60*24/1000*12/44  *30
	    if plot_freq=='month': xdata=tempdata3[std_variable]*60*60*24/1000*12/44     *30
	    
	    ax1=pl.subplot(1,1,1)
	    pl.title=('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID )
	    
	    #Define colours
	    if var_to_plot=='Fre_noct': plotcolour='#FF9966'
	    if var_to_plot=='Fre_NN': plotcolour='#FF6666'
	    if var_to_plot=='Fre_Con': plotcolour='#FF3366'
	    
	    ax1.plot(xdata,ydata, 'o', color=plotcolour) 
	
	    #Get regression stats
	    #In this function calculate the linear regression independantly rather than using the NN stats.
	    #Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	    slope, intercept, r_value, p_value, std_err = stats.linregress(xdata,ydata)
	    num_points=len(ydata)
	    graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
		           'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
		           'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
		           'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
		           'slope SE   ' + str("{0:.2f}".format(std_err)) +'\n' +
	                   'n points   ' + str("{0:.2f}".format(num_points))       )  
	    pl.figtext(0.7,0.20,graphtext1, bbox=dict(),size=16)
	
	    x = np.linspace(min(xdata),max(xdata))
	    y = slope * x + intercept
	    ax1.plot(x, y, linewidth = 2)    
	    
	    #plot linear regression 1:1
	    x2 = np.linspace(min(xdata),max(xdata))
	    y2 = 1 * x2 + 0
	    ax1.plot(x2, y2, linewidth = 1,label='1:1',color='k',linestyle='--') 
	    
	    print "Min X  " + str(min(xdata)) + "            Max X  " + str(max(xdata))
	    
	    ax1.legend(loc='upper left')
	    pl.xlabel(std_variable + ' (g C m-2 '+plot_freq+'-1)',size=16)
	    pl.ylabel(Lasslop_variable +' (g C m-2 '+plot_freq+'-1)',size=16)
	    
	    pl.suptitle('Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID +'_'+versionID,size=20)  
	    
	    pl.savefig(mypathforResults+'/Plot2 '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID+'_'+versionID+'.png') 
	    #pl.show()
	    pl.close()    
	    print 'Closed Plot '+std_variable+' vs '+Lasslop_variable+' at '+plot_freq+' for '+Site_ID+'_'+versionID
 

def cummulative_CO2(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill,versionID,cumm_var_to_plot):    
    print "Doing Cummulative CO2 plots for "+Site_ID
 
    #First plot use ALL data and years
    ###################################
    #Create the Cummulative variables
    #Define label for cummulative
    cumm_label=cumm_var_to_plot+'_Cumm'
    if cumm_var_to_plot=='Fc_ustar':
	cm = get_cmap('Dark2')
	var_input= 'Fc_ustar'
    elif cumm_var_to_plot=='GPP':
	cm = get_cmap('winter')
	var_input= 'GPP_Con'
    elif cumm_var_to_plot=='Fre':
	cm = get_cmap('hot')
	var_input= 'Fre_Con'    
    
    New_combined[cumm_label]=New_combined.groupby([lambda x: x.year])[var_input].cumsum()
    New_combined[cumm_label]=New_combined[cumm_label]
    
    #get to units of g.CO2.-m-2. multiply by 44 (MW of CO2) to get g CO2
    New_combined[cumm_label]=New_combined[cumm_label]*60*30*10000/1000000*44.01*12/44/1000000   

    #Group all data by Year
    tempdata_grouped=New_combined.groupby([lambda x: x.year])

    Num_plots=len(tempdata_grouped)
    
    #First plot Fc
   
    fig = figure(1,figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    ax  = fig.add_subplot(111)
    
    index=0
    for name, group in tempdata_grouped:
	plotyear=group.index[0].year
	index=index+1
	xdata=group[cumm_label]
	plotyear=group.index[0].year
	color_rgb = cm(1.*index/Num_plots)
	ax12=pl.plot(xdata.values, '-', linewidth=2,color=color_rgb,label=cumm_var_to_plot+' '+str(plotyear))

    index=0
    #for name, group in tempdata_grouped:
	#plotyear=group.index[0].year
	#index=index+1
	#xdata=group[cumm_label]
	#plotyear=group.index[0].year
	#color_rgb = cm(1.*index/Num_plots)
	##ax12=pl.plot(xdata, '-', linewidth=2,color=color_rgb,label='Fc '+str(plotyear))
	#x_coord_Fc=xdata[-1]
	#y_coord_Fc=len(xdata)-1    
	#ax.annotate('test', xy=(x_coord_Fc, y_coord_Fc),  xycoords='data',
	            #xytext=(-100, -100), textcoords='offset points',
	            #arrowprops=dict(arrowstyle="->")
	            #)	
	
    # Set the ticks and labels...
    #ticks = np.linspace(0, len(xdata), 12)
    #labels = range(ticks.size)
    #pl.xticks(ticks, labels)  
    
    pl.legend(loc='upper left')
    pl.ylabel('Cummulative Carbon flux (t C ha-1 y-1)')
    pl.xlabel('Time (30 min intervals)')
    pl.suptitle('Cummulative CO2 ' + cumm_var_to_plot + ' plot for '+Site_ID+'_'+versionID,size=20) 
    pl.savefig(mypathforResults+'/'+'Cummulative CO2 ' + cumm_var_to_plot + ' plot for '+Site_ID+'_'+versionID) 
    #pl.show()
    pl.close()    

def cummulative_H2O(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill,versionID):    
    print "Doing Cummulative H2O plots for "+Site_ID
 
    #First plot use ALL data and years
    ###################################
    #Create the Cummulative variables
    New_combined['Fc_Cumm']=New_combined.groupby([lambda x: x.year])['Fc_ustar'].cumsum()
    New_combined['Fc_Cumm']=New_combined['Fc_Cumm']
    New_combined['Fe_Cumm']=New_combined.groupby([lambda x: x.year])['Fe_Con'].cumsum()  
    New_combined['Precip_Cumm']=New_combined.groupby([lambda x: x.year])[Rain_Con_label_variable_to_fill].cumsum()
    
    #get to units of g.CO2.-m-2. multiply by 44 (MW of CO2) to get g CO2
    New_combined['Fc_Cumm']=New_combined['Fc_Cumm']*60*30*10000/1000000*44.01*12/44/1000000
    #Get from W.m-2 to mm
    New_combined['Fe_Cumm']=New_combined['Fe_Cumm']*60*30/1000000/2.45	#Get year for plot lables    

    #Group all data by Year
    tempdata_grouped=New_combined.groupby([lambda x: x.year])

    Num_plots=len(tempdata_grouped)
   
    cm = get_cmap('gist_rainbow')
    #Second plot Fe
    ax1=pl.figure(1, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    index=0
    for name, group in tempdata_grouped:
	plotyear=group.index[0].year
	index=index+1
	xdata1=group['Fe_Cumm']
	xdata2=group['Precip_Cumm']
	plotyear=group.index[0].year  
	color_rgb = cm(1.*index/Num_plots)
	ax1=pl.plot(xdata1.values, '-', linewidth=2,color=color_rgb,label='ET '+str(plotyear))
	ax1=pl.plot(xdata2.values, '--', linewidth=2,color=color_rgb,label='Precip '+str(plotyear))

	#index=0
	#for name, group in tempdata_grouped:
	    #plotyear=group.index[0].year
	    #index=index+1
	    #xdata1=group['Fe_Cumm']
	    #xdata2=group['Precip_Cumm']
	    #plotyear=group.index[0].year  
	    #color_rgb = cm(1.*index/Num_plots)
	    ##ax1=pl.plot(xdata1, '-', linewidth=2,color=color_rgb,label='ET '+str(plotyear))
	    ##ax1=pl.plot(xdata2, '--', linewidth=2,color=color_rgb,label='Precip '+str(plotyear))
	    #x_coord_et=xdata1[-1]
	    #y_coord_et=len(xdata1)-1
	    #x_coord_Precip=xdata2[-1]
	    #y_coord_Precip=len(xdata2)-1
	    #pl.annotate(name, xy=(x_coord_et, y_coord_et),  xycoords='data',
		        #xytext=(50, 0), textcoords='offset points',
		        #arrowprops=dict(arrowstyle="->")
		        #)	
	    #pl.annotate(name, xy=(x_coord_Precip, y_coord_Precip),  xycoords='data',
		        #xytext=(50, 0), textcoords='offset points',
		        #arrowprops=dict(arrowstyle="->")
		        #)	

    pl.legend(loc='upper left')
    pl.ylabel('Cummulative water flux (mm)')
    pl.xlabel('Time (30 minute periods)')
    pl.suptitle('Cummulative H2O plot for '+Site_ID+'_'+versionID,size=20) 
    pl.savefig(mypathforResults+'/'+'Cummulative H2O plot for '+Site_ID+'_'+versionID) 
    #pl.show()
    pl.close()   
    
def mintimeseries_plot(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,targets,output,item,versionID):    
    ANN_label=str(item+"_NN")    
    pl.plot( targets, 'b--' )
    pl.plot( output, 'k-' )
    pl.legend(('targets', 'output'))
    pl.xlabel('Time'); 
    pl.title('Outputs vs. target of trained network for '+item)
    pl.grid(True)
    pl.legend()
    pl.title('Tower vs ANN 30 min timeseries for '+item+' at ' +Site_ID)
    pl.ylabel(item + '('+units+')')
    pl.savefig(mypathforResults+'/'+'Tower vs ANN 30 min timeseries for '+item+' at ' +Site_ID+'_'+versionID) 
    #pl.show()
    pl.close()
    


def count_numbers(frame):
    return (frame.count())

def count_total(frame):
    return (frame.shape[0])


def do_nan_stats(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,versionID):

    for item in list_in:
	print "Doing Stats for Nans and pct filled"
	
	#Define the variable (item) for counts. 
	item_Con=item+'_Con'
	if item in ['Fc','Fe','Fh','Fg']:
	    item_Corr=item+'_NN'
	else:
	    item_Corr=item+'_Corr'
	
	#Do stats and counts for Sws
	Numbers_by_month= New_combined[[item,item_Con,item_Corr]].groupby([lambda x: x.year,lambda x: x.month]).apply(count_numbers)
	Numbers_by_month['total_n'] = New_combined[item_Con].groupby([lambda x: x.year,lambda x: x.month]).apply(count_total)  
	#Numbers_by_month['Pct_nan'] = (float(Numbers_by_month['Ta_EC'])/float(Numbers_by_month['total_n']))*100
	Numbers_by_month['Pct_nan']=Numbers_by_month.apply(lambda x: int(((np.float(x['total_n'])-np.float(x[item]))/np.float(x['total_n']))*100), axis=1)
	Numbers_by_month['Pct_notfilled']=Numbers_by_month.apply(lambda x: int(((np.float(x['total_n'])-np.float(x[item_Con]))/np.float(x['total_n']))*100), axis=1)
	Numbers_by_month.astype('int32')
	
	#Write out file
	Numbers_by_month.to_csv(mypathforResults+'/'+'Nan counts and Pct filled for '+item+' at ' +Site_ID+'_'+versionID+'.csv')	
	Numbers_by_month.to_pickle(mypathforResults+'/'+'Nan counts and Pct filled for '+item+' at ' +Site_ID+'_'+versionID+'.df') 	

def Dotables(mypathforResults,myBaseforResults,New_combined, Site_ID,subset_list_in,Ws_label,plot_freq,versionID):
    if plot_freq=="year":
	summary_df=New_combined.groupby([lambda x: x.year]).mean()
	summary_df_subset=New_combined[subset_list_in].groupby([lambda x: x.year]).mean()
    elif plot_freq=="month":
	summary_df=New_combined.groupby([lambda x: x.year,lambda x: x.month]).mean()
	summary_df_subset=New_combined[subset_list_in].groupby([lambda x: x.year,lambda x: x.month]).mean()	
    elif plot_freq=="week":
	summary_df=New_combined.groupby([lambda x: x.year,lambda x: x.week]).mean()
	summary_df_subset=New_combined[subset_list_in].groupby([lambda x: x.year,lambda x: x.week]).mean()
	
    #Write outputs for each site seperately
    summary_df.to_csv(myBaseforResults+"/"+"summary_for_freq_"+plot_freq+"_at_"+Site_ID+"_"+versionID+".csv")
    summary_df_subset.to_csv(myBaseforResults+"/"+"summary_subset_for_freq_"+plot_freq+"_at_"+Site_ID+"_"+versionID+".csv") 

###########################################################################################################
##                 START MAIN CODE
###########################################################################################################
def basic_diags(myBaseforResults,New_combined,Site_id,list_in,Ws_label,do_results,Rain_label_variable_to_fill,versionID):     
    global Site_ID
    Site_ID=Site_id
    
    print "============================================"
    print "DINGO: Starting Diagnostics"
    print "============================================"
    Rain_Con_label_variable_to_fill= Rain_label_variable_to_fill+'_Con'
    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/Diagnostics"):
	os.mkdir(myBaseforResults+"/Diagnostics")
    mypathforResults=myBaseforResults+"/Diagnostics"  
    
    number_of_inputs=len(list_in)
    
    startdate=New_combined.index[0]
    enddate=New_combined.index[len(New_combined)-1]
    
    list_in=['Fc','Fe','Fh','Fg','Ta','Ah',Rain_label_variable_to_fill,Ws_label,'Ts','Sws','Fsd','Fsu','Fld','Flu','Fn'] 
    do_nan_stats(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,versionID)
    
   
    ###############################
    #Do plots
    ###############################
    #Plots at specified fequencies.  Pass freq to function.  Can be 'day of year', 'week', 'month', 'year'
    #Plot timeseries of daily over all periods
    plot_freq = 'dayofyear'
    
    list_in=['Fc','Fe','Fh','Fg']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)
    list_in=['Fsd','Fsu','Fld','Flu','Fn']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)
    list_in=['Ta','Ah',Ws_label,'ps']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)    
    list_in=['Ts','Sws',Rain_label_variable_to_fill,'VPD']
    Doplots_at_daily(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)   
    
    #Note here that there are NEW variable names that will be dealt with in the function
    list_in=['BR','WUE','RUE','EBC','EF']
    Doplots_at_daily_other(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID) 

    plot_freq = 'week'
    
    list_in=['Fc','Fe','Fh','Fg']
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)
    list_in=['Fsd','Fsu','Fld','Flu','Fn']
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)
    list_in=['Ta','Ah',Ws_label,'ps']
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)    
    list_in=['Ts','Sws',Rain_label_variable_to_fill]
    Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)    

    plot_freq = 'month'
    
    #list_in=['Fc','Fe','Fh','Fg']
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)
    #list_in=['Fsu','Fsd','Fn']
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)
    #list_in=['Ta','Ah',Ws_label]
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)    
    #list_in=['Ts','Sws',Rain_label_variable_to_fill]
    #Doplots_at_freq(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID)       
    
    #Plot timeseries of monthly over all periods
    list_in=['Ta','Ah',Ws_label,'ps']
    Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,versionID)
    list_in=['Fc','Fe','Fh','Fg']
    Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,versionID)
    list_in=['Fsd','Fsu','Fld','Flu','Fn']
    Doplots_monthly_diff(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,versionID)    
    
    #Plot Nans and Pct data missing
    list_in=['Ta','Ah',Ws_label]	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,versionID)
    
    list_in=['Fc_ustar','Fe','Fh','Fg']	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,versionID)

    list_in=['Fsd','Fsu','Fld','Flu','Fn']	
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,versionID)    
    
    list_in=['Ts','Sws',Rain_label_variable_to_fill]
    Plot_Pct_nans(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,versionID)

    
    #Plots of EB closure
    ####################################################################################
    regressionEBclosure(mypathforResults,New_combined,Site_ID,startdate,enddate,versionID)
    
    # Do plots of Fre (Reichstein ustar approach vs Lasslop daytime Fre approach
    ####################################################################################
    #Define list of plots to make to process
    freq_list=["dayofyear", "week", "month"]
    list_in=["Fre","GPP","Fc"]	
    regressionLasslop(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in,versionID)    

    #Define list of plots to make to process
    freq_list=["dayofyear", "week", "month"]
    list_in=["Fre_noct","Fre_NN","Fre_Con"]	
    regressionFre(mypathforResults,New_combined,Site_ID,startdate,enddate,freq_list,list_in,versionID)   

    ##Plot time series of all 30 minute data
    #mintimeseries_plot(mypathforResults,predicted,observed,regress,item, Site_ID,units,targets,output,ANN_label,versionID)
    ##Plot regression of Tower versus ANN
    #regressionANN(mypathforResults,predicted,observed,regress,item, Site_ID,units,ANN_label,versionID)
    ##Plot diurnals for every second month 6 graphs
    #Doplots_diurnal(mypathforResults,New_combined,item, Site_ID,units,ANN_label,versionID)
    ##Plot timeseries of monthly over all periods
    #Doplots_monthly(mypathforResults,New_combined,item, Site_ID,units,ANN_label,versionID)
    
    #Do results plots
    if do_results==True:
	print "Starting Results"
	
	#Check for place to put results - does it exist? If not create
	if not os.path.isdir(myBaseforResults):
	    os.mkdir(myBaseforResults)
	#Then subdirectories
	if not os.path.isdir(myBaseforResults+"/Results"):
	    os.mkdir(myBaseforResults+"/Results")
	mypathforResults=myBaseforResults+"/Results"  	#Create directory for results

	#Do cummulative water and carbon plots
	cummulative_H2O(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill,versionID)
	
	cumm_var_to_plot="Fc_ustar"
	cummulative_CO2(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill,versionID,cumm_var_to_plot)
	cumm_var_to_plot="GPP"
	cummulative_CO2(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill,versionID,cumm_var_to_plot)
	cumm_var_to_plot="Fre"
	cummulative_CO2(mypathforResults,New_combined,Site_ID,startdate,enddate,Rain_Con_label_variable_to_fill,versionID,cumm_var_to_plot)
	
	
	#Do time series CO2 (GPP,Re plots)
	list_in=['Fc_ustar','GPP_Con','Fre_Con']
	plot_freq = 'dayofyear'
	Doplots_at_daily_carbon_umol(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID) 	
	Doplots_at_daily_carbon_g(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID) 	

	#Do time series Energy Balnce (Fn,Fe,Fh,Fg plots)
	list_in=['Fn_Con','Fe_Con','Fh_Con','Fg_Con']
	plot_freq = 'dayofyear'	
	Doplots_at_daily_EB(mypathforResults,myBaseforResults,New_combined, Site_ID,list_in,Ws_label,plot_freq,versionID) 
	
	#Do outputs of Tables.  Do this for monthly and annual.
	#Provide a subset list of variables you want.  The script will automatically generate an output of ALL variables as well
	#Do time series CO2 (GPP,Re plots).
	
	#test for CABLE reainfall elese table generation fails.  Just create and set to nans
	if ("Rainf_CABLE_mm" in New_combined.columns.values) == False: New_combined["Rainf_CABLE_mm"]=np.nan
	    
	
	subset_list_in=["Ta_Con","ps_Con","Ah_Con","VPD_Con","Precip_Con","Rainf_CABLE_mm","Fsd_Con","Fsu_Con","Flu_Con","Fld_Con","Fn_Con",
	                 "Fc","Fc_Con","Fc_ustar","Fe_Con","Fh_Con","Fg_Con","Fre_Con","GPP_Con","ustar_used","ustar_Reich_max", "ustar_Reich_var", "ustar_Barr", 
	                 "250m_16_days_EVI_new_smooth","250m_16_days_NDVI_new_smooth","Fpar_1km_new_smooth","Lai_1km_new_smooth","Gpp_1km_new_smooth"]
	#Define the freq at which to output.Can be "week", "month", "year"
	plot_freq = 'year'
	Dotables(mypathforResults,myBaseforResults,New_combined, Site_ID,subset_list_in,Ws_label,plot_freq,versionID) 
	
	#Define the freq at which to output.Can be "week", "month", "year"
	plot_freq = 'month'
	Dotables(mypathforResults,myBaseforResults,New_combined, Site_ID,subset_list_in,Ws_label,plot_freq,versionID) 
	
    print "============================================"
    print "DINGO: Finished Diagnostics " + Site_ID 
    print "============================================"	
    
    
    
    
    
