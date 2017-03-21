import constants as c
import datetime
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import meteorologicalfunctions as mf
import numpy
import qcck
import qcio
import qcts
import qcutils
import qcplot
import sys
import os
import pandas as pd
import matplotlib.dates as mdt
from configobj import ConfigObj



def plot_fingerprint(data,xlabel=None,sub=[1,1,1],extent=None,ticks=None):
    loc,fmt = qcplot.get_ticks(datetime.datetime.fromordinal(sd),datetime.datetime.fromordinal(ed))
    fmt = mdt.DateFormatter('%Y/%m')
    ax = plt.subplot(sub[0],sub[1],sub[2])
    
    cdict1 = {'red':   ((0.0, 0.0, 0.0),
	                (0.4, 0.2, 0.0),
                        (0.6, 0.8, 1.0),
	                (1.0, 1.0, 1.0)),
    
	      'green': ((0.0, 0.0, 0.0),
	                (0.4, 0.2, 0.0),
                        (0.6, 0.8, 1.0),
	                (1.0, 1.0, 1.0)),
    
	      'blue':  ((0.0, 0.0, 0.0),
	                (0.4, 0.0, 0.0),
                        (0.6, 0.0, 0.0),
	                (1.0, 0.0, 0.0))
	      }    

    cdict2 = {'blue':   ((0.0, 0.0, 0.0),
	                (0.666, 0.0, 1.0),
	                (1.0, 0.1, 1.0)),
    
	      'red': ((0.0, 0.0, 0.0),
	                (1.0, 0.0, 0.0)),
    
	      'green':  ((0.0, 0.0, 0.1),
	                (0.666, 1.0, 0.0),
	                (1.0, 0.0, 0.0))
	      }

    blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)
    blue_red2 = LinearSegmentedColormap('BlueRed2', cdict2)

    #To produce figures for caitlins plots we make this true.
    caitlin=False
    if caitlin==True: plt.imshow(data,extent=extent,aspect='auto',origin='lower',vmin=-30, vmax=15, cmap=blue_red2)
    if caitlin==False: plt.imshow(data,extent=extent,aspect='auto',origin='lower')
    
    ax.yaxis.set_major_locator(loc)
    ax.yaxis.set_major_formatter(fmt)
    plt.colorbar(orientation='horizontal',fraction=0.02,pad=0.075,ticks=ticks)
    plt.xticks([0,6,12,18,24])
    if xlabel != None: plt.xlabel(xlabel)
    if sub[2] != 1: plt.setp(ax.get_yticklabels(), visible=False)


def fingerprint_plots(myBaseforResults,New_combined,Site_ID,CFname,versionID):
    print "============================================="
    print "DINGO: STARTING FINGERPRINT PLOTS "
    print "============================================="    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
	os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/Fingerprints"):
	os.mkdir(myBaseforResults+"/Fingerprints")
    mypathforResults=myBaseforResults+"/Fingerprints"  
		
    #Returns CF object that we can get info from later
    cf = ConfigObj(CFname) 	    
    
    
    PlotWidth = float(cf['General']['PlotWidth'])
    PlotHeight = float(cf['General']['PlotHeight'])
    
    # get the datetime series
    DateTime = New_combined.index
    
    #Find the start and end date of the series
    global sd, ed
    sd = datetime.datetime.toordinal(DateTime[0])
    ed = datetime.datetime.toordinal(DateTime[-1])
    TitleStr = Site_ID +' from '+str(DateTime[0])+' to '+str(DateTime[-1])
    
    for nFig in cf['Plots'].keys():
	n = 0
	fig = plt.figure(nFig,figsize=[15,10])
	plt.figtext(0.5,0.95,TitleStr,horizontalalignment='center')
	SeriesList = qcutils.GetPlotVariableNamesFromCF(cf,nFig)
	list_string=''
	for z in SeriesList:
	    #create string list
	    list_string=list_string+' '+z	
	
	nPlots = len(SeriesList)
	for ThisOne in SeriesList:
	    n += 1
	    VarName = qcutils.GetAltNameFromCF(cf,ThisOne)
	    ticks = qcutils.GetcbTicksFromCF(cf,ThisOne)
	    lower, upper = qcutils.GetRangesFromCF(cf,ThisOne)
	    
	    #Round to nearest whole day at start and end
	    #sd_trunc=datetime.date(DateTime[0].year,DateTime[0].month,DateTime[0].day)+datetime.timedelta(days=1)
	    #ed_trunc=datetime.date(DateTime[-1].year,DateTime[-1].month,DateTime[-1].day)-datetime.timedelta(days=1)

	    sd_trunc=datetime.datetime(DateTime[0].year,DateTime[0].month,DateTime[0].day,0,0,0,0)+datetime.timedelta(days=1)
	    ed_trunc=datetime.datetime(DateTime[-1].year,DateTime[-1].month,DateTime[-1].day,0,0,0,0)-datetime.timedelta(days=1)
	    
	    New_combined=New_combined[sd_trunc:ed_trunc]
	    New_combined['DT']= New_combined.index
	    New_combined['DTindex']=(New_combined.loc[:,'DT'].apply(lambda x:int(datetime.datetime.strftime(x,'%Y')))*1000 + 
		                     New_combined.loc[:,'DT'].apply(lambda x:int(datetime.datetime.strftime(x,'%j')))     )
	    
	    by = lambda x: lambda y: getattr(y, x)
	    #data_30minDF=New_combined[ThisOne].groupby([['DTindex'],by('hour')]).mean()
	    data_30minDF=New_combined[ThisOne].groupby([by('year'),by('dayofyear'),by('hour')]).mean().unstack()
	    data_30min= numpy.array(data_30minDF)
	    nPerDay=len(New_combined[ThisOne].groupby([by('hour')]))
	    #nDays=len(New_combined[ThisOne].groupby([by('year'),by('dayofyear')]))
	    nDays=len(New_combined[ThisOne].groupby([by('year'),by('dayofyear')]))	    
	    data_daily = data_30min.reshape(nDays,nPerDay)

	    units=''
	    label = VarName + ' (' + units + ')'
	    plot_fingerprint(data_daily,xlabel=label,sub=[1,nPlots,n],extent=[0,24,sd,ed],ticks=ticks)
    
	pngname = mypathforResults+'/Fingerprint plots for '+list_string+' at '+Site_ID + '_' + versionID+'.png'
	fig.savefig(pngname,format='png')
    print "============================================="
    print "DINGO: FINISHED FINGERPRINT PLOTS "
    print "============================================="    
    plt.draw()
    #plt.show()
