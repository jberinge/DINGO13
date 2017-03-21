############################################################################
#Dynamic INtegrated Gap-filling and partitioning for OzFlux (DINGO)
############################################################################
# Reichstein et al module
# For the Reichstein et al. (2005a) approach, the non-gap-filled data
#set is split into 6 temperature classes of the same sample size (according 
#to quantiles) and for each temperature class the set is split into 20 u*-classes
#according to Papale et al., (2006). The threshold is defined as the u*-class 
#where the night-time flux reaches more than 95% of the average flux at the 
#higher u*-classes. The threshold is only accepted if the temperature class
#and u* are not or only weakly correlated (|r| < 0.3). The final threshold
#is defined as the median of the thresholds of the (up-to) six temperature 
#classes. This procedure is applied to the entire dataset, giving a maximum, 
#but conservative u* threshold.  The maximum value is saved as 'ustar_Reich_max'.
#In addition, the u* threshold is calculated continuously using a 3 month 
#moving window to account for seasonal variation of vegetation structure
#
# Programmed by Jason Beringer Updated 5/5/2016
############################################################################

import pandas as pd
import math as math
import numpy as np
import os
import datetime as dt
import meteorologicalfunctions as metfuncs
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages
import time

from pylab import *
from ffnet._tests import runtest
from ffnet import ffnet, mlgraph, readdata, tmlgraph, imlgraph
from numpy import array
    
def regressionUSTAR(xdata,ydata,regress, Site_ID,DFyear,DFmonth):    
    graphtext1=str('slope      ' + str("{0:.2f}".format(regress[0])) +'\n' +
                   'intercept  ' + str("{0:.2f}".format(regress[1])) +'\n' +
                   'r-value    ' + str("{0:.2f}".format(regress[2])) +'\n' +
                   'p-value    ' + str("{0:.2f}".format(regress[3])) +'\n' +
                   'slope SE   ' + str("{0:.2f}".format(regress[4])) +'\n' )  
    pl.figtext(0.7,0.6,graphtext1, bbox=dict())
    
    line = regress[0]*xdata+regress[1]
    pl.plot(xdata,line,'r-',xdata,ydata,'o')    
    
    #pl.plot(x, y, linewidth = 2, label = 'regression line')
    pl.title('Ustar threshold versus Temperature from '+str(DFyear)+ ' '+str(DFmonth) +' for '+Site_ID)
    pl.xlabel('Ustar threshold for Temp Bin (m.s-1)')
    pl.ylabel('Mean Temperature of Bin (oC)')
    pl.savefig(mypathforResults+'/Ustar threshold versus Temperature from '+str(DFyear)+ ' '+str(DFmonth) +' for '+Site_ID) 
    #pl.show()
    pl.close()    

def func_halfsat(x, a, b,c):
    return ((a * x)/(b+x))+c
def func(x, a, b):
    return a+b/x

              
def dotempbin(x,**kwargs):
    #Calculate 20 ustar quantiles
    ustar_classes=pd.qcut(x['ustar'], 20)
    ustar_grouped=x.groupby(ustar_classes)
    xdata= ustar_grouped['ustar'].mean()
    ydata= ustar_grouped['Fc'].mean() 
    MeanBinTemp= x['Ts_Con'].mean() 
    try:
        #print 'start'
        popt, pcov = curve_fit(func, xdata, ydata,maxfev=5000) 
        #print "finished fitting" 
        #print"calculate ustar at 95% flux contribution"
        #print "last ustar label",xdata[19]
        flux100pct=func(xdata[19], *popt)
        flux95pct=func(xdata[19], *popt)*.95
        #print "Flux at highest ustar ", flux100pct
        popt2, pcov2 = curve_fit(func, ydata, xdata,maxfev=5000)
        ustar_threshold=func(flux95pct, *popt2)
        #print "Ustar at 95% flux  ", ustar_threshold
        #print'double check % ', func(xdata[19], *popt)
    except:
        ustar_threshold=np.nan
 
    return pd.Series({'ustar_threshold': ustar_threshold, 'MeanBinTemp': MeanBinTemp})

###########################################################################################################
##                 START MAIN CODE
###########################################################################################################


def fill_ustar(x):

    DFyear=x.index[0].year
    DFmonth=x.index[0].month
    fillvalue=float( ResultsDF['ustar_final'][ResultsDF['year']==DFyear][ResultsDF['month']==DFmonth])
    x['ustar_Reich_var']=fillvalue
    return x

def ustar_filter(New_combined,myBaseforResults,Site_ID):
    print "========================================================"
    print "DINGO: STARTING REICHSTEIN USTAR THRESHOLD CALCULATION"
    print "========================================================"
    
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/Ustar"):
        os.mkdir(myBaseforResults+"/Ustar")
    global mypathforResults
    mypathforResults=myBaseforResults+"/Ustar"  
    
    #Define date range of file
    startdate_initial=New_combined.index[0]
    enddate_initial=New_combined.index[len(New_combined)-1]
    
    #Create empty list for results
    Resultslist=[]
    for index,group in New_combined[['Ts_Con','ustar']].groupby([lambda x: x.year,lambda x: x.month]):
        DFyear=int(index[0])
        DFmonth=int(index[1])
        
        #Create a 3 month window.  This will be a moving window as this gets indexed by the bygroup which is month
        middate=dt.datetime(DFyear,DFmonth,1,0,0,0)
        enddate=middate+ dt.timedelta(days=60)
        startdate=middate- dt.timedelta(days=30)      
        print "Processing ustar for date ", middate
        # Choose only data that has BOTH Fc and Ustar
        ustar_DF=New_combined[['Ts_Con','Fc','ustar','Fsd']].dropna(axis=0,how='any')
        
        #Choose date and only nightime where Fsd is less than 10 Wm-2
        ustar_DF=ustar_DF[startdate:enddate][ustar_DF['Fsd']<10]
        #Check number of samples > say 400 min for binning
        if len(ustar_DF)>400:
            #Calculate 6 temperature BINS quantiles equal samples
            temperature_classes=pd.qcut(ustar_DF['Ts_Con'], 6)#,labels=False)
            #Then groupby temperature classes
            temp_grouped=ustar_DF.groupby(temperature_classes)
            #Then apply the function 
            try:
                results=temp_grouped[['Fc','ustar','Ts_Con']].apply(dotempbin, args=(DFyear,DFmonth)).reset_index(drop=True)
                #Determine the final threshold for the period by taking the median ustar value in the temperature bins
                ustar_final=results['ustar_threshold'].median()
                #########################################################
                # Do linear regression of ustar and Temperature
                #########################################################
                #Check result to see if T and ustar are correlated. 
                #Do linear regression - returns  slope, intercept, r_value, p_value, std_err
                regress = stats.linregress(results['ustar_threshold'],results['MeanBinTemp'])
                regressionUSTAR(results['ustar_threshold'],results['MeanBinTemp'],regress, Site_ID,DFyear,DFmonth)    
                ustar_T_rsqu=regress[2]          
            except:
                ustar_final=np.nan
                ustar_T_rsqu=np.nan            
        else:
            ustar_final=np.nan
            ustar_T_rsqu=np.nan
            
        Resultslist.append([DFyear,DFmonth,ustar_final,ustar_T_rsqu])
        
    #########################################################
    # Outputs
    #########################################################
    #create a dataframe for the results
    global ResultsDF
    ResultsDF=pd.DataFrame(Resultslist,columns=['year', 'month','ustar_final', 'ustar_T_rsqu'])  
    
    Ustar_mean=ResultsDF['ustar_final'].mean()
    Ustar_max=ResultsDF['ustar_final'].max()
    print "The final mean Ustar for "+Site_ID+" is ",Ustar_mean
    print "The final max Ustar for "+Site_ID+" is ",Ustar_max
    
    ResultsDF_grouped=ResultsDF.groupby('month')
    ResultsDF_monthly=ResultsDF_grouped['ustar_final'].agg([ np.mean, np.std]).reset_index()
    
    ##################################
    #Final plot  before gap filling
    ##################################
    print "Doing Ustar plots"
    pdf = PdfPages(mypathforResults+'/Monthly Ustar threshold and SD for '+Site_ID+'.pdf')
    plt.plot(ResultsDF_monthly['month'], ResultsDF_monthly['mean'], linestyle="dashed", marker="o", color="green")
    #plot only errorbars
    plt.errorbar(ResultsDF_monthly['month'], ResultsDF_monthly['mean'], ResultsDF_monthly['std'], linestyle="None", marker="None", color="green")
    #plt.ylim(min(ResultsDF_monthly['mean']), (Ustar_max*1.05))
    plt.xlim(0,13)
    #create horizontal lines
    plt.axhline(y=Ustar_mean,linewidth=2, color='r')
    plt.axhline(y=Ustar_max,linewidth=2, color='b')
    # place a text box in upper left in axes coords
    textstr=str('Max  ustar (Blue) ' + str("{0:.2f}".format(Ustar_max))+'\n' +
                'Mean ustar (Red)  ' + str("{0:.2f}".format(Ustar_mean))             )
    plt.figtext(0.6,0.75,textstr, bbox=dict())
    
    fmt = '%b %Y '
    textstr2=str("From date  " +startdate_initial.strftime(fmt)+'\n' +
                 " to " +enddate_initial.strftime(fmt))
    plt.figtext(0.2,0.75,textstr2, bbox=dict())
    #labels
    plt.ylabel("Ustar threshold(m.s-1)")
    plt.xlabel("Month")
    #title
    plt.title("Monthly Ustar threshold and SD for "+Site_ID)
    #show plot
    plt.savefig(pdf, format='pdf')
    #plt.show() 
    pdf.close()
    
    #Final plot  before gap filling
    pdf = PdfPages(mypathforResults+'/Timeseries of monthly Ustar threshold and SD for '+Site_ID+'.pdf')
    plt.plot(ResultsDF['ustar_final'], linestyle="dashed", marker="o", color="blue")
    #create horizontal lines
    plt.axhline(y=Ustar_mean,linewidth=2, color='r')
    plt.axhline(y=Ustar_max,linewidth=2, color='b')
    # place a text box in upper left in axes coords
    textstr=str('Max  ustar (Blue) ' + str("{0:.2f}".format(Ustar_max))+'\n' +
                'Mean ustar (Red)  ' + str("{0:.2f}".format(Ustar_mean))             )
    plt.figtext(0.6,0.75,textstr, bbox=dict())
    fmt = '%b %Y '
    textstr2=str("From date  " +startdate_initial.strftime(fmt)+'\n' +
                 " to " +enddate_initial.strftime(fmt))
    plt.figtext(0.2,0.75,textstr2, bbox=dict())
    #labels
    plt.ylabel("Ustar threshold(m.s-1)")
    plt.xlabel("Date")
    #title
    plt.title("Timeseries of monthly Ustar threshold and SD for "+Site_ID)
    #show plot
    plt.savefig(pdf, format='pdf')
    #plt.show() 
    pdf.close()
    
    #Fill nans with the mean ustar
    ##############################
    ResultsDF['ustar_final'].fillna(Ustar_mean,inplace=True)
    ResultsDF['ustar_T_rsqu'].fillna(99,inplace=True)
    #Now output to the original dataframe
    #create a new column first
    New_combined['ustar_Reich_var']=np.nan
    Newest_combined=New_combined.groupby([lambda x: x.year,lambda x: x.month]).apply(fill_ustar)  
    #If value still missing assign the mean
    Newest_combined['ustar_Reich_var'].fillna(Ustar_mean)
 
    #Also have a conservative output (i.e. the max ustar over the whole period then 
    # make the ustar threshold equal to the max value
    Newest_combined['ustar_Reich_max'] = Ustar_max
    print "========================================================"
    print "DINGO: FINISH REICHSTEIN USTAR THRESHOLD CALCULATION"
    print "========================================================"
    return Newest_combined
       
       
       