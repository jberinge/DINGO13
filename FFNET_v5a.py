############################################################################
#Code for gap filling of flux data using ANN
#
#Following the meteorological drivers the fluxes of water (Fe), net ecosystem exchange of carbon (Fc), ecosystem respiration (Fre), sensible heat (Fh), and
#ground heat (Fg) are gap-filled using an artificial neural network (ANN) approach.  Essential background on ANNs can be found in Basheer et al. (2000) and 
#useful examples of its application can be found in (Moffat et al. (2010) and Bryant and Shreeve  (2002).  We implemented the ANN as follows using 
#gap-filling performed using the ffnet package (http://ffnet.sourceforge.net/) that is a fast and easy-to-use feed-forward neural network. The basic
#model we implemented is a multiple layer network with two layers of 24 and 16 nodes respectively, which creates a multilayer network that is fully connected.
#Training is done using a truncated Newton algorithm (TNC) to minimise a function with variables subject to bounds using gradient information.  The type of 
#ANN model and training algorithms can be altered in the code,
#
#
# Programmed by Jason Beringer (Last updated 5/5/2016)
#
#Note this module uses the ffnet ANN ((http://ffnet.sourceforge.net).
#To install on Windows I used Anaconda python distribution
#Download binary from http://www.lfd.uci.edu/~gohlke/pythonlibs/
#Use pip in anaconda to instal.
#
############################################################################

import pandas as pd
import numpy as np
import os
import datetime as dt
import time
import pylab as pl
import meteorologicalfunctions as metfuncs
from scipy import stats

from pylab import figure, ioff, clf, contourf, ion, draw, show
from ffnet import ffnet, mlgraph, readdata, tmlgraph, imlgraph
from numpy import array

def Doplots_diurnal_monthly(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,list_out,index_str,is_this_all):
    for index, item in enumerate(list_out):
        print "Doing diurnal plot for month and index ",index_str
        #Do Diurnal Plots for all 12 months
        #create an X axis series for all 24 hours
        t = np.arange(1, 25, 1)
        NN_label=str(item+"_NN")
        if is_this_all == True: NN_label=str(item+"_NN_all")
	Plottemp = PlottingDF[[NN_label,item]]

                
        fig=pl.figure(index, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('Monthly ANN ensemble diurnal average for variable '+item+' at '+Site_ID+ ' index '+index_str)
        pl.subplots_adjust(top=100)
       
        
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
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 
        
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'r',label=NN_label)
        pl.ylabel('Flux')    
     	pl.legend(loc='best')
	
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
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False
	
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'r',label=NN_label)
        pl.ylabel('Flux')      
    	pl.legend(loc='best')
	
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
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 
	
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'r',label=NN_label)
        pl.ylabel('Flux')      
	pl.legend(loc='best')
		
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
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 	    
	    
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'r',label=NN_label)
        pl.ylabel('Flux')  
	pl.legend(loc='best')
		
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
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 	    
	    
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'r',label=NN_label)
        pl.ylabel('Flux')  
	pl.legend(loc='best')
		
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
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 	    
	    
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'r',label=NN_label)
        pl.ylabel('Flux')  
	pl.legend(loc='best')
	
	pl.tight_layout()  
        pl.savefig(mypathforResults+'/Monthly ANN ensemble diurnal average for variable '+item+' at '+Site_ID+ ' index '+index_str)
        #pl.show() 
        pl.close(index)
	pl.close()
	time.sleep(2)

def Doplots_diurnal(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,list_out,index_str,frequency):
    for index, item in enumerate(list_out):
        print "Doing diurnal plot for month and index ",index_str
        #Do Diurnal Plots for all 12 months
        #create an X axis series for all 24 hours
        t = np.arange(1, 25, 1)
        NN_label=str(item+"_NN")
	if frequency == 'all': NN_label=str(item+"_NN_all")
        Plottemp = PlottingDF[[NN_label,item]]
        #Plottemp = PlottingDF[[NN_label,item]].dropna(how='any')
                
        figure(2)
        pl.title('Diurnal '+item)
        try:
            xdata1a=Plottemp[item].groupby([lambda x: x.hour]).mean()
            plotxdata1a=True
        except:
            plotxdata1a=False
        try:
            xdata1b=Plottemp[NN_label].groupby([lambda x: x.hour]).mean()
            plotxdata1b=True
        except:
            plotxdata1b=False 
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 	    
	    
        if plotxdata1a==True:
            pl.plot(t,xdata1a,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(t,xdata1b,'g',label=NN_label)
	pl.legend(loc='best')   
        pl.ylabel('Flux')  
        pl.suptitle('ANN ensemble diurnal average for variable '+item+' at '+Site_ID+ ' index '+index_str,fontsize=14)
        pl.subplots_adjust(top=0.85)
        pl.tight_layout()  
        pl.savefig(mypathforResults+'/ANN ensemble diurnal average for variable '+item+' at '+Site_ID+ ' index '+index_str)
        #pl.show() 
        pl.close(2)
	time.sleep(2)
	
def Doplots_monthly(mypathforResults,PlottingDF,variable_to_fill, Site_ID,units,list_out,index_str,frequency):   
    for index, item in enumerate(list_out):
	
	ANN_label=str(item+"_NN")     #Do Monthly Plots
        print "Doing  Monthly  plot index", index_str
        #t = arange(1, 54, 1)
        NN_label=str(item+"_NN")
	if frequency == 'all': NN_label=str(item+"_NN_all")
        Plottemp = PlottingDF[[NN_label,item]]
        #Plottemp = PlottingDF[[NN_label,item]].dropna(how='any')
        fig=pl.figure(3, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
        pl.title('ANN and Tower plots by year and month for variable '+item+' at '+Site_ID+ ' index '+index_str)
    
        try:
            xdata1a=Plottemp[item].groupby([lambda x: x.year,lambda x: x.week]).mean()
            plotxdata1a=True
        except:
            plotxdata1a=False
        try:
            xdata1b=Plottemp[NN_label].groupby([lambda x: x.year,lambda x: x.week]).mean()
            plotxdata1b=True
        except:
            plotxdata1b=False
	if (len(xdata1a)==0): plotxdata1a=False
	if (len(xdata1b)==0): plotxdata1b=False 	    
	    
        if plotxdata1a==True:
            pl.plot(xdata1a.values,'b',label=item) 
        if plotxdata1b==True:
            pl.plot(xdata1b.values,'DarkOrange',label=NN_label)
        pl.ylabel('Flux')    
        pl.xlabel('Year - Month')       
        pl.legend()
        pl.savefig(mypathforResults+'/ANN and Tower plots by year and month for variable '+item+' at '+Site_ID+ ' index '+index_str)
        #pl.show()
        try: 
	    time.sleep(2)
	    pl.close(3)
	    time.sleep(2)
	except:
	    print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
	    print "&&&&&&&&&& Graph aborted  Monthly  plot index " + 'ANN and Tower plots by year and month for variable '+item+' at '+Site_ID+ ' index '+index_str
	    pass

	
def regressionANN(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,list_out,index_str):    
    for index, item in enumerate(list_out):
	fig=pl.figure(4, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	ANN_label=str(item+"_NN") 
        graphtext1=str('slope      ' + str("{0:.2f}".format(regress[index][0])) +'\n' +
                       'intercept  ' + str("{0:.2f}".format(regress[index][1])) +'\n' +
                       'r-value    ' + str("{0:.2f}".format(regress[index][2])) +'\n' +
                       'p-value    ' + str("{0:.2f}".format(regress[index][3])) +'\n' +
                       'slope SE   ' + str("{0:.2f}".format(regress[index][4])) +'\n' +
                       'estim. SE  ' + str("{0:.2f}".format(regress[index][5]))         )  
        pl.figtext(0.7,0.6,graphtext1, bbox=dict())
        pl.plot(observed[:,index], predicted[:,index], 'o', label='targets vs. outputs')    
        slope = regress[index][0]; intercept = regress[index][1]
    
        x = np.linspace(min(observed[:,index]),max(observed[:,index]))
        y = slope * x + intercept
        pl.plot(x, y, linewidth = 2, label = 'regression line')
        pl.legend()
        pl.title('Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str)
        pl.xlabel('Tower ' + '('+units+')')
        pl.ylabel('ANN ' + '('+units+')')
        pl.legend(shadow=True, fancybox=True,loc='best')
        
        pl.savefig(mypathforResults+'/'+'Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str)
        #pl.show()
        pl.close(4)    
	time.sleep(2)

def regressionANN2(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,list_out,index_str):    
    for index, item in enumerate(list_out):
	
	ANN_label=str(item+"_NN") 

	predicted1=np.reshape(observed[:,index],(len(predicted)))
	observed1=np.reshape(predicted[:,index],(len(observed)))
	
	#In this function calculate the linear regression independantly rather than using the NN stats.
	#Use slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	slope, intercept, r_value, p_value, std_err = stats.linregress(predicted1,observed1)
    
	graphtext1=str('slope      ' + str("{0:.2f}".format(slope)) +'\n' +
	               'intercept  ' + str("{0:.2f}".format(intercept)) +'\n' +
	               'r-value    ' + str("{0:.2f}".format(r_value)) +'\n' +
	               'p-value    ' + str("{0:.2f}".format(p_value)) +'\n' +
	               'slope SE   ' + str("{0:.2f}".format(std_err))       )  

	figure(5)        
	fig=pl.figure(4, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
	pl.figtext(0.7,0.6,graphtext1, bbox=dict())
        pl.plot(observed[:,index], predicted[:,index], 'o', label='targets vs. outputs')    
        slope = regress[index][0]; intercept = regress[index][1]
    
        x = np.linspace(min(observed[:,index]),max(observed[:,index]))
        y = slope * x + intercept
        pl.plot(x, y, linewidth = 2, label = 'regression line')
        pl.legend()
        pl.title('Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str)
        pl.xlabel('Tower ' + '('+units+')')
        pl.ylabel('ANN ' + '('+units+')')
        pl.legend(shadow=True, fancybox=True,loc='best')
        
        pl.savefig(mypathforResults+'/'+'Tower vs ANN for '+item+' at ' +Site_ID+ ' index '+index_str)
        #pl.show()
        pl.close(4)    
	time.sleep(2)
	
def mintimeseries_plot(mypathforResults,predicted,observed,regress,variable_to_fill, Site_ID,units,targets,output,list_out,index_str):    
    for index, item in enumerate(list_out):
        ANN_label=str(item+"_NN")    
	fig=pl.figure(5, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
        pl.plot( targets[:,index], 'b--' )
        pl.plot( output[:,index], 'k-' )
        pl.legend(('targets', 'output'))
        pl.xlabel('Time'); 
        pl.title('Outputs vs. target of trained network for '+item)
        pl.grid(True)
        pl.legend()
        pl.title('Tower vs ANN 30 min timeseries for '+item+'  at ' +Site_ID+ ' index '+index_str)
        pl.ylabel(item + '('+units+')')
        pl.savefig(mypathforResults+'/'+'Tower vs ANN 30 min timeseries for '+item+'  at ' +Site_ID+ ' index '+index_str)
        #pl.show()
        pl.close(5)
	time.sleep(2)
	
def ANN_gapfill_func(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,index_str,is_this_all,ANN_label_all,ANN_label,frequency,Use_Fc_Storage):

    ###########################################################################################################
    ##                 START MAIN CODE
    ###########################################################################################################
    if 'Fc' in list_out:
        units="umol.m-2.s-1"
    elif ('Fe' or 'Fh' or 'Fg') in list_out:
        units="W.m-2"
    else:
        units=" "
    
    ###### User-set IO file locations ######
    
    print "Starting ANN gap filling"
    #Check for place to put results - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/ANN"):
        os.mkdir(myBaseforResults+"/ANN")
    mypathforResults=myBaseforResults+"/ANN"  
      
    #We need to update VPD for input here so also need e and es
    # Calculate vapour pressure from absolute humidity and temperature
    #  Ah - absolute humidity, g/m3
    #  Ta - air temperature, C
    New_combined['VPD_Con']=(metfuncs.es(New_combined['Ta_Con']))-(metfuncs.vapourpressure(New_combined['Ah_Con'],New_combined['Ta_Con']))
    
    number_of_inputs=len(list_in)
    number_of_outputs=len(list_out)
    #startdate=dt.date(2008,7,1)
    #enddate=dt.date(2008,8,1)
    alllist=list_in + list_out
    xnow=New_combined[alllist]                         #[startdate:enddate]
    xnow=xnow.dropna(how='any')
    #Drop nans and missing values so that Good data only is used in the training
    xarray=np.array(xnow.dropna().reset_index(drop=True))
    #Define inputs and targets for NN from DF
    inputs =  xarray[:, :number_of_inputs] #first 2 columns
    lastcolums=(-1*number_of_outputs)
    targets = xarray[:, lastcolums:] #last column
    
    # Generate standard layered network architecture and create network
    #different network architectures avaiable
    #conec = mlgraph((number_of_inputs,24,16,number_of_outputs))  # Creates standard multilayer network architecture
    conec = tmlgraph((number_of_inputs,24,16,number_of_outputs))  # Creates multilayer network full connectivity list
    #conec = imlgraph((number_of_inputs,24,16,number_of_outputs))  # Creates multilayer architecture with independent outputs
    
    net = ffnet(conec)
    
    print "TRAINING NETWORK..."
    net.train_tnc(inputs, targets, maxfun = iterations, messages=1)
    #net.train_rprop(inputs, targets, maxiter=iterations)
    #net.train_momentum(inputs, targets, maxfun = iterations, messages=1)
    #net.train_genetic(inputs, targets, maxfun = iterations, messages=1)
    #net.train_cg(inputs, targets, maxfun = iterations, messages=1)
    #net.train_bfgs(inputs, targets, maxfun = iterations, messages=1)
    
    
    # Test network
    print "TESTING NETWORK..."
    output, regression = net.test(inputs, targets, iprint = 0)
    
    print "R-squared:           %s  " %str(regression[0][2])
    #print "max. absolute error: %s  " %str(abs( array(output).reshape( len(output) ) - array(targets) ).max())
    output, regress = net.test(inputs, targets)
    
    #Create array for results. Then loop through elements on the original data to predict the ANN value
    predicted=np.empty((len(xarray),number_of_outputs))
    observed=np.empty((len(xarray),number_of_outputs))
    for index,rowdata in enumerate(xarray):
        predicted[index]=net([rowdata[0:number_of_inputs]])
        observed[index]=np.array(rowdata[-1.0*number_of_outputs : ])
	#observed[index]=np.array(rowdata[(-1.0*number_of_outputs)])
	
    ############################################    
    # Generate output and return new variables
    ############################################
    #Create a new variable called '_NN'
    for index, item in enumerate(list_out):
	ANN_label=str(item+"_NN")
	ANN_label_all=str(item+"_NN_all")
	if is_this_all == True:
	    New_combined[ANN_label_all]=net.call(New_combined[list_in])[:,index] 
	else:
	    New_combined[ANN_label]=net.call(New_combined[list_in])[:,index]    
    
    for index, item in enumerate(list_out):   
        #####################################################
        #  Plots 
        #####################################################
        #Plot time series of all 30 minute data
        mintimeseries_plot(mypathforResults,predicted,observed,regress,item, Site_ID,units,targets,output,list_out,index_str)
        #Plot regression of Tower versus ANN
        regressionANN2(mypathforResults,predicted,observed,regress,item, Site_ID,units,list_out,index_str)
        #Plot diurnals for every second month 6 graphs - only when enough months so all or annual
        if frequency=="all" or frequency=="annual" or is_this_all==True:
	    Doplots_diurnal_monthly(mypathforResults,New_combined,item, Site_ID,units,list_out,index_str,is_this_all)
        #Plot diurnals for every second month 6 graphs
        Doplots_diurnal(mypathforResults,New_combined,item, Site_ID,units,list_out,index_str,frequency)	
        #Plot timeseries of monthly over all periods
        Doplots_monthly(mypathforResults,New_combined,item, Site_ID,units,list_out,index_str,frequency)
                        
    return New_combined

###########################################################################################################
##                 START MAIN CODE
###########################################################################################################
#Not efficient but have written a wrapper around to get the grouping.
#Start by dividing the data up into slices.  Here can be all, annual, monthly or a column with a categorical variable (i.e. stage of LULCC or burning)

def ANN_gapfill(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,frequency,Use_Fc_Storage):  
    
    #here if the configuration file option to include Fcstorage then apply that here
    #This is the first time that Fc is used so we make a new variable.  Later we just
    #write that to the 'Fc_Con' variable.
    #ONly do this when this ANN is called wit Fc in the ouput
       
    if (Use_Fc_Storage=='Yes') and ('Fc'  in list_out):
	New_combined['Fc_inc_store']=New_combined['Fc']+New_combined['Fc_storage']
	list_out=['Fc_inc_store']
	
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
	
    #create and reset the variables
    for index, item in enumerate(list_out):
	ANN_label=str(item+"_NN")
	New_combined[ANN_label]=np.nan
	ANN_label_all=str(item+"_NN_all")
	New_combined[ANN_label_all]=np.nan
    
    if frequency == 'annual':
	New_combined_grouped=New_combined.groupby([lambda x: x.year])
    
    elif frequency == 'monthly':
	New_combined_grouped=New_combined.groupby([lambda x: x.year,lambda x: x.month])
    elif frequency=="all":
	pass
    else:
	New_combined_grouped=New_combined.groupby([frequency])
    
    #Do the ANN for the grouped analysis first. Then do 'ALL' data.
    #Need to fill in blocks with the annual ANN for blocks where there may not be enough data points
    #start with new column.  #Create a new variable called '_NN'
    #Fill with nans first
    #Only do for GROUPED data otherwise skip
    if frequency!="all":
	i=0
	for index,group in New_combined_grouped:
	    is_this_all=False
	    #Define a string for the index.  if it is monthly then need to make a string includes year and month 2 dimensions
	    #Otherwise just use the DF index 1 dimension.
	    if frequency=='monthly':
		DFyear=int(index[0])
		DFmonth=int(index[1])  
		index_str=str(DFyear*1000+DFmonth)
	    else:
		index_str=str(index)
	    print "Performing ANN on index ",index_str
	    
	    #Here try ANN on the block.  If the block doesnt have enough values it will return an error and that 
	    #block is effectively skipped and will be blank
	    try:
		temp = 	ANN_gapfill_func(myBaseforResults,group,Site_ID,list_in,list_out,iterations,index_str,is_this_all,ANN_label_all,ANN_label,frequency,Use_Fc_Storage)
	    except:
		print "***************************Something wrong pass back unchanged for index", index_str
		temp = group
		
		pass
	    #after the first block then start to concatenate parts together
	    if i==0:
		temp_concat=temp
	    else:
		temp_concat=pd.concat([temp_concat,temp],axis=0)
	    i=i+1
    
    #NOW do for  'ALL' data.
    #Need to fill in blocks with the annual ANN for blocks where there may not be enough data points
    #start with new column.  #Create a new variable called '_NN'  #Fill with nans first
    #Define a string for the index.
    index_str='all'
    is_this_all=True
    
    print "Performing ANN on index ",index_str
    #If using the grouped data then use the concatentaed file otherwise if using ALL then use the original DF
    if frequency=="all":
	temp_concat = ANN_gapfill_func(myBaseforResults,New_combined,Site_ID,list_in,list_out,iterations,index_str,is_this_all,ANN_label_all,ANN_label,frequency,Use_Fc_Storage)
    else:
	temp_concat = ANN_gapfill_func(myBaseforResults,temp_concat,Site_ID,list_in,list_out,iterations,index_str,is_this_all,ANN_label_all,ANN_label,frequency,Use_Fc_Storage)
    #Now put together the results from the groups and the all ANN
    #If user choose all data then simply copy 'all' column to the NN column
    #Otherwise keep column from groups and check for nans and fill them with the 'all' column
    for index, item in enumerate(list_out):
	ANN_label=str(item+"_NN")
	ANN_label_all=str(item+"_NN_all")
	if frequency == 'all':
	    temp_concat[ANN_label]=temp_concat[ANN_label_all]
	else:
	    temp_concat[ANN_label][temp_concat[ANN_label].isnull()]=temp_concat[ANN_label_all]
	 
     
    ############################################    
    # Generate output and return new variables
    ############################################
    #Create a new variable called '_NN'
    for index, item in enumerate(list_out):
	ANN_label=str(item+"_NN")
	ANN_label_all=str(item+"_NN_all")
	#Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
	construct_label=str(item+"_Con")
	temp_concat[construct_label]=np.nan
	#add a column for the constructed QC flag
	#This will be 1 if valid data from the tower else 99
	construct_flag_label=str(item+"_Con_QCFlag")
	temp_concat[construct_flag_label]=np.nan
	#Start by filling with valid  values from tower
	temp_concat[construct_flag_label][temp_concat[item].notnull()]=1
	temp_concat[construct_label][temp_concat[item].notnull()]=temp_concat[item]
	
	temp_concat[construct_flag_label][temp_concat[construct_flag_label].isnull()]=99
	temp_concat[construct_label][temp_concat[construct_flag_label]==99]=temp_concat[ANN_label]
	temp_concat[construct_label][temp_concat[construct_flag_label]==99]=temp_concat[ANN_label]
     
	
    #Clean up columns
    for index, item in enumerate(list_out):
	ANN_label_all=str(item+"_NN_all")
	del temp_concat[ANN_label_all]
    
    #Now if Fc storage is needed then write the temporary variable back to the Fc_Con
    if (Use_Fc_Storage=='Yes') and ('Fc_inc_store'  in list_out):
	New_combined['Fc_NN']=temp_concat['Fc_inc_store_NN']
	New_combined['Fc_Con']=temp_concat['Fc_inc_store_Con']
	New_combined['Fc_Con_QCFlag']=temp_concat['Fc_inc_store_Con_QCFlag']
	temp_concat=New_combined
	
    print"Completely Finished ANN outputs of "+str(list_out)+" at "+ Site_ID           
    
    return (temp_concat)




