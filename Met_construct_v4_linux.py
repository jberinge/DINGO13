############################################################################
# This script correlates tower data with external meteorology for a variable.  Then calulates the correlation coefficient.  Then adjusts new time series and gap fills 
#
# Programmed by Jason (Dec 1, 2012) to June 2013
############################################################################
import pandas as pd
from pandas.tools.plotting import scatter_matrix
import os
import datetime as dt
import xlrd
import string
import numpy as np
import netCDF4
import math
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['legend.fancybox'] = True
from matplotlib.backends.backend_pdf import PdfPages
from pylab import * 
import plotstuff

from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
import reportlab.pdfgen 
from reportlab.lib.units import inch

#Import OxFlux python modules
import meteorologicalfunctions as mf
#Import custom code for processing
import Gap_Fill_climatology_v2 as gap_fill_climatology

#import AWS_Flux_Correlate_plot_V1
from scipy import stats

# Functions here
#====================================================================

def excel_to_pydate(exceldate):
    datemode=0           # datemode: 0 for 1900-based, 1 for 1904-based
    pyear, pmonth, pday, phour, pminute, psecond = xlrd.xldate_as_tuple(exceldate, datemode)
    py_date = dt.datetime(pyear, pmonth, pday, phour, pminute, psecond)
    return(py_date)


def is_number(s):
    try:
        float(s)
        s=float(s)
        return s
    except Exception:
        s=np.nan
        return s
    else:
        s=np.nan
        return s

def is_nan(obj):
    test_number=math.isnan(obj)
    return test_number

def regress_func(x,xlabel,ylabel,newcolumn):
    #Do the regression.  Start by subsetting the two columns required.
    #Then drop any NaN case wise
    #reset (get rid of the index) and drop the index rather than keeping it as a column
    #so it can be passed to linregress
    xnow=x[[xlabel,ylabel]]
    xnow=xnow.dropna(how='any')
    xdata=xnow[xlabel].dropna().reset_index(drop=True)
    ydata=xnow[ylabel].dropna().reset_index(drop=True)   

   #Check to see if minumum number of samples before writing the results.
    if xnow[xlabel].count()>10 and xnow[ylabel].count()>10:
	slope, inter, rsqu, pval, se= stats.linregress(xdata,ydata)
	#print "stats:",slope, inter, rsqu, pval, se
	#Here use the original column to do apply the lin regresssion as 
	#values had been dropped previously
	x[newcolumn]=slope*x[xlabel]+inter
    else:
	x[newcolumn]=np.NaN
     #Return the new column to then join to existing file
    return x[newcolumn]

	
def construct_data(ALL_combined, VarToCorrelate, AWSVarToCorrelate, bestAWS_ID, Site_ID, corr_freq, myBaseforResults):
    #================================================================
    #                Main code started here
    #================================================================
       
    #If L3 then change Ta label so the code can run
    #ALL_combined = ALL_combined.rename(columns={'Ta': 'Ta_EC'})
    print "Correlation frequency ", corr_freq
    #Check for place to put resuolts - does it exist? If not create
    if not os.path.isdir(myBaseforResults):
        os.mkdir(myBaseforResults)
    #Then subdirectories
    if not os.path.isdir(myBaseforResults+"/AWS"):
        os.mkdir(myBaseforResults+"/AWS")
    if not os.path.isdir(myBaseforResults+"/AWS/"+VarToCorrelate):
	os.mkdir(myBaseforResults+"/AWS/"+VarToCorrelate)
    mypathforResults=myBaseforResults+"/AWS/"+VarToCorrelate

    #Do correlation and plots for AWS variables
    #-------------------------------------------
    ID1=bestAWS_ID[0]
    ID2=bestAWS_ID[1]
    ID3=bestAWS_ID[2]
    #Initialise the list empty with 4 items
    Labels=[None]*4
    Labels[0]=VarToCorrelate           #Flux tower label
    Labels[1]=AWSVarToCorrelate+"_"+ID1    #AWS first tower label
    Labels[2]=AWSVarToCorrelate+"_"+ID2    #AWS second tower label
    Labels[3]=AWSVarToCorrelate+"_"+ID3    #AWS third tower label
    
            
    #AWS_Flux_Correlate_plot_V1.AWS_correlate(ALL_combined,corr_freq,varname,bestAWS_ID)
    
    ## Temp statement to read in DAta that we have saved already.  Otherwise we process the whole thing again
    #ALL_combined= pd.DataFrame.load('ALL_combined.df')


    #Go through the list of VARIABLES of interest and ID's
    #and make sure each value in each Column is a number so it can be converted
    # to a float.  Then convert to a float so that it can be applied in stats.linregress
    #Problem is that the series imported into the data frame are obj not float and cant be used in array operations
    #These are the variables from teh AWS file that we want to import at some stage
    AWSvariables=['Ta','RH','Tw','DP','WS','Rain','GUST','P']
    #Loop through the AWS ID.  There are 3 AWS files in the ALL_combined dataframe
    for index1,ID in enumerate(bestAWS_ID):
	#Loop through the variables list
	for index2,AWSvariable in enumerate(AWSvariables):
	    tempAWS_label=AWSvariable+'_'+ID
	    newAWS_label=AWSvariable+'_'+ID+'_NEW'
	    ALL_combined[tempAWS_label]=ALL_combined[tempAWS_label].map(is_number)
	    ALL_combined[tempAWS_label]=ALL_combined[tempAWS_label].astype('float64')
	    
    #Do for the same for current variables in the Flux data (i.e. Ta, Ah, etc).
    ALL_combined[VarToCorrelate]=ALL_combined[VarToCorrelate].map(is_number)
    ALL_combined[VarToCorrelate]=ALL_combined[VarToCorrelate].astype('float64')    
    
    # Do calculations on met variables depending on variable type
    if VarToCorrelate=='Ah':
	for index,ID in enumerate(bestAWS_ID):
	    Ta_temp_label='Ta_'+ID
	    RH_temp_label='RH_'+ID	    
	    #Here call function from Ozflux meteorolgoical functions script
	    ALL_combined[Labels[index+1]]=mf.absolutehumidityfromRH(ALL_combined[Ta_temp_label],ALL_combined[RH_temp_label] )				    
	    #Replace -9999 in returned data to nan's.  Leave gaps and do regressions.  Then fill
	    ALL_combined[Labels[index+1]].replace(-9999,value=np.nan,method='None', inplace=True)

    if VarToCorrelate=='Ws':
	for index,ID in enumerate(bestAWS_ID):	    
	    #Here convert Bom km/h to flux m/s
	    ALL_combined[Labels[index+1]]=ALL_combined[Labels[index+1]]*0.27778				    
	    #Replace -9999 in returned data to nan's.  Leave gaps and do regressions.  Then fill
	    ALL_combined[Labels[index+1]].replace(-9999,value=np.nan,method='None', inplace=True)    

    if VarToCorrelate=='P':
	for index,ID in enumerate(bestAWS_ID):	    
	    #Here convert Bom hPa to kPa as per flux
	    ALL_combined[Labels[index+1]]=ALL_combined[Labels[index+1]]/10				    
	    #Replace -9999 in returned data to nan's.  Leave gaps and do regressions.  Then fill
	    ALL_combined[Labels[index+1]].replace(-9999,value=np.nan,method='None', inplace=True)    
    
    #subset the large dataframe to variables we want for this process.  One at a time	
    SubsetDF=ALL_combined[[Labels[0],Labels[1],Labels[2],Labels[3]]]
    
                    
    #Get things to report. Number of samples in each variable. Number of Nan's
    #This Gives the SIZE = total records in the by GROUP
    #The Count gives non Nan values by Group
    stats_size=SubsetDF.groupby([lambda x: x.year,lambda x: x.month]).size()
    stats_count=SubsetDF.groupby([lambda x: x.year,lambda x: x.month]).count()

    #Do  FIRST stats correlation
    #Initialise the list empty with 3 items
    Statsresults=[]   
    SubsetDF1=SubsetDF[[Labels[0],Labels[1]]]
    
    SubsetDF1_withNans = SubsetDF1
    SubsetDF1 = SubsetDF1.dropna(how='any')
    #Check to see what basis we want the correlation to be performed. All, annual or monthly
        #Use SCiPy stats scipy.stats.linregress(x, y)
    if corr_freq=='monthly':
        SubsetDF1_grouped=SubsetDF1.groupby([lambda x: x.year,lambda x: x.month])
        SubsetDF1_withNans_grouped=SubsetDF1_withNans.groupby([lambda x: x.year,lambda x: x.month])
	Corr_grouped1=SubsetDF1_grouped.apply(lambda x: stats.linregress(x [Labels[1]], x [Labels[0]]))
    elif corr_freq=='annual':
        SubsetDF1_grouped=SubsetDF1.groupby([lambda x: x.year])
        SubsetDF1_withNans_grouped=SubsetDF1_withNans.groupby([lambda x: x.year])
	Corr_grouped1=SubsetDF1_grouped.apply(lambda x: stats.linregress(x [Labels[1]], x [Labels[0]]))
    else:
        SubsetDF1_grouped=SubsetDF1
        SubsetDF1_withNans_grouped=SubsetDF1_withNans
	Corr_grouped1=(stats.linregress(SubsetDF1[Labels[1]], SubsetDF1[Labels[0]]))
    
    #Use SCiPy stats scipy.stats.linregress(x, y)
    Corr_all1=stats.linregress(SubsetDF1[Labels[1]], SubsetDF1[Labels[0]])
    
    #Do  SECOND stats correlation
    #Initialise the list empty with 3 items 
    SubsetDF2=SubsetDF[[Labels[0],Labels[2]]]
    SubsetDF2 = SubsetDF2.dropna(how='any')
    #Check to see what basis we want the correlation to be performed. All, annual or monthly
    if corr_freq=='monthly':
        SubsetDF2_grouped=SubsetDF2.groupby([lambda x: x.year,lambda x: x.month])
	Corr_grouped2=SubsetDF2_grouped.apply(lambda x: stats.linregress(x [Labels[2]], x [Labels[0]]))
    elif corr_freq=='annual':
        SubsetDF2_grouped=SubsetDF2.groupby([lambda x: x.year])
	Corr_grouped2=SubsetDF2_grouped.apply(lambda x: stats.linregress(x [Labels[2]], x [Labels[0]]))
    else:
        SubsetDF2_grouped=SubsetDF2
	Corr_grouped2=stats.linregress(SubsetDF2[Labels[2]], SubsetDF2[Labels[0]])
    
    Corr_all2=stats.linregress(SubsetDF2[Labels[2]], SubsetDF2[Labels[0]])

    
    #Do  THIRD stats correlation
    #Initialise the list empty with 3 items
    SubsetDF3=SubsetDF[[Labels[0],Labels[3]]]
    SubsetDF3 = SubsetDF3.dropna(how='any')
    #Check to see what basis we want the correlation to be performed. All, annual or monthly
    if corr_freq=='monthly':
        SubsetDF3_grouped=SubsetDF3.groupby([lambda x: x.year,lambda x: x.month])
	Corr_grouped3=SubsetDF3_grouped.apply(lambda x: stats.linregress(x [Labels[3]], x [Labels[0]]))
    elif corr_freq=='annual':
        SubsetDF3_grouped=SubsetDF3.groupby([lambda x: x.year])
	Corr_grouped3=SubsetDF3_grouped.apply(lambda x: stats.linregress(x [Labels[3]], x [Labels[0]]))
    else:
        SubsetDF3_grouped=SubsetDF3
	Corr_grouped3=stats.linregress(SubsetDF3[Labels[3]], SubsetDF3[Labels[0]])
    Corr_all3=stats.linregress(SubsetDF3[Labels[3]], SubsetDF3[Labels[0]])
  
    
    #Output Panda DF to list so it can be used to create a tabel in ReportLab
    #Create list for EACH AWS Site ID and output BY gruop for YEAR and MONTH
    #totable2=Corr_grouped1.unstack
    #totable3=Corr_grouped2.unstack
    #totable4=Corr_grouped3.unstack 

    
    #Create a Pandas DF of the ALL data stats.  A linregression of all data. 
    #Returns'slope', 'intercept', 'r_value', 'p_value', 'std_err'
    #That is then put in a Pandas DF and output to list
    d1=[Corr_all1,Corr_all2, Corr_all3]
    Corr_data=pd.DataFrame(d1,columns=['slope', 'intercept', 'r_value', 'p_value', 'std_err'])
    #Create a DF with just the site ID's with column names ID
    d2=pd.Series([ID1, ID2, ID3],name='ID')
    #Join the two DF together
    d3=Corr_data.join(d2)
    #Set the DF Index to ID
    new= d3.set_index('ID')
    #Sort the DF by best (highest r quared value .  This finds the order in which to use sites in the data gapfill.
    new=new.sort(columns='r_value',ascending=False)
    #Output Panda DF to list so it can be used to create a tabel in ReportLab
    #Create list of the Site ID correlations for ALL data.
    totable1=list(new.T.itertuples())
    
    #Outputs (Tables and PDF)
    #===============================
    #Output PDF Canvas
    doc1name=(mypathforResults+'/Table Tower and AWS correlation_%s_%s.pdf' % (VarToCorrelate,Site_ID))
    doc1 = SimpleDocTemplate(doc1name,pagesize=A4, rightMargin=40,leftMargin=40, topMargin=20,bottomMargin=20)
    styles=getSampleStyleSheet()  
    width, height = A4       
    # container for the 'Flowable' objects
    #Then put objects in there to fill and generate PDF later
    container1 = []
    #Create first table of best ALL data correlations
    table1=Table(totable1)
    ptext1 = (' <font size=12>The data ALL data correlation stats for variable %s at AWS site ID %s</font>' % (VarToCorrelate, Site_ID)) 
    ptext1a = (' <font size=12>Sites are         %s                     %s                 %s</font>' % (ID1,ID2,ID3))   
    container1.append(Paragraph(ptext1, styles["Normal"]))
    container1.append(Paragraph(ptext1a, styles["Normal"]))
    container1.append(Spacer(1, 12))
    #table1.setStyle(TableStyle([('BACKGROUND',(0,0),(0,rowsto_output1),colors.green),('BACKGROUND',(1,0),(1,rowsto_output1),colors.yellow)]))
    #set_column_titles(array('ID','Dist', 'Have Data','Lat','Long','Name','Start', 'End','Elev'))    
    container1.append(table1)
    container1.append(Spacer(1, 20))
    doc1.build(container1)
    

    #NOTE#
    ######
    #I couldnt get anything appropriate to pas to the table constructor here.  So I have commented it out
    #It works for a DF converted to LIST above for small table
    #Problem here is we want to pass the DF grouped object.  This has to to_list method
    
    ##Output PDF Canvas 2
    #doc2name=("Table Tower and AWS correlationwith BoM ID %s variable %s at %s.pdf" % (ID1,VarToCorrelate,Site_ID))
    #doc2 = SimpleDocTemplate(doc2name,pagesize=A4, rightMargin=40,leftMargin=40, topMargin=20,bottomMargin=20)
    #styles=getSampleStyleSheet()  
    #width, height = A4       
    ## container for the 'Flowable' objects
    ##Then put objects in there to fill and generate PDF later
    #container2 = []
    #ptext2 = (' <font size=12>The data shows a breakdown by YEAR and MONTH of the correlation stats for variable %s at AWS site ID %s</font>' % (VarToCorrelate, Site_ID))     
    #table2=Table(totable2)
    #container2.append(Paragraph(ptext2, styles["Normal"])) 
    #container2.append(table2) 
    #container2.append(Spacer(1, 20))
    #doc2.build(container2)
    
    ##Output PDF Canvas 3
    #doc3name=("Table Tower and AWS correlationwith BoM ID %s variable %s at %s.pdf" % (ID2,VarToCorrelate,Site_ID))
    #doc3 = SimpleDocTemplate(doc3name,pagesize=A4, rightMargin=40,leftMargin=40, topMargin=20,bottomMargin=20)
    #styles=getSampleStyleSheet()  
    #width, height = A4       
    ## container for the 'Flowable' objects
    ##Then put objects in there to fill and generate PDF later
    #container3 = []    
    #ptext3 = (' <font size=12>The data shows a breakdown by YEAR and MONTH of the correlation stats for variable %s at AWS site ID %s</font>' % (VarToCorrelate, Site_ID))     
    #table3=Table(totable3)
    #container3.append(Paragraph(ptext3, styles["Normal"]))
    #container3.append(table3)
    #container3.append(Spacer(1, 20))
    #doc3.build(container3)   
    
    ##Output PDF Canvas 4
    #doc4name=("Table Tower and AWS correlationwith BoM ID %s variable %s at %s.pdf" % (ID3,VarToCorrelate,Site_ID))
    #doc4 = SimpleDocTemplate(doc4name,pagesize=A4, rightMargin=40,leftMargin=40, topMargin=20,bottomMargin=20)
    #styles=getSampleStyleSheet()  
    #width, height = A4       
    ## container for the 'Flowable' objects
    ##Then put objects in there to fill and generate PDF later
    #container4 = []      
    #ptext4 = (' <font size=12>The data shows a breakdown by YEAR and MONTH of the correlation stats for variable %s at AWS site ID %s</font>' % (VarToCorrelate, Site_ID))     
    #table4=Table(totable4)
    #container4.append(Paragraph(ptext4, styles["Normal"])) 
    #container4.append(table4)     
    #container4.append(Spacer(1, 20))
    #doc4.build(container4)      
    
    # Write output files
    ################################################################
    #when write file at this stage when ALL then the output is just a list and cant be written using the 
    #Panadas file output
    if corr_freq != 'all':
	outputfilename1=str((mypathforResults+'/'+'Subset grouped AWS ID %s variable %s site %s.csv' % (ID1, VarToCorrelate, Site_ID)))
	Corr_grouped1.to_csv(outputfilename1, sep=',')  
	outputfilename2=str((mypathforResults+'/'+'Subset grouped AWS ID %s variable %s site %s.csv' % (ID2, VarToCorrelate, Site_ID)))
	Corr_grouped2.to_csv(outputfilename2, sep=',')  
	outputfilename3=str((mypathforResults+'/'+'Subset grouped AWS ID %s variable %s site %s.csv' % (ID3, VarToCorrelate, Site_ID)))
	Corr_grouped3.to_csv(outputfilename3, sep=',') 
	outputfilename4=str((mypathforResults+'/'+'Correlation stats all data  variable %s site %s.csv' % ( VarToCorrelate, Site_ID)))
	new.to_csv(outputfilename4, sep=',')     
    
    # Produce graphs
    ################################################################
    # Produce ALL data graph first
    # The DF 'new' has info in this format ID as index then 'slope', 'intercept', 'r_value', 'p_value', 'std_err'
    #plotdata = pd.DataFrame(SubsetDF)
    #Test plot using Pandas plots tools
    #Produce a 4 way correlation plot.
    #fig1 = plt.figure()
    #scatter_matrix(SubsetDF, alpha=0.2, figsize=(8, 8), diagonal='kde')
    #plt.show()
    
    #Calculate some things
    n_datapoints=len(SubsetDF)
    startdate= SubsetDF.index[0]
    enddate= SubsetDF.index[n_datapoints-1]

   
    
    #fig2 = plt.figure()
    #plt.scatter(SubsetDF[Labels[0]], SubsetDF[Labels[1]],'yo')
    ##plt.scatter(SubsetDF[Labels[0]], SubsetDF[Labels[2]],'g+')
    ##plt.scatter(SubsetDF[Labels[0]], SubsetDF[Labels[3]],'r^')
    #plt.show()
    
    #PLOT Three dta sets on Here
    #Loop through each ID in the list withing the DF 'new'
    #Remember that this list has been sorted so top has highest correlation
    slopetemp=[None]*3
    intercepttemp=[None]*3
    tempx_line0=[]
    tempy_line0=[]
    tempx_line1=[]
    tempy_line1=[]
    tempx_line2=[]
    tempy_line2=[]

    for index, IDx in enumerate(new.index):
	#set some temporary varibales to pass plot variables to
	plotXcolumn=AWSVarToCorrelate+"_"+IDx
	plotYcolumn=VarToCorrelate
    
	#For plotting we want to find tha range of varaible values across all sites and tower to get entire range to plot
	DFmins= SubsetDF.min()
	DFmaxs= SubsetDF.max()
	scale_min= int(min(DFmins[Labels[0]],DFmins[Labels[1]],DFmins[Labels[2]],DFmins[Labels[3]]))-1
	scale_max= int(max(DFmaxs[Labels[0]],DFmaxs[Labels[1]],DFmaxs[Labels[2]],DFmaxs[Labels[3]]))+1
	
	#create series to plot line
	slopetemp[index]=new.lookup([IDx], ['slope'])
	intercepttemp[index]=new.lookup([IDx], ['intercept'])

	for increment in range(scale_min,scale_max):
	    if index==0:
		tempx_line0.append(increment)
		tempy_line0.append(slopetemp[index]*increment+intercepttemp[index])
	    elif index==1:
		tempx_line1.append(increment)
		tempy_line1.append(slopetemp[index]*increment+intercepttemp[index])
	    elif index==2:
		tempx_line2.append(increment)
		tempy_line2.append(slopetemp[index]*increment+intercepttemp[index])	

	#Produce these plots in reverse order 2 to 0 so that plot with least correlations is on the bottom
	if index==2:
	    plt.plot(SubsetDF[plotXcolumn], SubsetDF[plotYcolumn], 'b+',tempx_line2, tempy_line2, ':b' ,label=IDx,linewidth=2) 
	elif index==1:
	    plt.plot(SubsetDF[plotXcolumn], SubsetDF[plotYcolumn], 'ro' ,tempx_line1, tempy_line1, '--r',label=IDx,linewidth=2) 
	elif index==0:
	    plt.plot(SubsetDF[plotXcolumn], SubsetDF[plotYcolumn], 'y>',tempx_line0, tempy_line0, '-y',label=IDx,linewidth=2) 
		
    plt.xlim(scale_min, scale_max)
    plt.ylim(scale_min, scale_max)  
    #create text for ID and r2 box
    graphtext1=('Station ID      R2 \n' +str(new.index[0]) +"     "+ "{0:.3f}".format(float(new.lookup([ID1], ['r_value'])))+"\n"
	+str(new.index[1]) +"     "+ "{0:.3f}".format(float(new.lookup([ID2], ['r_value'])))+"\n"
	+str(new.index[2]) +"     "+ "{0:.3f}".format(float(new.lookup([ID3], ['r_value']))))
    #create text for start and end dates
    graphtext2=('Data start date: '+str(startdate)+'\n'
                +'End date: '+str(enddate)+'\n'
                +'Number records: '+str(n_datapoints))   
    units=' oC'
    plt.figtext(0.5,0.25,graphtext1, bbox=dict())
    plt.figtext(0.5,0.1,graphtext2, bbox=dict())
    plt.title('Variable  '+VarToCorrelate+ ' for ALL data' + ' at ' +Site_ID)
    plt.xlabel('BoM AWS stations ' + '('+units+')')
    plt.ylabel(Labels[0]+ '   ' + '('+units+')')
    plt.legend(shadow=True, fancybox=True,loc='best')
    
    #Output to PDF using PdfPages a backend for MatPlotLib
    fname_graph=mypathforResults+'/'+'Linear Plot Tower vs BoM AWSs for ALL data - Variable  '+VarToCorrelate+ ' at ' +Site_ID+'.pdf'
    # Create the PdfPages object to which we will save the pages:
    pdf = PdfPages(fname_graph)
    savefig(pdf, format='pdf',facecolor='w', edgecolor='w') # note the format='pdf' argument!
    #show()	
    close()		
    pdf.close()    
    
    #Now try and produce a panel of graph of 6 
    ##############################################
    #This just finds the number of Year and Month groups so that we can use it later for plotting
    #Where we have 6 graphs to a page.  How many graphs do we need?
    subsets=len(SubsetDF1_grouped)
    
    if corr_freq=='monthly':
	for (k1, k2), group in SubsetDF1_grouped:
	    #Calculate some things
	    n_datapoints=len(group)
	    startdate= group.index[0]
	    enddate= group.index[n_datapoints-1]
	    #Set temp variables and lists
	    slopetemp=[]
	    intercepttemp=[]
	    tempx_line=[]
	    tempy_line=[]
    
	    #set some temporary varibales to pass plot variables to
	    plotXcolumn=VarToCorrelate+"_"+IDx
	    plotYcolumn=VarToCorrelate
	
	    #For plotting we want to find the range of varaible values across all sites and tower to get entire range to plot
	    DFmins= group.min()
	    DFmaxs= group.max()
	    scale_min= int(min(DFmins[Labels[0]],DFmins[Labels[1]]))-1
	    scale_max= int(max(DFmaxs[Labels[0]],DFmaxs[Labels[1]]))+1
	    
	    #create series to plot line
	    #Need to extract the linear regression stats done bygroup earlier
    
	    slopetemp, intercepttemp, r_valuetemp, p_valuetemp, std_errtemp = stats.linregress(group[Labels[1]],group[Labels[0]])
	    
	    for increment in range(scale_min,scale_max):
		    tempx_line.append(increment)
		    tempy_line.append(slopetemp*increment+intercepttemp)
		    
	    ## Could work for later  pd.merge(df, k1_means, left_on='key1', right_index=True	
	    #Produce the plot 
	    plt.plot(group[Labels[1]], group[Labels[0]], 'go',tempx_line, tempy_line, ':b' ,label=IDx,linewidth=2) 
	    #Set the scale mins and maxs
	    plt.xlim(scale_min, scale_max)
	    plt.ylim(scale_min, scale_max)  
	    #create text for ID and r2 box
	    graphtext1=     str('intercept   ' + str("{0:.2f}".format(intercepttemp) +'\n')
		              + 'slope       ' + str("{0:.2f}".format(slopetemp)) +'\n'
	                      + 'r value     ' + str("{0:.2f}".format(r_valuetemp)) +'\n'
	                      + 'p_value     ' + str("{0:.2f}".format(p_valuetemp)) +'\n'
		              + 'std_err     ' + str("{0:.2f}".format(std_errtemp)) +'\n')  
	    #create text for start and end dates
	    graphtext2=('Data start date: '+str(startdate)+'\n'
		        +'End date: '+str(enddate)+'\n'
		        +'Number records: '+str(n_datapoints))   
	    units=' oC'
	    plt.figtext(0.7,0.3,graphtext1, bbox=dict())
	    plt.figtext(0.5,0.13,graphtext2, bbox=dict())
	    plt.title('Tower vs Best BoM AWS  Year '+ str(k1) +'  Month '+ str(k2) 
	              +' Variable  '+VarToCorrelate+ ' at ' +Site_ID
	              +'\n'+FLUXfilename)
	    plt.xlabel('BoM AWS station ' + '('+units+')')
	    plt.ylabel(Labels[0]+ '   ' + '('+units+')')
	    plt.legend(shadow=True, fancybox=True,loc='best')
	    
	    #Output to PDF using PdfPages a backend for MatPlotLib
	    fname_graph=mypathforResults+'/'+'Linear Plot Tower vs BoM AWSs for Year '+ str(k1) +' and Month '+ str(k2) +' - Variable  '+VarToCorrelate+ ' at ' +Site_ID+'.pdf'
	    # Create the PdfPages object to which we will save the pages:
	    pdf = PdfPages(fname_graph)
	    savefig(pdf, format='pdf',facecolor='w', edgecolor='w') # note the format='pdf' argument!
	    #show()	
	    close()		
	    pdf.close() 	
    elif corr_freq=='annual':
	for k1 , group in SubsetDF1.groupby([lambda x: x.year]):
	    
	    #Calculate some things
	    n_datapoints=len(group)
	    startdate= group.index[0]
	    enddate= group.index[n_datapoints-1]
	    #Set temp variables and lists
	    slopetemp=[]
	    intercepttemp=[]
	    tempx_line=[]
	    tempy_line=[]
    
	    #set some temporary varibales to pass plot variables to
	    plotXcolumn=VarToCorrelate+"_"+IDx
	    plotYcolumn=VarToCorrelate
	
	    #For plotting we want to find the range of varaible values across all sites and tower to get entire range to plot
	    DFmins= group.min()
	    DFmaxs= group.max()
	    scale_min= int(min(DFmins[Labels[0]],DFmins[Labels[1]]))-1
	    scale_max= int(max(DFmaxs[Labels[0]],DFmaxs[Labels[1]]))+1
	    
	    #create series to plot line
	    #Need to extract the linear regression stats done bygroup earlier
    
	    slopetemp, intercepttemp, r_valuetemp, p_valuetemp, std_errtemp = stats.linregress(group[Labels[1]],group[Labels[0]])
	    
	    for increment in range(scale_min,scale_max):
		    tempx_line.append(increment)
		    tempy_line.append(slopetemp*increment+intercepttemp)
	    ## Could work for later  pd.merge(df, k1_means, left_on='key1', right_index=True	
	    #Produce the plot 
	    plt.plot(group[Labels[1]], group[Labels[0]], 'go',tempx_line, tempy_line, ':b' ,label=IDx,linewidth=2) 
	    #Set the scale mins and maxs
	    plt.xlim(scale_min, scale_max)
	    plt.ylim(scale_min, scale_max)  
	    #create text for ID and r2 box
	    graphtext1=str('intercept  ' + str("{0:.2f}".format(intercepttemp) +'\n')
	                      + 'r value      ' + str("{0:.2f}".format(r_valuetemp)) +'\n'+'p_value      ' + str("{0:.2f}".format(p_valuetemp)) +'\n'
	                      + 'std_err      ' + str("{0:.2f}".format(std_errtemp)) +'\n')  
	    #create text for start and end dates
	    graphtext2=('Data start date: '+str(startdate)+'\n'
	                +'End date: '+str(enddate)+'\n'
	                +'Number records: '+str(n_datapoints))   
	    units=' oC'
	    plt.figtext(0.7,0.3,graphtext1, bbox=dict())
	    plt.figtext(0.5,0.13,graphtext2, bbox=dict())
	    plt.title('Tower vs Best BoM AWS  Year '+ str(k1) +' Variable  '+VarToCorrelate+ ' at ' +Site_ID)
	    plt.xlabel('BoM AWS station ' + '('+units+')')
	    plt.ylabel(Labels[0]+ '   ' + '('+units+')')
	    plt.legend(shadow=True, fancybox=True,loc='best')
	    
	    #Output to PDF using PdfPages a backend for MatPlotLib
	    fname_graph=mypathforResults+'/'+'Linear Plot Tower vs BoM AWSs for Year '+ str(k1) +' - Variable  '+VarToCorrelate+ ' at ' +Site_ID+'.pdf'
	    # Create the PdfPages object to which we will save the pages:
	    pdf = PdfPages(fname_graph)
	    savefig(pdf, format='pdf',facecolor='w', edgecolor='w') # note the format='pdf' argument!
	    #show()	
	    close()		
	    pdf.close() 	
	
    ###########################################################
    # Do the correlation analysis and APPLY it to the columns 
    ##########################################################
    
    #Apply depending on the frequency required
    # 'slope', 'inter', 'rsqu', 'pval', 'se'
    # suffix with _all, _yr or _mon to indicate
    #Suffux with site ID i.e. _045623
    #STEP 1. Create Pandas DF for each of the Yearly and Monthly breakdowns to get stats
    #But statslinreg wont work with Nans so start with DF without Nans.
    #Then later step 2 come back to fill the Nans in the original DF = SubSetDF
    
    #Create a new label for pandas df column for the contructed variable (and the QC flag) and column to fill label
    construct_label=str(VarToCorrelate+"_Con")
    fill_label=str(VarToCorrelate)
    #start by copying the existing data from the tower to the construct column, then fill the missing bits
    SubsetDF[construct_label]=SubsetDF[fill_label]  

    #Also later we need the following columns alrady defined in the DF so here goes
    for ID in new.index:
	corr_label_All=str(VarToCorrelate+"_AllCorr"+"_"+ID)
	SubsetDF[corr_label_All]=zeros
	corr_label_Ann=str(VarToCorrelate+"_AnnCorr"+"_"+ID)
	SubsetDF[corr_label_Ann]=zeros
	corr_label_Mon=str(VarToCorrelate+"_MonCorr"+"_"+ID)
	SubsetDF[corr_label_Mon]=zeros	

    #Do ALL correlation
    slope_all={}
    inter_all={}
    rsqu_all ={}
    pval_all ={}
    se_all={}
    
    for ID in new.index:
	xlabel=str(AWSVarToCorrelate+"_"+ID)
	ylabel=str(VarToCorrelate)
	
	temp=SubsetDF.dropna(how='any')
	xvalues=temp[xlabel]
	yvalues=temp[ylabel]

	slope_all[ID], inter_all[ID], rsqu_all[ID], pval_all[ID], se_all[ID]= stats.linregress(xvalues,yvalues)
	print "ID : ",ID, slope_all[ID]
    
    for key, value in slope_all.iteritems():
	print key, value
    
    
    #Do ANNUAL correlation
    #prepare the datasets
    #Setup DF grouped by year
    temp_a1=SubsetDF
    temp_a2=SubsetDF.dropna(how='any')
    tempAnnualgrouped=temp_a1.groupby([lambda x: x.year])
    tempAnnualgrouped_noNans=temp_a2.groupby([lambda x: x.year]) 

    #loop through the BoM sites Best to worst correlation over 3 sites
    #setup for a new Pandas datatable to put results in
    #Create a list fits
    AnnualStatsList=[]
    for ID in new.index:
	xlabel=str(AWSVarToCorrelate+"_"+ID)
	ylabel=str(VarToCorrelate)

	for a1 , group in tempAnnualgrouped_noNans:
	    xvalues=group[xlabel]
	    yvalues=group[ylabel]
	    
	    slope_yr, inter_yr, rsqu_yr, pval_yr, se_yr= stats.linregress(xvalues,yvalues)
	    #output to list and later make into DF
	    AnnualStatsList.append([a1,ID, slope_yr, inter_yr, rsqu_yr, pval_yr, se_yr])
	    
    AnnualStats=pd.DataFrame(AnnualStatsList, columns=['year','ID', 'slope', 'inter', 'rsqu', 'pval', 'se'])
    #AnnualStats=temp1.set_index(['year','ID'])
    
    print AnnualStats.head(5)
	
    #Do MONTHLY correlation
    #prepare the datasets
    #Setup DF grouped by year
    temp_m1=SubsetDF
    temp_m2=SubsetDF.dropna(how='any')
    
    tempMonthlyGrouped_noNans = temp_m2.groupby([lambda x: x.year,lambda x: x.month])
    #loop through the BoM sites Best to worst correlation over 3 sites
    #setup for a new Pandas datatable to put results in
    #Create a list fits
    MonthlyStatsList=[]
    for ID in new.index:
	xlabel=str(AWSVarToCorrelate+"_"+ID)
	ylabel=str(VarToCorrelate)

	for (m1 , m2) , group in tempMonthlyGrouped_noNans:
	    xvalues=group[xlabel]
	    yvalues=group[ylabel]
	    
	    slope_yr, inter_yr, rsqu_yr, pval_yr, se_yr= stats.linregress(xvalues,yvalues)
	    #output to list and later make into DF
	    MonthlyStatsList.append([m1, m2 ,ID, slope_yr, inter_yr, rsqu_yr, pval_yr, se_yr])
	    
    MonthlyStats=pd.DataFrame(MonthlyStatsList, columns=['year','month','ID', 'slope', 'inter', 'rsqu', 'pval', 'se'])
    #MonthlyStats=temp2.set_index(['year','month','ID'])
    print MonthlyStats.head(5)	
	
    #Now step2 Apply the stats back to original dataframe
    #=====================================================
    #Create the series first
  
    #Do MONTHLY series  
    for ID in new.index:
	#print "Applying correlations for Monthly" + ID
	#for ID in new.index:
	xlabel=str(AWSVarToCorrelate+"_"+ID)
	ylabel=str(VarToCorrelate)
	#Create new variable
	corr_label=str(VarToCorrelate+"_MonCorr"+"_"+ID)
	temp_m1[corr_label].iloc[:]=np.nan
	testshape=len(temp_m1.groupby([lambda x: x.year,lambda x: x.month]))
	MonCorResults = temp_m1.groupby([lambda x: x.year,lambda x: x.month], group_keys=False, as_index=False).apply(regress_func,xlabel,ylabel,corr_label)
	#Do a fill of data where one value is missing.  This solves the problem of when the AWS data is 
	#60 minutes but other AWS data is 30 minutes.  Under that case the selected AWS ID oscillates 
	#between station IDs
	MonCorResults.fillna(method='ffill', limit=2, inplace=True)
	temp_m1[corr_label]=MonCorResults



    #Do ALL series  
    for ID in new.index:
	#print "Applying correlations for ALL" + ID
	#for ID in new.index:
	xlabel=str(AWSVarToCorrelate+"_"+ID)
	ylabel=str(VarToCorrelate)
	#Create new variable
	corr_label=str(VarToCorrelate+"_AllCorr"+"_"+ID)
	temp_m1[corr_label].iloc[:]=np.nan
	#AllCorResults = temp_m1.apply(regress_func, axis=0, broadcast=False, raw=False, args=(xlabel,ylabel,corr_label))
	#Do the regression.  Start by subsetting the two columns required.
	#Then drop any NaN case wise
	#reset (get rid of the index) and drop the index rather than keeping it as a column
	#so it can be passed to linregress
	xnow=temp_m1[[xlabel,ylabel]]
	xnow=xnow.dropna(how='any')
	xdata=xnow[xlabel].dropna().reset_index(drop=True)
	ydata=xnow[ylabel].dropna().reset_index(drop=True)   
	slope, inter, rsqu, pval, se= stats.linregress(xdata,ydata)
	
	print "stats:",slope, inter, rsqu, pval, se
	#Here use the original column to do apply the lin regresssion as 
	#values had been dropped previously
	temp_m1[corr_label]=slope*temp_m1[xlabel]+inter  
    	#Do a fill of data where one value is missing.  This solves the problem of when the AWS data is 
	#60 minutes but other AWS data is 30 minutes.  Under that case the selected AWS ID oscillates 
	#between station IDs
	temp_m1[corr_label].fillna(method='ffill', limit=2, inplace=True)	
	
	
    #Do ANNUAL series  
    for ID in new.index:
	#print "Applying correlations for Annual" + ID
	#for ID in new.index:
	xlabel=str(AWSVarToCorrelate+"_"+ID)
	ylabel=str(VarToCorrelate)
	#Create new variable
	corr_label=str(VarToCorrelate+"_AnnCorr"+"_"+ID)
	corr_label_all=str(VarToCorrelate+"_AllCorr"+"_"+ID)
	temp_m1[corr_label].iloc[:]=np.nan
	testshape=len(temp_m1.groupby([lambda x: x.year]))
	
	#test to see if more than one year.  Otherwise this will give an error
	#Instead use the ALL series instead and just copy.
	#this is calculated in previous block
	#pdb.set_trace()
	if testshape >  1:
         AnnCorResults = temp_m1.groupby([lambda x: x.year]).apply(regress_func,xlabel,ylabel,corr_label)
         AnnCorResults.index=AnnCorResults.index.droplevel(0)
         temp_m1[corr_label]=AnnCorResults	
	else:
	    temp_m1[corr_label]=temp_m1[corr_label_all]
	#Do a fill of data where one value is missing.  This solves the problem of when the AWS data is 
	#60 minutes but other AWS data is 30 minutes.  Under that case the selected AWS ID oscillates 
	#between station IDs
	temp_m1[corr_label].fillna(method='ffill', limit=2, inplace=True)	


    #########################################################################################
    #STEP3 Now we have created a number of Time series that have we can use for the gap filling
    #Now apply the gap filling!
    ##########################################################################################
    #Apply this based on the frequency variable passed along All, Annual or Monthly   
    #create different columns for each frequency of gap filling
    #################################################################
    #add a column for the constructed QC flag
    #This will be 1 if valid data from the tower else the AWS ID or
    construct_flag_label=str(VarToCorrelate+"_Con_QCFlag")
    temp_m1[construct_flag_label]=np.nan 

   #Set Construct flag equal to 1 to say that EC data was used and QC is OK. Later add to this  the AWS ID based on missing data 
    temp_m1[construct_flag_label][((temp_m1[VarToCorrelate])>-50) & ((temp_m1[VarToCorrelate])<100)]=1
    
    if corr_freq=='all':
	for row in range(0,(temp_m1.shape[0])):
	    # Check to see if the Current row of variable to contrustct a variable is a number. 
	    # Call function and fill if necessary 
	    #Get the column number from the name
	    construct_col_number=temp_m1.columns.get_loc(construct_label)    
	    #Use iloc indexing to get the points in the dataframe		
	    if is_nan(temp_m1.iloc[row,construct_col_number])==True:
		#Loop through BoM ID's first
		#Fill in using each of the ID stations in the list	
		#This should fill using best station in the list.  If any Nans then
		#they should be filled by the next statio, and then the next
		#Assign the AWS station ID as the flagfor ID in new.index:
		#Define Flag variable name and column number
		corr_label_All_Flag=str(VarToCorrelate+"_Con_QCFlag")
		corr_Flag_col_number=temp_m1.columns.get_loc(corr_label_All_Flag)
		#Loop through again with Correlate ALL 
		for ID in new.index:		
		    corr_label_All=str(VarToCorrelate+"_AllCorr"+"_"+ID)
		    corr_All_col_number=temp_m1.columns.get_loc(corr_label_All)
		    if (is_nan(temp_m1.iloc[row,corr_All_col_number])==False) and (is_nan(temp_m1.iloc[row,construct_col_number])==True):		    
			temp_m1.iloc[row,construct_col_number]=temp_m1.iloc[row,corr_All_col_number]  
			temp_m1.iloc[row,corr_Flag_col_number]=999

	    
	print "All Stats"
	print "Mean ", temp_m1[construct_label].mean()
	print "Count ", temp_m1[construct_label].count()		

    elif corr_freq=='annual':
	for row in range(0,(temp_m1.shape[0])):
	    # Check to see if the Current row of variable to contrustct a variable is a number. 
	    # Call function and fill if necessary 
	    #Get the column number from the name
	    construct_col_number=temp_m1.columns.get_loc(construct_label)
	     
	    #Use iloc indexing to get the points in the dataframe
		
	    if is_nan(temp_m1.iloc[row,construct_col_number])==True:
		#Loop through BoM ID's first
		#Fill in using each of the ID stations in the list	
		#This should fill using best station in the list.  If any Nans then
		#they should be filled by the next statio, and then the next
		#Assign the AWS station ID as the flagfor ID in new.index:
		#Define Flag variable name and column number
		corr_label_Ann_Flag=str(VarToCorrelate+"_Con_QCFlag")
		corr_Flag_col_number=temp_m1.columns.get_loc(corr_label_Ann_Flag)
		
		#If Construct is missing AND BoM is present then fill and and flag value
		for ID in new.index:
		    corr_label_Ann=str(VarToCorrelate+"_AnnCorr"+"_"+str(ID))
		    corr_col_number=temp_m1.columns.get_loc(corr_label_Ann)		    
		    if (is_nan(temp_m1.iloc[row,corr_col_number])==False) and (is_nan(temp_m1.iloc[row,construct_col_number])==True):		    
			temp_m1.iloc[row,construct_col_number]=temp_m1.iloc[row,corr_col_number]
			temp_m1.iloc[row,corr_Flag_col_number]=100
	    #If still missing then use ALL fill
	    #Loop through again with Correlate ALL 
	    for ID in new.index:
		if is_nan(temp_m1.iloc[row,construct_col_number])==True:
		   #Create a new variable		
		    corr_label_All=str(VarToCorrelate+"_AllCorr"+"_"+ID)
		    corr_All_col_number=temp_m1.columns.get_loc(corr_label_All)
		    temp_m1.iloc[row,construct_col_number]=temp_m1.iloc[row,corr_All_col_number]  
		    temp_m1.iloc[row,corr_Flag_col_number]=999

	print "Annual Stats"
	print "Mean ", temp_m1[construct_label].mean()
	print "Count ", temp_m1[construct_label].count()
		
    elif corr_freq=='monthly':

	for row in range(0,(temp_m1.shape[0])):
	    # Check to see if the Current row of variable to contrustct a variable is a number. 
	    # Call function and fill if necessary 
	    #Get the column number from the name
	    construct_col_number=temp_m1.columns.get_loc(construct_label)
	     
	    #Use iloc indexing to get the points in the dataframe
		
	    if is_nan(temp_m1.iloc[row,construct_col_number])==True:
		#Loop through BoM ID's first
		#Fill in using each of the ID stations in the list	
		#This should fill using best station in the list.  If any Nans then
		#they should be filled by the next statio, and then the next
		#Assign the AWS station ID as the flagfor ID in new.index:
		#Define Flag variable name and column number
		corr_label_Mon_Flag=str(VarToCorrelate+"_Con_QCFlag")
		corr_Flag_col_number=temp_m1.columns.get_loc(corr_label_Mon_Flag)
		
		#If Construct is missing AND BoM is present then fill and and flag value
		for ID in new.index:
		    corr_label_Mon=str(VarToCorrelate+"_MonCorr"+"_"+str(ID))
		    corr_col_number=temp_m1.columns.get_loc(corr_label_Mon)		    
		    if (is_nan(temp_m1.iloc[row,corr_col_number])==False) and (is_nan(temp_m1.iloc[row,construct_col_number])==True):		    
			temp_m1.iloc[row,construct_col_number]=temp_m1.iloc[row,corr_col_number]
			temp_m1.iloc[row,corr_Flag_col_number]=30
	    #If still missing then use ALL fill
	    #Loop through again with Correlate ALL 
	    for ID in new.index:
		if is_nan(temp_m1.iloc[row,construct_col_number])==True:
		   #Create a new variable		
		    corr_label_All=str(VarToCorrelate+"_AllCorr"+"_"+ID)
		    corr_All_col_number=temp_m1.columns.get_loc(corr_label_All)
		    temp_m1.iloc[row,construct_col_number]=temp_m1.iloc[row,corr_All_col_number]  
		    temp_m1.iloc[row,corr_Flag_col_number]=999
		

	print "Monthly Stats"
	print "Mean ", temp_m1[construct_label].mean()
	print "Count ", temp_m1[construct_label].count()
	
    #Do ALL Counts	
    #=================================================
    #Counts now done in diagnostics
	
    if VarToCorrelate=='Ta':
	yunits='oC' 
    elif VarToCorrelate=='Ah':
        yunits='g m-3'
    elif VarToCorrelate=='Ws':
        yunits='m s-1' 
    elif VarToCorrelate=='Ws_CSAT':
        yunits='m s-1' 
    elif VarToCorrelate=='P':
	yunits='kPa'
    elif VarToCorrelate=='ps':
	yunits='kPa'
    else:
	    yunits=' '   	
    #Do some plots.  Call the routines
    #Plot the 30 minute data
    title4plot="Tower and Constructed 30 minute"
    xdata1=temp_m1[construct_label];     xdata2=temp_m1[fill_label]
    ylabel=VarToCorrelate ;     xlabel="Time" ;     xunits="months"
    plotstuff.plotandPDF2timeseries(title4plot,xdata1,xdata2,xlabel,ylabel,xunits,yunits,Site_ID,VarToCorrelate,mypathforResults)
    
    #Plot the monthly
    title4plot="Tower and Constructed monthly averages"
    xdata1=temp_m1[construct_label].groupby([lambda x: x.year,lambda x: x.week]).mean()     
    xdata2=temp_m1[fill_label].groupby([lambda x: x.year,lambda x: x.week]).mean()
    ylabel=VarToCorrelate ;     xlabel="Time" ;    xunits="months"
    plotstuff.plotandPDF2timeseries(title4plot,xdata1,xdata2,xlabel,ylabel,xunits,yunits,Site_ID,VarToCorrelate,mypathforResults)

    ############################################################
    # Create a variable that is the best correlated output. OUtput in entirety                                #
    ############################################################
    #Create a new label for pandas df column for the contructed variable 
    corr_label=str(VarToCorrelate+"_Corr")
    temp_m1[corr_label]=np.nan
    for ID in new.index:
	corr_label_Ann=str(VarToCorrelate+"_AnnCorr"+"_"+str(ID))
	temp_m1[corr_label][temp_m1[corr_label].isnull()]=temp_m1[corr_label_Ann]    

    ###########################################
    # Call function to do climatology gap fill#
    ###########################################
    #If there are any values still missing use climatology to gap fill
    temp_m1=gap_fill_climatology.climatology_monthly_diurnal(temp_m1,VarToCorrelate)

    
    #Write out file
    print "Writing out files for "+VarToCorrelate+" at "+Site_ID+ " at frequency "+ corr_freq
    fname_results=mypathforResults+'/'+'Results Tower and BoM AWS Variable '+VarToCorrelate+ ' at ' +Site_ID+'.csv'
    temp_m1.to_csv(fname_results)
	
    print "Finished Meteorological Gap filling for "+VarToCorrelate+" at "+Site_ID+ " at frequency "+ corr_freq
    
    con_label_Flag=str(VarToCorrelate+"_Con_QCFlag")
	
    #return a dataframe with two variables
    return temp_m1[[construct_label,con_label_Flag,corr_label]]