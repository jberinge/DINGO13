# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 17:00:07 2014

Pass dictionary containing the following parameters:
    - 'xlab' - xlabel
    - 'ylab' - ylabel
    - 'title'

@author: imchugh
"""

import matplotlib.pyplot as plt
import numpy as np
import pdb

def line_plot(x_series, y_series, **options):
    """
    Produces single axis plot;
    x_series must be a 1D array (or pandas series)
    y_series may be 2D, but first dimension must match x_series
    Keyword args are:
    'colors': specify color cycle (otherwise default) (list of color strings)
    'xlim': limits of x axis (list of ints)
    'ylim': limits of y axis (list of ints)
    'title': plot title (string)
    'vert_line': vertical lines on plot
    'hor_line': horizontal lines on plot
    'var_labels': variable names to be shown in legend
    'line_style': line styles to be used (must be of same length as # of series)
    'line_width': line widths to be used (must be of same length as # of series)
    """
    
    # Instantiate plot
    fig = plt.figure(figsize=(16,8))
    fig.patch.set_facecolor('white')
    ax = plt.gca()

    # Set colors
    if 'colors' in options.keys():
        ax.set_color_cycle(options['colors'])
    
    # Generate and plot line objects
    lineObjects = plt.plot(x_series,y_series)

    # Set limits
    if 'xlim' in options.keys():
        ax.set_xlim(options['xlim'])
    if 'ylim' in options.keys():
        ax.set_ylim(options['ylim'])

    # Set title layout    
    if 'title' in options.keys():    
        plt.title(options['title'],fontsize=30,y=1.03)
        
    # Set axis labels        
    if 'xlab' in options.keys():
        ax.set_xlabel(options['xlab'], fontsize = 30)
        ax.xaxis.set_label_coords(0.5, -0.065)
    if 'ylab' in options.keys():
        ax.set_ylabel(options['ylab'], fontsize = 30)
        ax.yaxis.set_label_coords(-0.06, 0.5)
        
    # Set reference lines
    if 'vert_line' in options.keys():
        for i in options['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    if 'hor_line' in options.keys():
        for i in options['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')

    # Set line properties
    if 'line_style' in options.keys():
        if len(options['line_style']) != len(lineObjects):
            print 'Wrong number of line styles supplied for assignment to series'
        else:
            for i, line in enumerate(lineObjects):
                plt.setp(line, linestyle = options['line_style'][i])

    if 'line_width' in options.keys():
        if len(options['line_width']) != len(lineObjects):
            print 'Wrong number of line widths supplied for assignment to series'
        else:
            for i, line in enumerate(lineObjects):
                plt.setp(line, linewidth = options['line_width'][i])

    # Set legend
    if 'var_labels' in options.keys():
        if len(options['var_labels']) != len(lineObjects):
            print 'Wrong number of variable labels supplied for assignment to series'
        else:
            for i, line in enumerate(lineObjects):
                plt.setp(line, label = options['var_labels'][i])
            plt.legend(fontsize=16, loc = [0.9,0.7], frameon=False)
    
    # Other    
    plt.tick_params(labelsize=18)
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.show()

def bar_plot(x_series, y_series, **options):
    
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111)
    
    # Set bar width
    if 'width' in options.keys():
        width = options['width']
    else:
        width = 0.3
    
    # Create list for bar object containers
    barObjects = []
    
    # Do the plotting
    if len(np.shape(y_series)) != 1:
        num_series = np.shape(y_series)[1]
        for n in range(num_series):
            barObjects.append(ax.bar(x_series + width * n, y_series[:, n], width))
    else:
        barObjects.append(ax.bar(x_series + width, y_series, width))

    # Set title layout    
    if 'title' in options.keys():    
        ax.set_title(options['title'],fontsize=30,y=1.03)
        
    # Set axis labels        
    if 'xlab' in options.keys():
        ax.set_xlabel(options['xlab'], fontsize = 20)
        ax.xaxis.set_label_coords(0.5, -0.065)
    if 'ylab' in options.keys():
        ax.set_ylabel(options['ylab'], fontsize = 20)
        ax.yaxis.set_label_coords(-0.07, 0.5)

    # Set category labels
    if 'xticklabs' in options.keys():
        ax.set_xticks(x_series + width)
        ax.set_xticklabels(options['xticklabs'], fontsize = 12)

    # Set bar properties    
    if 'color' in options.keys():
        for n, series in enumerate(barObjects):
            for obj in series:
                obj.set_color(options['color'][n])
    
    # Set legend
    if 'var_labels' in options.keys():
        if len(options['var_labels']) != len(barObjects):
            print 'Wrong number of variable labels supplied for assignment to series'
        else:
            for i, series in enumerate(barObjects):
                series.set_label(options['var_labels'][i])
            plt.legend(fontsize=12, loc = [0.8,0.7], frameon=False)
    
    plt.show()

def bar_example(x_series, y_series, **options):

    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111)
    
    N = 5
    menMeans = y_series[:,0]
    menStd =   (2, 3, 4, 1, 2)
    
    width = 0.35       # the width of the bars
    
    if len(np.shape(y_series)) != 1:
        num_series = np.shape(y_series)[1]
        for n in range(num_series):
            barObjects.append(ax.bar(x_series + width * n, y_series[:, n]))
    else:
        barObjects.append(ax.bar(x_series + width, y_series))    
    
#    pdb.set_trace()
    rects1 = ax.bar(x_series, menMeans, width, color='r', yerr=menStd)

    
    
    womenMeans = y_series[:,1]
    womenStd =   (3, 5, 2, 3, 3)
    rects2 = ax.bar(x_series+width, womenMeans, width, color='y', yerr=womenStd)
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Scores')
    ax.set_title('Scores by group and gender')
    ax.set_xticks(x_series+width)
    ax.set_xticklabels( ('G1', 'G2', 'G3', 'G4', 'G5') )
    
    ax.legend( (rects1[0], rects2[0]), ('Men', 'Women') )
    
    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
    
    autolabel(rects1)
    autolabel(rects2)
    
    plt.show()


def line_plot_b(df,d):
    fig=plt.figure(figsize=(16,8))
    fig.patch.set_facecolor('white')
    for i,var in enumerate(df.columns):
        if 'colors' in d:
            plt.plot(df.index,df[var],label=var,color=d['colors'][i])
        else:
            plt.plot(df.index,df[var],label=var,linewidth=1)
    plt.title(d['title'],fontsize=24)
    plt.xlabel(d['xlab'],fontsize=18)
    plt.ylabel(d['ylab'],fontsize=18)
    plt.tick_params(labelsize=14)
    plt.xticks(rotation=45)
    if 'vert_line' in d:
        for i in d['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    if 'hor_line' in d:
        for i in d['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')
    plt.legend(loc='lower left')
    plt.show()

def line_plot_w_points(df,d):
    
    fig=plt.figure(figsize=(16,8))
    fig.patch.set_facecolor('white')
    mark_count=0
    for i,var in enumerate(df.columns):
        if not var in d['as_markers']:
            plt.plot(df.index,df[var],label=var,color=d['colors'][i],linewidth=1)
        else:
            plt.plot(df.index,df[var],marker=d['marker_style'][mark_count],label=var,color=d['colors'][i],markersize=10)
            mark_count +=1
    plt.title(d['title'],fontsize=24)
    plt.xlabel(d['xlab'],fontsize=18)
    plt.ylabel(d['ylab'],fontsize=18)
    plt.tick_params(labelsize=14)
    plt.xticks(rotation=45)
    for i in d['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    for i in d['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')
    plt.legend(loc='lower right')
    plt.show()    

def one_ax_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}

    d['title']='Whroo Cumulative NEE'
    d['xlab']='\nDate'
    d['ylab']= '$NEE\/(gC\/m^{-2}d^{-1})$'
    d['colors'] = ['b','r','g','m','c','y','k'] # ['b','g','r','c','m','y','k'] # matplotlib default
    # d['as_markers']=['Lai_Canon_hemi','Lai_LAI2000']
    d['marker_style']=['o','s']
    d['vert_line']=[]#['2012-04-01','2012-07-01','2012-10-01','2013-01-01','2013-04-01','2013-07-01','2013-10-01']
    d['hor_line']=[0]
    
    return d

def fill_line_plot(df,d):
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    x=df.index    
    y1=df[df.columns[0]]
    y2=df[df.columns[1]]
    plt.xlim((0,24))
    plt.plot(x,y1,color='black',linestyle='-.',linewidth=2,label=df.columns[0])
    plt.plot(x,y2,color='black',linestyle='-',linewidth=2,label=df.columns[1])
    plt.fill_between(x, y1, y2, where=y2>=y1, facecolor='blue', edgecolor='None',interpolate=True)
    plt.fill_between(x, y1, y2, where=y1>=y2, facecolor='red', edgecolor='None',interpolate=True)
    plt.xlabel(d['xlab'],fontsize=18)
    plt.ylabel(d['ylab'],fontsize=18)
    plt.title(d['title'],fontsize=24)
    plt.xticks([0,4,8,12,16,20,24])
    if 'vert_line' in d:
        plt.axvline(x=d['vert_line'],color='black',linestyle='-')
    if 'hor_line' in d:
        plt.axhline(y=d['hor_line'],color='black',linestyle='-')
    plt.legend(loc='lower right')
    
def two_ax_line_plot(df,d):
    """
    Pass dataframe (df) and dictionary (d) specified by function 
    'two_ax_dict_set' - see below;
    returns open plot
    """
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    ax = plt.gca()
    ax2 = ax.twinx()
    y1=d['y1']
    y2=d['y2']
    for i,var in enumerate(y1['vars']):
        if 'colors' in y1:
            ax.plot(df.index,df[var],label=var,color=y1['colors'][i],linewidth=1)
        else:
            ax.plot(df.index,df[var],label=var,linewidth=1)
    for i,var in enumerate(y2['vars']):
        if 'colors' in y2:
            ax2.plot(df.index,df[var],label=var,color=y2['colors'][i],linewidth=1)
    if 'vert_line' in d['globals']:
        for i in d['globals']['vert_line']:
            plt.axvline(x=i,color='black',linestyle='dotted',linewidth=0.5)
    if 'hor_line' in d['globals']:
        for i in d['globals']['hor_line']:
            plt.axhline(y=i,color='black',linestyle='-')        
    plt.title(d['globals']['title'],fontsize=24)
    ax.set_xlabel(d['globals']['xlab'],fontsize=18)
    ax.set_ylabel(y1['lab'],fontsize=18)
    ax2.set_ylabel(y2['lab'],fontsize=18)
    plt.tick_params(labelsize=14)
    ax.legend(loc='upper left',fontsize=14)
    ax2.legend(loc='upper right',fontsize=14)
    plt.show()

def two_ax_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}

    d['globals']={}
    d['globals']['title']='GPP and ET Riggs\n'
    d['globals']['xlab']='\nDate'
    d['globals']['vert_line']=['2012-04-01','2012-07-01','2012-10-01','2013-01-01','2013-04-01','2013-07-01','2013-10-01']

    d['y1']={}
    d['y1']['vars']=['r_GPP_rm']
    d['y1']['colors']=['green']
    d['y1']['lab']='$GPP\/(gC\/m^{2}d^{-1})$'

    d['y2']={}
    d['y2']['vars']=['r_ET_rm']
    d['y2']['colors']=['blue']
    d['y2']['lab']='$ET\/(gH_{2}O\/m^{-2}d^{-1})$'
        
    return d

def stacked_bar_plot(df,d):
    """    
    Pass dataframe (df) and dictionary (d) specified by function 
    'stacked_bar_dict_set' - see below;
    returns open plot
    """
    fig=plt.figure(figsize=(12,8))
    fig.patch.set_facecolor('white')
    ax = plt.gca()
    width=0.35
    plt.xlim(0.85,1.5)
    sum_vars=[]
    for i,var in enumerate(df.columns):
        sum_vars==[] if i==0 else sum_vars.append(df.ix[0,i-1])
        if 'colors' in d:
            plt.bar(df.index+1, df[var], width, color=d['colors'][i], label=var,bottom=sum(sum_vars))
    ax.axes.get_xaxis().set_ticks([])        
    plt.title(d['title'],fontsize=24)
    plt.ylabel(d['ylab'],fontsize=18)
    if 'xlab' in d: plt.ylabel(d['xlab'],fontsize=18)
    if 'ylab' in d: plt.ylabel(d['ylab'],fontsize=18)    
#    plt.legend(loc='upper right')
    plt.text(1.175,36,'Aboveground Vegetation:\n 37.75', fontsize=18,horizontalalignment='center',verticalalignment='center')
    plt.text(1.175,15,'Litter: 5.8', fontsize=18,horizontalalignment='center',verticalalignment='center')
    plt.text(1.175,5,'Belowground vegetation:\n 10.74', fontsize=18,horizontalalignment='center',verticalalignment='center')
    bbox_props = dict(boxstyle="larrow,pad=0.3", ec="black", fc='white', lw=2)
    plt.text(1.38,11,'Soil: 1.69',fontsize=18,bbox=bbox_props)
    fig.show()
    
def stacked_bar_dict_set():
    """
    Returns dictionary required to run the two y axis algorithm above
    """
    d={}
    
    d['title']='Whroo carbon pools\n'
    d['colors']=['blue','brown','grey','green']
    d['ylab']='$C\/storage\/(tC\/ha^{-1})$'
    
    return d