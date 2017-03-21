# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:41:33 2015

@author: imchugh
"""

import windrose
import matplotlib.pyplot as plt
from matplotlib import *
import numpy as np

def new_axes():
    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    rect = [0.1, 0.1, 0.8, 0.8]
    ax = windrose.WindroseAxes(fig, rect, axisbg='w')
    fig.add_axes(ax)
    return ax

def set_legend(ax):
    l = ax.legend(borderaxespad=-0.10)
    plt.setp(l.get_texts(), fontsize=8)
    
def do_plot(wd, ws):
    ax = new_axes()
    ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
    set_legend(ax)

def do_contour_plot(wd, ws):
    ax = new_axes()
    ax.contourf(wd, ws, bins = np.arange(0,8,1), cmap=cm.hot)
    set_legend(ax)