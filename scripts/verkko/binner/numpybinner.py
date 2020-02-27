from __future__ import absolute_import

import numpy as np
import binner


def getLimits2D(x, y=None, xscale='log', yscale='log'):
    xidx = (x != np.nan)
    if xscale == 'log':
        xidx *= (x != 0)
    if y is not None:
        yidx = (y != np.nan)
        if yscale == 'log':
            yidx = (y != 0)
        return xidx, yidx
    return xidx


def densityPlot(x,y,bins2D=None,bins=None,xscale='log',yscale='log',nbinsX=1.5,nbinsY=1.5,weights=None):
    """ 
    2-Plot static distance vs temporal distance
    Plot fraction of Finite distances vs temporal distance
    Plot the distance distribution
    """
    x=np.array(x) 
    y=np.array(y)
    if bins2D==None or bins==None :
        xidx,yidx=getLimits2D(x,y,xscale,yscale)
        minValueX, maxValueX = min(x[xidx]), max(x[xidx]) 
        minValueY, maxValueY = min(y[yidx]), max(y[yidx])
        print minValueX, maxValueX, minValueY, maxValueY
        bins2D = binner.Bins2D(float, minValueX, maxValueX, xscale, nbinsX, float, minValueY, maxValueY, yscale, nbinsY)
        bins=binner.Bins(float, minValueX, maxValueX, xscale, nbinsX)
    
    binned_data,b1,b2=np.histogram2d(x,y,bins=bins2D.bin_limits,weights=weights)
    binned_data=binned_data.astype(np.float64)
    X, Y = bins2D.edge_grids
    z,b=np.histogram(x,bins=bins.bin_limits,weights=weights)
    for i in np.nonzero(z)[0]:
        binned_data[i,:]=binned_data[i,:]/z[i]
    binned_data = binned_data.transpose()
    binned_data = np.ma.masked_array(binned_data, binned_data == 0)
    if weights==None :
        z1,b=np.histogram(x,bins=bins.bin_limits,weights=y)
    else :
        z1,b=np.histogram(x,bins=bins.bin_limits,weights=y*weights)
    z1=z1.astype(np.float64)
    for i in np.nonzero(z)[0]:
        z1[i]=z1[i]/z[i] 
    return X, Y, binned_data, bins, z1

def createBins(x,minValueX=None,maxValueX=None,xscale='log',nbinsX=1.5):
    xidx=getLimits2D(x,xscale=xscale)
    if minValueX==None :
        minValueX=min(x[xidx])
    if maxValueX==None :
        maxValueX=max(x[xidx])
    bins=binner.Bins(float, minValueX, maxValueX, xscale, nbinsX)
    return bins

def sumDivide(x,y,minValueX=None,maxValueX=None,xscale='log',nbinsX=1.5):
    """ 
    1-Plot 
    Bins.bin_sum_divide
    """
    x=np.array(x) 
    y=np.array(y)
    bins=createBins(x,minValueX=minValueX,maxValueX=maxValueX,xscale=xscale,nbinsX=nbinsX)
    z,b=np.histogram(x,bins=bins.bin_limits,weights=y,density=True)
    return bins, z

def average(x,y,bins=None,minValueX=None,maxValueX=None,xscale='log',nbinsX=1.5):
    """ 
    1-Plot 
    Bin data and return the average of data points in each bin
    Bins.bin_average
    """
    x=np.array(x) 
    y=np.array(y)
    if bins==None:
        bins=createBins(x,minValueX=minValueX,maxValueX=maxValueX,xscale=xscale,nbinsX=nbinsX)
    z,b=np.histogram(x,bins=bins.bin_limits,weights=y)
    z1,b=np.histogram(x,bins=bins.bin_limits)
    z=z.astype('float')
    for i in np.nonzero(z1)[0]:
        z[i]=z[i]/z1[i] 
    return bins, z

def sum(x,y,bins=None,minValueX=None,maxValueX=None,xscale='log',nbinsX=1.5):
    """ 
    1-Plot 
    Bin data and return the sum of data points in each bin
    Bins.bin_sum
    """
    x=np.array(x) 
    y=np.array(y)
    if bins==None:
        bins=createBins(x,minValueX=minValueX,maxValueX=maxValueX,xscale=xscale,nbinsX=nbinsX)
    z,b=np.histogram(x,bins=bins.bin_limits,weights=y)
    return bins, z

def countDivide(x,bins=None,minValueX=None,maxValueX=None,xscale='log',nbinsX=1.5):
    """ 
    1-Plot 
    Bins.bin_count_divide
    """
    x=np.array(x) 
    if bins==None:
        bins=createBins(x,minValueX=minValueX,maxValueX=maxValueX,xscale=xscale,nbinsX=nbinsX)
    z,b=np.histogram(x,bins=bins.bin_limits,density=True)
    return bins, z


def count(x,bins=None,minValueX=None,maxValueX=None,xscale='log',nbinsX=1.5):
    """ 
    1-Plot 
    Bins.bin_count
    """
    x=np.array(x) 
    if bins==None:
        bins=createBins(x,minValueX=minValueX,maxValueX=maxValueX,xscale=xscale,nbinsX=nbinsX)
    try :
        binLimits=bins.bin_limits
    except AttributeError :
        binLimits=bins
    z,b=np.histogram(x,bins=binLimits)
    return bins, z
