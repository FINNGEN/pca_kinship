import pickle
import matplotlib
import pylab as p

from itertools import imap as imap_itertools
from itertools import izip as izip_itertools

import binner

class imap:
    def __init__(self,f,i):
        self.f=f
        self.i=i
    def __iter__(self):
        return imap_itertools(self.f,self.i)
        
class izip:
    def __init__(self,*args):
        self.args=args
    def __iter__(self):
        return apply(izip_itertools,self.args)


def applySave(function,args,filename,keywords={}):
    """
    Applies the function with arguments, saves the results and
    returns them.
    If the results are alredy saved, the function is not run,
    but only the results are returned.
    """
    try:
        result=pickle.load(open(filename,'r'))        
    except Exception:
        result=function(*args,**keywords)
        outputfile=open(filename,'w')
        pickle.dump(result,outputfile)
        outputfile.close()
    return result

def evalSave(evalstr,filename,globals=None):
    """
    Evaluates the evalstr, saves the results and
    returns them.
    If the results are alredy saved, the evalstr is not evaluated,
    but only the results are returned.
    """
    try:
        result=pickle.load(open(filename,'r'))        
    except Exception:
        if globals!=None:
            result=eval(evalstr,globals)
        else:
            result=eval(evalstr)
        outputfile=open(filename,'w')
        pickle.dump(result,outputfile)
        outputfile.close()
    return result
    


def binData(data,binfunction,bintype,binvalue,minValue=None,maxValue=None,dataType=None):
    #if data is tuple, the first value is the key
    if minValue==None or maxValue==None:
        dataIterator=data.__iter__()
        firstValue=dataIterator.next() #assume that data has more than 1 value
        try:
            firstKey=firstValue[0]
            #data is tuples
            if minValue==None:
                minValue=min(imap(lambda x:x[0],data))
            if maxValue==None:
                maxValue=max(imap(lambda x:x[0],data))
        except Exception:
            #data is not tuples
            firstKey=firstValue
            if minValue==None:
                minValue=min(data)
            if maxValue==None:
                maxValue=max(data)
    if dataType==None:
        dataType=type(firstKey)
    if bintype=="linlog":
        minValue=int(minValue)
    bins=binner.Bins(dataType,minValue,maxValue,bintype,binvalue)
    function=getattr(bins,'bin_'+binfunction)
    binned=function(data)
    return bins,binned

def binDataSave(data,binfunction,bintype,binvalue,name,minValue=None,maxValue=None,dataType=None,dir=''):
    filename=dir+'binned_'+name+'_bintype='+bintype+"_binvalue="+str(binvalue)+".txt"
    return applySave(binData,(data,binfunction,bintype,binvalue),filename,keywords={'minValue':minValue,'maxValue':maxValue,'dataType':dataType})

def binNodePropertyDist(net,function,propertyName,bintype,binvalue):
    values=imap(function,net)
    minValue,maxValue=min(values),max(values)
    
def plotNodeProperties(net,function,propertyName1,propertyName2):
    pass
    
def plotNodePropertyDist(net,function,propertyName,bintype,binvalue,ytype):
    bins,binned=binNodePropertyDist(net,function,propertyName,bintype,binvalue)
    fig=p.figure();ax=fig.add_subplot(111)
    ax.loglog(inDegBins.centers,binnedInDeg,'.')
    

def plotDensityMap(bins,binnedData,xscale='log',yscale='log',normaliseX=True,logScale=True):
    if logScale:
        l, b, r, t = 0.1, 0.12, 1.0, 0.970
    else:
        l, b, r, t = 0.1, 0.12, 1.05, 0.970
    axes_rect = [l,b,r-l,t-b]
    fig=p.figure()
    fig.subplots_adjust(left=0.01, bottom=.05, right=.985, top=.95, wspace=.005, hspace=.05)
    ax = fig.add_axes(axes_rect)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    X,Y=bins.edge_grids
    if normaliseX:
        ySum=p.sum(binnedData,1)
        ySum=p.ma.masked_array(ySum,ySum==0)
        Z=binnedData.transpose()/ySum
    else:
        Z=binnedData.transpose()

    if logScale:
        mappable=ax.pcolor(X,Y,p.ma.array(Z,mask=Z==0),cmap=p.get_cmap("jet"),norm=matplotlib.colors.LogNorm())
    else:
        mappable=ax.pcolor(X,Y,p.ma.array(Z,mask=Z==0),cmap=p.get_cmap("jet"))

    markersize=5.0
    linewidth=2.0

    fig.colorbar(mappable)

    #ax.set_ylim((T/2.0/bins2D.centers[0][-1],T/2.0*2))
    #ax.set_xlim(bins2D.centers[0][0]*day,bins2D.centers[0][-1]*day)
    return fig,ax

