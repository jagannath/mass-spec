#! /home/jaggu/anaconda/bin/python2.7

# AUTHOR : JAGANNATH S
# CREATED : 2013-06-17 

from __future__ import division
import sys
import os
import time
import numpy as np
import pylab as p
import scipy.stats
import collections
from operator import itemgetter
import pandas as pd
import matplotlib.pyplot as plt
check = "Import diffSpectra successful"


class MassSpecData(object):
    """
    This class takes in the input of the Mass spec data as a csv file (format
    temperory!). Also the number of lines to be eliminated, considered header
    can be provided.A roundoff value is provided for. 
    Generates a dictionary of peaklist (mz:int).
    """
    def __init__(self,fname,DIR=os.getcwd()+'/data/',headerLines=6,roundoff=0):
        self.fname = DIR+fname
        self.hLines = headerLines
        self.header,self.lines = self.readFile()
        self.roundoff = roundoff #Default roundoff=0
    def readFile(self):
        with open(self.fname,'r') as f:
            lines = f.readlines()
        return lines[0:self.hLines],lines[(self.hLines+1):]
    def peakList(self,roundoff=0):
        mzInt_list_DICT = collections.defaultdict(list)
        mzInt_DICT = collections.defaultdict(float)
        self.mzInt_DICT = collections.defaultdict(float)
        ROUNDOFF = "{0:."+str(roundoff)+"f}"
        for line in self.lines[1:]:
            try:mz,intensity = map(float,line.rstrip().split('\t'))
            except ValueError:mz,intensity = map(float,line.rstrip().split(','))
            else:print "Something wrong with csv file"
            trunc_mz = float(ROUNDOFF.format(mz))
            mzInt_list_DICT[trunc_mz].append(intensity)
            mzInt_DICT[trunc_mz]+=intensity
        return mzInt_DICT

class MassSpecDICT(object):
    """
    This takes in a mz:intensity DICT as an input and rounds off the mz
    values.It can do more like truncating out regions of mz or providing
    statistics in general.
    Returns back a dictionary
    """
    def __init__(self,DICT):
        self.DICT = DICT
    def roundoff(self,roundoff = 2):
        ROUNDOFF = "{0:."+str(roundoff)+"f}"
        roffDICT = collections.defaultdict(float)
        for mz, intensity in self.DICT.items():
            trunc_mz = float(ROUNDOFF.format(mz))
            roffDICT[trunc_mz]+=intensity
        return roffDICT

class ManageData(object):
    def __init__(self,cid,uvpd):
        """The class handles cid and uvpd peak list data and performs various
        operations on them. On calling it, it returns the original and final
        dataframe """
        self.cid = cid
        self.uvpd = uvpd
    def __call__(self):
        self.mergedDF = self.mergeDF()
        self.normalizeDF()
        self.subtractBaselineDF()
        self.oddsRatioDF()
        self.distinguishIonsDF()
        df_final = self.remBgrnd()
        return self.mergedDF, df_final

    def mergeDF(self):
        cidSer = pd.Series(self.cid.values(),index=self.cid.keys())
        uvpdSer = pd.Series(self.uvpd.values(),index=self.uvpd.keys())
        self.df = pd.DataFrame({'CID':cidSer,'UVPD':uvpdSer})
        self.df = self.df.fillna(0) #All NaN is filled to 0
        return self.df
    
    def normalizeDF(self,inputColNames=['CID','UVPD'],method='canonical'):
        """
        There are a number of methods of normalization that can be applied. The
        types are - 
        canonical = I(jnorm) = Ij/Sum(Ij)
        direct = I(jnorm) = Ij - Imin/(Imax - Imin)
        log = I(jnorm) = log(I)
        returns : pd.Series, pd.Series : Normalized Series
        """
        outColNames = ['norm'+x for x in inputColNames]
        def _normalize(col,method):
            if method == 'canonical':
                outCol = (self.df[col]/np.sum(self.df[col]))*100
            elif method == 'max':
                outCol = (self.df[col]/np.amax(self.df[col]))*100
            elif method == 'direct':
                outCol = ( self.df[col]-np.amin(self.df[col]) )/ ( np.ptp(self.df[col]) )* 100
            elif method == 'log':
                outCol = np.log10(self.df[col]) #Note if -ve value then gives NaN
            else: raise SystemExit("Invalid normalize method ")
            return outCol

        outputCols = [_normalize(col,method=method) for col in inputColNames]
        self.df[outColNames[0]],self.df[outColNames[1]] = outputCols
        return outputCols

    def subtractBaselineDF(self,inputColNames=['normCID','normUVPD'],method='mean',std = 1):
        def _subBaseline(col,method='mean',std=1):
            """
            This function calculates the base line for the column. The method
            could be the mean, median or a normal distribution with a sigma
            value mean for the intensity or mz. 
            returns : pd.Series : Series with background subtracted values
            """
            if method == 'mean':
                baseline = self.df[col] - np.mean(self.df[col])
            elif method == 'median':
                baseline = self.df[col] - np.median(self.df[col])
            elif method == 'normAvgINT':
                # This is the normalizing of the background based on a std; The mean
                # is the mean for the column intensity. The value is higher at
                # mean and tapers at the ends.So high intensity gets lower
                # background subtraction and ones around the mean are higher.  
                mu = np.mean(self.df[col])
                distr = scipy.stats.norm(mu,std)
                def _criteria(s): return s[col] - distr.pdf(s[col])
                baseline = self.df.apply(_criteria,axis=1)
            elif method == 'normAvgMZ':
                # This is the background more in lines of how the spectrum is.
                # Low intensity at lower mz and higher at the center
                xaxis = self.df.index
                zincr = 4/len(xaxis)
                zlist = np.arange(-2,2,zincr)
                mu = np.mean(self.df[col])
                distr = scipy.stats.norm(mu,std)
                xvallist = [distr.pdf(z) for z in zlist]
                baseline = self.df[col] - xvallist
            else:
                raise SystemExit("Invalid Baseline subtraction method")
            return baseline

        outColNames = [x+'_diff' for x in inputColNames]
        outputCols = [_subBaseline(col,method=method) for col in inputColNames]
        self.df[outColNames[0]],self.df[outColNames[1]] = outputCols
        return outputCols
    
    def oddsRatioDF(self):
        self.df['oddsRatio'] = self.df['normCID_diff']/self.df['normUVPD_diff']
    
    def distinguishIonsDF(self):
        def criteria(s):
            cdif,udif = 'normCID_diff','normUVPD_diff'
            if ((s[cdif]>0)&(s[udif]<0)): return 'BION'
            elif ((s[cdif]>0)&(s[udif]>0)): return 'YION'
            elif ((s[cdif]<0)&(s[udif]>0)): return 'INTERESTING'
            else : return 'BGRND'
        self.df['Ions'] = self.df.apply(criteria,axis=1)
   
    def negOdds(self):
        myfunc = lambda s:-1*s['oddsRatio'] if (s['Ions']=='BION') else s['oddsRatio']
        self.df['NegOdds'] = self.df.apply(myfunc,axis=1)

    def remBgrnd(self):
        df2 = self.df[~(self.df['Ions']=='BGRND')]
        return df2

class Output(object):
    """ For now this takes in the final df. And generates (a) a spectra (only
    bions) and (b) xlsx output file
    """
    def __init__(self,df,dirname):
#        self.inputs = inputHandles
        self.df = df
        self.cols = ['CID','UVPD','normCID','normUVPD','normCID_diff','normUVPD_diff','oddsRatio','Ions']
        self.idx = df.index
        self.outputDir = os.getcwd()+'/results/'+dirname
        if not os.path.exists(self.outputDir): os.makedirs(self.outputDir)

    def statistics(self):
        """ Makes a statistics dataframe. Wanted to write a separate sheet, but
        seems difficult """
        df = self.df
        cols = ['CID','UVPD','normCID','normUVPD']
        ions = df['Ions']
        [cidMean,uvpdMean,normCIDMean,normUVPDMean] = map(lambda x:np.mean(df[x]),cols)
        [cidMedian, uvpdMedian,normCIDMedian,normUVPDMedian] = map(lambda x:np.median(df[x]),cols)
        [ctbions,ctyions,ctintrg,ctbgrnd] = map(lambda x: ions[ions==x].count(),['BION','YION','INTERESTING','BGRND'])
        d = ({'CID':pd.Series([cidMean,cidMedian],index=['Mean','Median']),
              'UVPD':pd.Series([uvpdMean,uvpdMedian],index=['Mean','Median']),
              'normCID':pd.Series([normCIDMean,normUVPDMean],index=['Mean','Median']),
              'normUVPD':pd.Series([normUVPDMedian,normUVPDMedian],index=['Mean','Median']),
              'BION':pd.Series([ctbions],index=['Counts']),
              'YION':pd.Series([ctyions],index=['Counts']),
              'INTERESTING': pd.Series([ctintrg],index=['Counts']),
              'BGRND':pd.Series([ctbgrnd],index=['Counts'])})
        return pd.DataFrame(d)

    def spectra(self,fname='spectra.png',col='oddsRatio',bion=None):
        meanCID,meanUVPD  = map(lambda x:np.mean(self.df[x]),['normCID','normUVPD'])
        if bion:
            myfunc = lambda s:-1*s['oddsRatio'] if (s['Ions']=='BION') else s['oddsRatio']
            self.df['BooleanOdds'] = self.df.apply(myfunc,axis=1)
            self.xaxis = self.idx
            self.yaxis = self.df['BooleanOdds']
        else:
            self.xaxis = self.idx
            self.yaxis = self.df['oddsRatio']
            print self.yaxis
        fig = plt.figure()
        ax = fig.add_subplot(111,xlabel='mz',ylabel='oddsRatio')
        rects = ax.bar(self.xaxis,self.yaxis,width=0.35)
        ax.set_xlim([0,1200])
        ax.axhline(y=0)
        ax.axhline(y=meanCID,color='red')
        ax.axhline(y=-1*meanUVPD,color='magenta')
        #self.autolabel(rects,ax,self.labelList())
        fig.savefig(self.outputDir+'/'+fname)
        return True
    
    def writeCSV(self,fname,bion=None):
        # It writes the output csv file to results directory in the current
        # working directory
        path = self.outputDir+'/'+fname+'.xlsx'
        writer = pd.ExcelWriter(path)
        statdf = self.statistics()
        self.df.to_excel(writer,sheet_name='SUMMARY',cols = self.cols,float_format="%0.2f",header=True,index_label='m/z')
        if bion: self.df[self.df['Ions']=='BION'].to_excel(writer,sheet_name='BIONS',float_format="%0.2f",header=True,index_label='m/z')
        #statdf.to_excel(writer,sheet_name='STATISTICS',header=True)
        writer.save()
        #self.df.to_csv(self.outputDir+'/'+fname+'.csv',sep='\t',header=True,index=True,index_label='m/z')
        #self.df[self.df['Ions']=='BION'].to_csv(self.outputDir+'/'+fname+'.bion.csv',sep='\t',header=True,index=True,index_label='m/z')
        return True

    def autolabel(self,rects,ax,reducedLabels):
        # Attaches some text labels
        for rect in rects:
            height = rect.get_height()
            mz = rect.get_x()
            if height>10 and mz in reducedLabels :
                label = ax.text(rect.get_x()+rect.get_width()/2.,1.1*height,'%d'%int(mz),ha='center',va='bottom',color='r')
                label.set_size(8)
    def labelList(self):
        # Creates and returns a reduced label list
        i = 0
        reducedLabels = list()
        xlist,ylist = self.xaxis,self.yaxis
        def window(pos,neighbours):
            try:
                if xlist[pos+1]-xlist[pos]<2:
                    neighbours.append([xlist[pos+1],ylist.iloc[pos+1]])
                    return window(pos+1,neighbours)
                else:
                    return neighbours,pos+1
            except IndexError:
                return neighbours, pos+1
        while i < len(xlist):
            neighbours = [(xlist[i],ylist.iloc[i])]
            neighbours,i = window(i,neighbours)
            sortedList = sorted(neighbours, key=itemgetter(1),reverse=True)
            (maxMZ,maxINT) = sortedList[0]
            reducedLabels.append(maxMZ)
        return reducedLabels


def main():
    name = 'LGEYGFQNALIVR'
    dir = '/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/jennyBroadbelt/scott/data/LGEYGFQNALIVR/'
    uvpd = MassSpecData(name+'.UVPD.csv',DIR=dir,headerLines=9)
    uvpd_DICT = uvpd.peakList()
    cid = MassSpecData(name+'.CID.csv',DIR=dir,headerLines=9)
    cid_DICT = cid.peakList()
    data = ManageData(cid_DICT,uvpd_DICT)
    data_original,data_final = data()
      
    msout = Output(data_final,name)
    msout.spectra(fname='spectra.'+name+'.png')
    msout.writeCSV('result.mean.'+name)
    
if __name__ == '__main__':
    #Global Definition; These must be defined and in this format else programme fails.                                          
    roundoff = 0
    ROUNDOFF = "{0:."+str(roundoff)+"f}"
    blist = [428.1598,542.2027,671.2453,784.3293,885.3770,1014.4196,1161.4880,1232.5251]
    ylist = [260.1499,407.2183,536.2609,637.3086,750.3927,879.4353,993.4782,1092.5466]
   
    bions = [float(ROUNDOFF.format(ion)) for ion in blist]
    yions = [float(ROUNDOFF.format(ion)) for ion in ylist]
    t0 = time.clock()
    main()
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))





