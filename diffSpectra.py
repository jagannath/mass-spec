#! /home/jaggu/anaconda/bin/python2.7

# AUTHOR : JAGANNATH S
# CREATED : 2013-06-17 

from __future__ import division
import sys
import os
import time
import numpy as np
import pylab as p
import collections
from operator import itemgetter
import pandas as pd
import matplotlib.pyplot as plt

class MassSpecData:
    """
    This class takes in the input of the Mass spec data as a csv file (format
    temperory!). Also the number of lines to be eliminated, considered header
    can be provided.A roundoff value is provided for. 
    Generates a dictionary of peaklist (mz:int).
    """
    def __init__(self,fname,headerLines=6,roundoff=0):
        default_DIR = os.getcwd()+'/data/'
        self.fname = default_DIR+fname
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
            mz,intensity = map(float,line.rstrip().split('\t'))
            trunc_mz = float(ROUNDOFF.format(mz))
            mzInt_list_DICT[trunc_mz].append(intensity)
            mzInt_DICT[trunc_mz]+=intensity
        return mzInt_DICT

class ManageData:
    def __init__(self,cid,uvpd):
        """The class handles cid and uvpd peak list data and performs various
        operations on them. On calling it, it returns the original and final
        dataframe """
        self.cid = cid
        self.uvpd = uvpd
    def __call__(self):
        df_original = self.mergeDF()
        self.normalizeDF()
        self.diffDF()
        self.oddsRatioDF()
        self.diffCasesDF()
        df_final = self.remBgrnd()
        return df_original, df_final
    def mergeDF(self):
        cidSer = pd.Series(self.cid.values(),index=self.cid.keys())
        uvpdSer = pd.Series(self.uvpd.values(),index=self.uvpd.keys())
        self.df = pd.DataFrame({'CID':cidSer,'UVPD':uvpdSer})
        self.df.fillna(0) #All NaN is filled to 0
        return self.df
    def normalizeDF(self):
        df = self.df
        self.df['normCID'],self.df['normUVPD'] = map(lambda x:
                                                     (self.df[x]/np.sum(self.df[x])*100),['CID','UVPD'])
    def diffDF(self):
        self.df['normCID_diff'],self.df['normUVPD_diff'] = map(lambda
                                                               x:self.df[x]-np.mean(self.df[x]),['normCID','normUVPD'])
    def oddsRatioDF(self):
        self.df['oddsRatio'] = self.df['normCID_diff']/self.df['normUVPD_diff']

    def diffCasesDF(self):
        def criteria(s):
            cdif,udif = 'normCID_diff','normUVPD_diff'
            if ((s[cdif]>0)&(s[udif]<0)): return 'BION'
            elif ((s[cdif]>0)&(s[udif]>0)): return 'YION'
            elif ((s[cdif]<0)&(s[udif]>0)): return 'INTRSTNG'
            else : return 'BGRND'
        self.df['Ions'] = self.df.apply(criteria,axis=1)
    def remBgrnd(self):
        df2 = self.df[~(self.df['Ions']=='BGRND')]
        return df2

class Output:
    """ For now this takes in the final df. And generates (a) a spectra (only
    bions) and (b) only csv output
    """
    def __init__(self,df,name,inputHandles):
        self.inputs = inputHandles
        self.df = df
        self.idx = df.index
        self.outputDir = os.getcwd()+'/results/'+name
        if not os.path.exists(self.outputDir): os.makedirs(self.outputDir)
    def spectra(self,fname='spectra.png',col='oddsRatio'):
        myfunc = lambda s:-1*s['oddsRatio'] if (s['Ions']=='BION') else s['oddsRatio']
        self.df['BooleanOdds'] = self.df.apply(myfunc,axis=1)
        self.xaxis = self.idx
        self.yaxis = self.df['BooleanOdds']
        fig = plt.figure()
        ax = fig.add_subplot(111,xlabel='mz',ylabel='oddsRatio')
        rects = ax.bar(self.xaxis,self.yaxis,width=0.35)
        ax.axhline(y=0)
        self.autolabel(rects,ax,self.labelList())
        fig.savefig(self.outputDir+'/'+fname)
        return True
    def writeCSV(self,fname):
        # It writes the output csv file to results directory in the current
        # working directory
        path = self.outputDir+'/'+fname+'.xlsx'
        writer = pd.ExcelWriter(path)
        self.df.to_excel(writer,sheet_name='SUMMARY',float_format="%0.2f",header=True,index_label='m/z')
        self.df[self.df['Ions']=='BION'].to_excel(writer,sheet_name='BIONS',float_format="%0.2f",header=True,index_label='m/z')
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
    name = 'LVNELTEFAK'
    uvpd = MassSpecData('UVPD_LVNELTEFAK.csv')
    uvpd_DICT = uvpd.peakList()
    cid = MassSpecData('CID_LVNELTEFAK.csv')
    cid_DICT = cid.peakList()
    data = ManageData(cid_DICT,uvpd_DICT)
    data_original,data_final = data()
      
    msout = Output(data_final,name,[uvpd,cid])
    msout.spectra(fname='spectra.png')
    msout.writeCSV('result')
    
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





