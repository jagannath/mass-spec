#! home/jaggu/anaconda/python2.7

"""
A python script entirely devoted to drawing spectra. I can't say how well
contained it will be and how much script repetion is needed, but it definitely
cleans up a lot of long unnecessary scripts that can be kept separate.
"""

from __future__ import division
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import pylab
import os
import sys
from operator import itemgetter 

check = "drawSpectra.py successfully imported"

class Figure(object):
    """
    >>spectra = Figure(DataFrame,params)
    >>spectra.draw(fname,type=['all','normCID','normUVPD','NegOdds'])
    >>spectra.baseline(type=['mean','median','normAvg'])
    >>spectra.label()
    >>
    if all - Odds Ratio
            -------------
            CID | UVPD
            
    """
    def __init__(self,df,size=( 20,20)):
        self.df = df
        self.fig = plt.figure(1, figsize=(size))
        self.fig.subplots_adjust(hspace=0.45,wspace=0.3)
        self.xaxis = df.index
        self.xlim = [self.xaxis[0],self.xaxis[-1]]

    def draw(self,col='all'):
        self.choice = list()
        xlabel = 'mz'
        ylabel = 'Normalized Intensity'
        def _subfigure(ax,ycol):
            yaxis = self.df[ycol]
            rects = ax.bar(self.xaxis,yaxis,width=0.35)
            ax.set_xlim(self.xlim)
            ax.set_ylim([0,np.max(yaxis)*1.2])
            ax.set_title(ycol)
            ax.axhline(y=0,color='black')
            return ax,rects

        if col == 'all':
            ax1 = self.fig.add_subplot( 211,xlabel=xlabel,ylabel='negative Odds Ratio')
            ax1,rects1 = _subfigure(ax1,'NegOdds')
            
            ax2 = self.fig.add_subplot( 223,xlabel=xlabel,ylabel=ylabel)
            ax2,rects2 = _subfigure(ax2,'normCID')

            ax3 = self.fig.add_subplot( 224,xlabel='mz',ylabel=ylabel)
            ax3,rects3 = _subfigure(ax3,'normUVPD')
            self.choice = [col, [(ax1,rects1),(ax2,rects2),(ax3,rects3)]]

        elif col in ['normCID','normUVPD','NegOdds']:
            ax = self.fig.add_subplot( 111,xlabel=xlabel,ylabel=ylabel)
            ax,rects = _subfigure(ax,col)
            self.choice  = [col, [(ax,rects)]]

        else: raise SystemExit("Invalid type (use only NegOdds,normCID,normUVPD,all")

    def baseline(self,method='normAvgMZ',std=1):
        assert self.choice[0] in ['all','normCID','normUVPD','NegOdds']
        xlen = self.df['CID'].count()
        def _drawLine(ax,col,std=std):
            if method == 'mean':
                meancol = np.mean(self.df[col])
                ax.axhline(y=meancol,color='red')
            elif method == 'median':
                mediancol = np.median(self.df[col])
                ax.axhline(y=meancol,color='red')
            elif method == 'normAvgMZ':
                mu = np.mean(self.df[col])
                std = std
                distr = scipy.stats.norm(mu,std)
                zlist = np.arange(-2,2,( 4/xlen))
                yaxis = [distr.pdf(z) for z in zlist]
                line = lines.Line2D(self.xaxis,yaxis,color='red')
                ax.add_line(line)
            elif method == 'normAvgINT':
                mu = np.mean(self.df[col])
                distr = scipy.stats.norm(mu,std)
                def _criteria(s): return distr.pdf(s)
                yaxis = self.df.apply(_criteria,axis=1)
                line = lines.Line2D(self.xaxis,yaxis,color='red')
                ax.add_line(line)
            else: raise SystemExit("Invalid method (use only mean,median or normAvg")
            return ax

        if self.choice[0] == 'all':
            [(ax1,rects1),(ax2,rects2),(ax3,rects3)] = self.choice[1]
            ax2 = _drawLine(ax2,'normCID',std=std)
            ax3 = _drawLine(ax3,'normUVPD',std=std)
        else:
            col = self.choice[0]
            [(ax,rects)] = self.choice[1]
            if not (col is 'NegOdds'): ax = _drawLine(ax,col,std=std)


    def label(self,mustLabel=list(),cutoff=10):
        # mustLabel can be passed as a list of mz which then has to be labelled.
        assert self.choice[0] in ['all','normCID','normUVPD','NegOdds']
        def _reducelabel(xylist):
            def _window(pos,neighbours):
                try:
                    if (xylist[pos+1][0]-xylist[pos][0]) < 5:
                        neighbours.append(xylist[pos+1])
                        return _window(pos+1,neighbours)
                    else: 
                        return neighbours,pos+1
                except IndexError:
                    return neighbours, pos+1
            i = 0
            redLabels = list()
            while i < len(xylist):
                neighbours = [xylist[i]]
                neighbours,i = _window(i,neighbours)
                sortedList = sorted(neighbours,key=itemgetter(1),reverse=True)
                (maxMZ,maxINT) = sortedList[0]
                redLabels.append(maxMZ)
            return redLabels

        def _autolabel(ax,rects):
            valList = [(rect.get_x(),rect.get_height()) for rect in rects if rect.get_height()>10]
            mincut = (max(rect.get_height() for rect in rects)*cutoff)/100
            reducedLabels = _reducelabel(valList)            
            for rect in rects:
                ht = rect.get_height()
                mz = rect.get_x()
                if (ht>mincut and mz in reducedLabels):
                    label = ( ax.text(mz+rect.get_width()/2.,1.1*ht,"%d"%(int(mz)),ha='center',va='bottom',color='r') )
                    label.set_size(8)
            return ax

        if self.choice[0] == 'all':
            [(ax1,rects1),(ax2,rects2),(ax3,rects3)] = self.choice[1]
            [_autolabel(ax,rects) for ax,rects in self.choice[1]]
        else:
            col = self.choice[0]
            [(ax,rects)] = self.choice[1]
            _autolabel(ax,rects)


if __name__ == '__main__':
    print (__doc__)
