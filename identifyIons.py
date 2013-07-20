#! /home/jaggu/anaconda/python2.7

"""
Very ambitious with no real knowledge of what mathematical functions to use. But
I am just going to go through a peak list of b and y ions and try saying which
looks like a better candidate. It is more likely an important ion is - (a) It has some adjoining peaks
with a mass difference (like 18 for water, 17 for Ammonia, 27 for HCN and 28 for
CO). (b) Not throwing it off based on its intensity. I hope to assign a
probability scores for each of the cases and a combined score over all. 
"""

import sys
import numpy as np
import pandas as pd
import diffSpectra as ds
import drawSpectra 
import time


def score(bion,bionDICT):
    """
    A delta list is provided - Neutral losses. Checks for key in DICT after
    subtracting the value. If KeyError then nothing, else adds a score and
    returns a list of  [Score, [(ion,intensity,delta value),...]]
    """
    score = 0
    deltaList = [1,17,18,27,28]
    bionList = bionDICT.keys()
    intensityList = bionDICT.values()
    meanInt = np.mean(intensityList)
    proofList = list()
    scoreList = list()
    if bionDICT[bion]>meanInt:
        score+=1
        proofList.append([bion,bionDICT[bion],0])
        scoreList = [score,proofList]
    for d in deltaList:
        bdiff = bion - d
        if bdiff in bionList:
            if bionDICT[bdiff]>meanInt:
                proofList.append([bdiff,bionDICT[bdiff],d])
                score+=1
    scoreList = [score,proofList]
    return [score,proofList]

def action(data,normalize='canonical',baseline='mean',std=1):
    name = 'LGEYGFQNALIVR'
    mergedDF = data.mergeDF()
    data.normalizeDF(method=normalize)
    data.subtractBaselineDF(method=baseline)
    data.oddsRatioDF()
    data.distinguishIonsDF()
    data.negOdds()
    df_final = data.remBgrnd()
    fname ='LGEYGFQNALIVR'+'_normalize='+normalize+'_baseline='+baseline+'_std='+str(std)
    out = ds.Output(data.df,dirname='TEST_'+name)
    out.writeCSV(fname,bion=True)
    return mergedDF, data.df, df_final


def main():
    baselineList = ['mean','median','normAvg']
    normalizeList = ['canonical','max','direct','log']

    name = 'LGEYGFQNALIVR'
    dir = '/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/jennyBroadbelt/scott/data/LGEYGFQNALIVR/'
    uvpd = ds.MassSpecData(name+'.UVPD.csv',DIR=dir,headerLines=9)
    uvpd_DICT = uvpd.peakList()
    cid = ds.MassSpecData(name+'.CID.csv',DIR=dir,headerLines=9)
    cid_DICT = cid.peakList()
           
    data = ds.ManageData(cid_DICT,uvpd_DICT)
    #for std in [0.01,0.1,1,10,100]:
    #    print "std %d"%(std)
    #    mergedDF, df, df_final = action(data,normalize='direct',baseline='normAvgMZ',std=std)
    mergedDF, df, df_final = action(data,normalize='direct',baseline='normAvgINT')
    fig = drawSpectra.Figure(data.df)
    fig.draw()
    fig.baseline()
    fig.label()
    print fig
    fig.fig.savefig('spectra_test.png')

#    for baseline in baselineList:
#        for normalize in normalizeList:
#            action(data,normalize=normalize,baseline=baseline)
    
    sys.exit(1)
    iondf = data_final[['CID','Ions']]
    bdf = iondf[iondf['Ions']=='BION']
    bionData =list(bdf.itertuples())
    bionDICT = dict()
    for mz,CID,type in bionData:
        bionDICT[mz] = CID

    bions = [x[0] for x in bionData]
    for b in bions:
        scoreList = score(b,bionDICT)
        if scoreList[0] > 0:
            print b, scoreList
    
if __name__ == '__main__':
    #Global Definition; These must be defined and in this format else programme fails.                                          
    roundoff = 0
    ROUNDOFF = "{0:."+str(roundoff)+"f}"
    blist = [386.1128,515.1544,678.2187,735.2402,882.3086,1010.367,1124.41,1195.447,1308.531,1421.615,1520.684]
    ylist = [260.1499,407.2183,536.2609,637.3086,750.3927,879.4353,993.4782,1092.5466]
   
    bions = [float(ROUNDOFF.format(ion)) for ion in blist]
    yions = [float(ROUNDOFF.format(ion)) for ion in ylist]
    t0 = time.clock()
    main()
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))

