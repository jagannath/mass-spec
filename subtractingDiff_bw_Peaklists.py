#! usr/bin/python

"""
Scott of Dr Broadbelt's group has two peaklists for the same peptide, one originating from a CID which produces both b and y ions and another from UVPD that produces only the y ions.
The idea is to see differences in the peak lists and the ones that are more are the b ions. 
"""
from __future__ import division
import sys
import os
import time
import numpy as np
import pylab as p
import collections
from operator import itemgetter


def getPeakList(f):
    ifile = open(f,'r')
    lines = ifile.readlines()
    ifile.close()
    mzInt_list_DICT = collections.defaultdict(list)
    mzInt_DICT = dict()
    normMz_DICT = dict()
    for line in lines:
	if not line.startswith("\""):
	    mz, intensity = map(float,line.rstrip().split('\t'))
	    trunc_mz = float("{0:.0f}".format(mz))
	    mzInt_list_DICT[trunc_mz].append(intensity)

    for k, v in mzInt_list_DICT.items():
	mzInt_DICT[k] = sum(v) 
	if k in yions: print "Y ion: %d - %f"%(k,sum(v))
	if k in bions: print "B ion: %d - %f"%(k,sum(v))

    maxVal = max(mzInt_DICT.values())
    for k, v in mzInt_DICT.items():
	normInt = (v/maxVal)*100
	normMz_DICT[k] = normInt
    return normMz_DICT

def window(lst,pos,neighbors):
    try:
	if lst[pos+1] - lst[pos] < 2:
	    neighbors.append(lst[pos+1])
	    return window(lst,pos+1,neighbors)
	else:
	    return neighbors,pos+1
    except IndexError:
	return neighbors, pos+1

def reduceLabels(mzlst,pksList_DICT):
    lst = sorted(mzlst)
    neighborPeaks = list()
    i = 0
    lastpos = len(lst)
    reducedLabels = list()
    while i<lastpos:
	neighbors = [lst[i]]
	neighbors,i = window(lst,i,neighbors)
	mzIntList = [(mz,pksList_DICT[mz]) for mz in neighbors]
	sortedLst = sorted(mzIntList, key=itemgetter(1),reverse=True)
	(maxMZ,maxINT) = sortedLst[0]
	reducedLabels.append(maxMZ)

    return reducedLabels
    

def autolabel(rects,ax,pksList_DICT):
    # attach some text labels
    allLabels = [rect.get_x() for rect in rects if rect.get_height()>5]
    labels = reduceLabels(allLabels,pksList_DICT)
    for rect in rects:
	height = rect.get_height()
	if height> 10 and rect.get_x() in labels:
	    mz = rect.get_x()
	    if mz in yions:
		ax.text(rect.get_x()+rect.get_width()/2.,1.5*height, '%d'%(mz),ha='center', va='bottom',color='green')
	    elif mz in bions:
		ax.text(rect.get_x()+rect.get_width()/2.,1.5*height, '%d'%(mz),ha='center', va='bottom',color='blue')
	    else:
		ax.text(rect.get_x()+rect.get_width()/2.,1.2*height, '%d'%(mz),ha='center', va='bottom',color='red')
	    print rect.get_x(), height
    return allLabels

def main():
    #f_CID = 'UVPD_2pulses.csv'
    f_CID = 'CID_LVNELTEFAK.csv'
    f_UVPD = 'UVPD_LVNELTEFAK.csv'
    print "CID"
    pksList_CID_DICT = getPeakList(f_CID)
    print "UVPD"
    pksList_UVPD_DICT = getPeakList(f_UVPD)
    mz_CID = pksList_CID_DICT.keys()
    mz_UVPD = pksList_UVPD_DICT.keys()
    sumAllINT_CID = sum(pksList_CID_DICT.values())
    sumAllINT_UVPD = sum(pksList_UVPD_DICT.values())
    #sumAllINT_CID = max(pksList_CID_DICT.values())
    #sumAllINT_UVPD = max(pksList_UVPD_DICT.values())
    print sumAllINT_CID, sumAllINT_UVPD
    print len(mz_CID)
    print len(mz_UVPD)
    counter = 0
    allDiff = list()
    oddsDICT = dict()


    bionDICT = dict()
    yionDICT = dict()
    fracINT_CID = dict()
    fracINT_UVPD = dict()
    notYorBIONS = dict()
    notYorBIONS_ODDS = dict()
    for k, cidINT in pksList_CID_DICT.items():
	try:
	    uvpdINT = pksList_UVPD_DICT[k]
	except KeyError:
	    uvpdINT = 0
	diff = cidINT - uvpdINT
	prob_cid = cidINT/sumAllINT_CID
	prob_uvpd = uvpdINT/sumAllINT_UVPD
	fracINT_CID[k] = prob_cid
	fracINT_UVPD[k] = prob_uvpd
	#prob_cid = cidINT
	#prob_uvpd = uvpdINT
	try: oddsRatio = (prob_cid)/(prob_uvpd)
	except ZeroDivisionError: oddsRatio = 0 # This shouldnt happen, but occurs when there is no peak for either cid or uvpd
	#print k, cidINT, uvpdINT, prob_cid, prob_uvpd
	if cidINT < 0.1 or uvpdINT < 0.1 : oddsRatio = 0
	diff = oddsRatio
	allDiff.append(diff)
	if k not in yions and k not in bions:
	    notYorBIONS[k] = [prob_cid*100, prob_uvpd*100, oddsRatio]
	    notYorBIONS_ODDS[k] = oddsRatio

	#if diff < 20: diff = 0
	oddsDICT[k] = diff


    asort = sorted(notYorBIONS_ODDS.iteritems(),key=itemgetter(1),reverse=True)[0:30]
    #notYorBIONS_keys = sorted(list(i[0] for i in asort))
    notYorBIONS_keys = list(i[0] for i in asort)
    yionCID_INT = list(fracINT_CID[k] for k in yions)
    yionUVPD_INT = list(fracINT_UVPD[k] for k in yions)
    bionCID_INT = list(fracINT_CID[k] for k in bions)
    bionUVPD_INT = list(fracINT_UVPD[k] for k in bions)

    print " Yion mz \t NORM_CID_INT \t NORM_UVPD_INT \t ODDS"
    for i,k in enumerate(yions):
	print "%d \t \t %f \t %f \t %f"%(k,yionCID_INT[i]*100,yionUVPD_INT[i]*100,(yionCID_INT[i]*100)/(yionUVPD_INT[i]*100))

    print " Bion mz \t NORM_CID_INT \t NORM_UVPD_INT \t ODDS"
    for i,k in enumerate(bions):
	print "%d \t \t %f \t %f \t %f"%(k,bionCID_INT[i]*100,bionUVPD_INT[i]*100,(bionCID_INT[i]*100)/(bionUVPD_INT[i]*100))

    print "Not YorBion mz \t NORM_CID_INT \t NORM_UVPD_INT \t ODDS"
    for i,k in enumerate(notYorBIONS_keys):
	print "%d \t \t %f \t %f \t %f"%(k,notYorBIONS[k][0],notYorBIONS[k][1],notYorBIONS[k][2])
    
    #print list(yionCID_INT)
    #print list(yionUVPD_INT)

    #sys.exit(1)   
    fig = p.figure(figsize=(20,20))
    ax = fig.add_subplot(1,1,1)
    #bins = 100
    #n,bins,patches = p.hist(list(yionCID_INT),100,color='g')
    #p.show()
    rects = ax.bar(oddsDICT.keys(),allDiff,log=True,bottom=10)
    allLabels = autolabel(rects,ax,oddsDICT)
    
    print allLabels
    p.xlabel('m/z ')
    p.ylabel('Odds Ratio - CID Peak vs UVPD Peak ')
    #ax.set_yscale('log')
    #p.ylim([10**-4,10**2])
    p.xlim([0,max(allLabels)])
    p.title('Difference in Peak list between CID and 5 pulses UVPD. Peptide : LVNELTEFAK with AMCA modification')
    fig_name = 'difference between CID and UVPD peaks_temp'
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'CIDandUVPD_LVNELTEFAK_difference.png'
    p.savefig(fig_name,format = "png", orientation='landscape')
    #p.show()
    
if __name__ == '__main__':
    #Global Definition; These must be defined and in this format else programme fails.
    bions = [428,542,671,784,885,1014,1161,1232]
    yions = [260,407,536,637,750,879,993,1092]
    t0 = time.clock()
    main()
    t1 = time.clock()
    print "Script - %s \t Completed in %s secs \t %s"%(sys.argv, t1-t0, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    
