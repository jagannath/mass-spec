#! /home/anaconda/python2.7

"""
This script is aimed to parse the ms2 file; This file is obtained from the mzXML
file for only the ms2 spectra. The file is organized such that the first
precursor mass is used for CID and the next UVPD.  
"""

import sys
import os
import collections
import diffSpectra

class MS2(object):
    def __init__(self,ms2fname):
        self.ms2fname = ms2fname
        self.header = str()
        self.ms2DICT = collections.defaultdict(list)
    def stripHeader(self,lines):
        counter = 0
        for line in lines:
            if line.startswith('H'):
                counter+=1
                self.header+=line
            else:
                return self.header,lines[counter:]
    def breakSpectra(self):
        with open(self.ms2fname) as f:
            lines = f.readlines()
        self.header,lines = self.stripHeader(lines)
        precursor, spec  = str(), list()
        #ms2DICT[precursorMass]=[(RTime,#),(BPI,#),(BPM,#),(TIC,#),[(mz,intensity)...]]
        for line in lines:
            if line.startswith('S'):
                if (precursor) and (spec): 
                    specValues = [(k,v) for k,v in specDetails.items()]
                    specValues.append(spec)
                    self.ms2DICT[precursor].append(specValues)
                    print len(self.ms2DICT[precursor]),specValues[0]
                    if len(self.ms2DICT[precursor])==2: yield (self.ms2DICT)
                else: 
                    precursor = str()
                    spec = list()
                    inSpec = False
                    specDetails = dict()
            elif line.startswith('I'):
                inSpec = True
                id,nbr = line[1:].strip().split('\t')
                specDetails[id]=nbr
            else:
                precursor = specDetails['BPM']
                mz,intensity = map(float,line.split(' '))
                spec.append((mz,intensity))

class MGF2(object):
    def __init__(self,mgf2fname):
        self.mgf2fname = mgf2fname
        self.ms2DICT = collections.defaultdict(list)
    def breakSpectra(self):
        with open(self.mgf2fname) as f:
            lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if line.startswith('BEGIN IONS'):
                inspec = True
                specDetails = dict()
                spec = list()
            elif line.startswith('END IONS'):
                inspec = False
                specValues = [(k,v) for k,v in specDetails.items()]
                specValues.append(spec)
                self.ms2DICT[precursor].append(specValues)
#                if len(self.ms2DICT[precursor])==2: yield (self.ms2DICT)
            elif '=' in line:
                if line.startswith('TITLE'):
                    lst = [line.split(' ')[i] for i in [0,2,3]]
                    for l in lst:
                        id,val = l.split('=')
                        specDetails[id]=val
                else:
                    id,val = line.split('=')
                    specDetails[id]=val
            else:
                precursor = specDetails['PEPMASS']
                mz,intensity = map(float,line.split(' '))
                spec.append((mz,intensity))



def makeDICT(twoSpec):
    # Order is assumed to be cid followed by uvpd Maybe later put some condition
    
    try:
        cid,uvpd = map(lambda x: dict(iter(x)),[twoSpec[0][-1],twoSpec[1][-1]])
    except IndexError:
        print twoSpec
        sys.exit(1)
    return cid,uvpd



def main():
    DIR ='/project/marcotte/jagannath/projectfiles/proteomics_single_molecule/jennyBroadbelt/scott/data/BSA_AMCA/'
    msfile = DIR + 'BSA_AMCA_CIDUVPD.mgf'
    #ms2 = MS2(msfile)
    ms2 = MGF2(msfile)
    ms2.breakSpectra()
    DICT = ms2.ms2DICT

    dirName = "BSA_AMCA"
    counter = 0
    for k,v in DICT.items():
        print "Processing %d Precursor Mass : %s ..."%(counter,k)
        if not len(v)==2: 
            print "Skipping this ..."
            continue
        counter+=1
        fname = "PEPMASS=%s"%(k)
        cid_DICT,uvpd_DICT = makeDICT(v)
        cid,uvpd = map(diffSpectra.MassSpecDICT,[cid_DICT,uvpd_DICT])
        cid,uvpd = cid.roundoff(0),uvpd.roundoff(0)
        combinedSpec = diffSpectra.ManageData(cid,uvpd)
        df_summary,df_final = combinedSpec()
        if df_final.empty:
            out = diffSpectra.Output(df_summary,dirName)
            out.writeCSV(fname)
            out.spectra(fname=fname+'.png')
            print "In summary"
        else:
            out = diffSpectra.Output(df_final,dirName)
            out.writeCSV(fname,bion=True)
            out.spectra(fname=fname+'.png',bion=True)
            print df_final
            print "In Final"



main()







