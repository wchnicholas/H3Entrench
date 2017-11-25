#!/usr/bin/python
import os
import sys
import string
import operator
from Bio import SeqIO
from itertools import imap
from collections import Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def ProcessTlib(R1file, strain, resi190, WT):
  R2file = R1file.replace('_R1_','_R2_')
  R1fixpos    = 26
  R2fixpos    = 141
  R1RandPos = 131
  R2RandPos = 27
  Length  = 12
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  muts = []
  for R1record in R1records:
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1fix  = R1seq[R1fixpos:R1fixpos+3]
    R2fix  = R2seq[R2fixpos:R2fixpos+3]
    R1roi  = R1seq[R1RandPos:R1RandPos+Length]
    R2roi  = R2seq[R2RandPos:R2RandPos+Length]
    R1flank5 = R1seq[R1RandPos-3:R1RandPos]
    R1flank3 = R1seq[R1RandPos+Length:R1RandPos+Length+3]
    R2flank5 = R2seq[R2RandPos-3:R2RandPos]
    R2flank3 = R2seq[R2RandPos+Length:R2RandPos+Length+3]
    if R1roi != rc(R2roi) or R1fix != rc(R2fix): continue
    if translation(R1fix) != resi190: continue
    if strain == 'HK68':
      if R1flank5 != 'AGG' or R1flank3 != 'AGA' or R2flank5 != 'TCT' or R2flank3 != 'CCT': continue
    elif strain == 'Vic11':
      if R1flank5 != 'AGG' or R1flank3 != 'AGA' or R2flank5 != 'TCT' or R2flank3 != 'CCT': continue
    else: 
      print "Strain name not found: %s" % strain; sys.exit()
    if 'N' in R1roi or 'N' in R2roi: continue
    #print translation(R1fix), translation(rc(R2fix)), translation(R1roi), translation(rc(R2roi))
    roi  = R1roi[0:6]+R1roi[9:12]
    Tmut = translation(roi)
    if Tmut == WT:
      muts.append('WT')
    else: 
      muts.append(Tmut)
  R1records.close()
  R2records.close()
  print "Finish reading %s and %s" % (R1file, R2file)
  return Counter(muts)

def Output(InputDict, R1Rep1Dict, R1Rep2Dict, R1Rep3Dict, outfile, WT):
  print "Compiling results into %s" % outfile
  outfile = open(outfile,'w')
  muts = list(set(InputDict.keys()+R1Rep1Dict.keys()+R1Rep2Dict.keys()+R1Rep3Dict.keys()))
  WT_R1Rep1EnrichRatio = float(R1Rep1Dict['WT']+1)/float(InputDict['WT']+1)
  WT_R1Rep2EnrichRatio = float(R1Rep2Dict['WT']+1)/float(InputDict['WT']+1)
  WT_R1Rep3EnrichRatio = float(R1Rep3Dict['WT']+1)/float(InputDict['WT']+1)
  outfile.write("\t".join(['mut','mutclass','InputCount',
                           'R1Rep1Count','R1Rep2Count','R1Rep3Count',
                           'R1Rep1EnrichRatio','R1Rep2EnrichRatio','R1Rep3EnrichRatio',
                           'R1Rep1Fitness','R1Rep2Fitness','R1Rep3Fitness'])+"\n")
  for mut in muts:
    R1Rep1EnrichRatio = float(R1Rep1Dict[mut]+1)/float(InputDict[mut]+1)
    R1Rep2EnrichRatio = float(R1Rep2Dict[mut]+1)/float(InputDict[mut]+1)
    R1Rep3EnrichRatio = float(R1Rep3Dict[mut]+1)/float(InputDict[mut]+1)
    R1Rep1Fitness     = float(R1Rep1EnrichRatio)/float(WT_R1Rep1EnrichRatio)
    R1Rep2Fitness     = float(R1Rep2EnrichRatio)/float(WT_R1Rep2EnrichRatio)
    R1Rep3Fitness     = float(R1Rep3EnrichRatio)/float(WT_R1Rep3EnrichRatio)
    mutclass    = 0 if 'WT' == mut or 'WT-' in mut else hamming(WT,mut)
    mut         = WT if 'WT' == mut else mut
    outfile.write("\t".join(map(str,[mut,mutclass,InputDict[mut],
                                     R1Rep1Dict[mut],R1Rep2Dict[mut],R1Rep3Dict[mut],
                                     R1Rep1EnrichRatio,R1Rep2EnrichRatio,R1Rep3EnrichRatio,
                                     R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness]))+"\n")
  outfile.close()

def main():
  HK68WT_InputDict   = ProcessTlib('fastq/Tlib190-13_S13_L001_R1_001.fastq', 'HK68','E','GLS')
  HK68WT_R1Rep1Dict  = ProcessTlib('fastq/Tlib190-1_S1_L001_R1_001.fastq', 'HK68','E','GLS')
  HK68WT_R1Rep2Dict  = ProcessTlib('fastq/Tlib190-2_S2_L001_R1_001.fastq', 'HK68','E','GLS')
  HK68WT_R1Rep3Dict  = ProcessTlib('fastq/Tlib190-3_S3_L001_R1_001.fastq', 'HK68','E','GLS')
  Output(HK68WT_InputDict,HK68WT_R1Rep1Dict,HK68WT_R1Rep2Dict,HK68WT_R1Rep3Dict,'result/HK68WT_Tlib.count','GLS') 
  HK68E190D_InputDict   = ProcessTlib('fastq/Tlib190-14_S14_L001_R1_001.fastq', 'HK68','D','GLS')
  HK68E190D_R1Rep1Dict  = ProcessTlib('fastq/Tlib190-4_S4_L001_R1_001.fastq', 'HK68','D','GLS')
  HK68E190D_R1Rep2Dict  = ProcessTlib('fastq/Tlib190-5_S5_L001_R1_001.fastq', 'HK68','D','GLS')
  HK68E190D_R1Rep3Dict  = ProcessTlib('fastq/Tlib190-6_S6_L001_R1_001.fastq', 'HK68','D','GLS')
  Output(HK68E190D_InputDict,HK68E190D_R1Rep1Dict,HK68E190D_R1Rep2Dict,HK68E190D_R1Rep3Dict,'result/HK68E190D_Tlib.count','GLS') 
  Vic11WT_InputDict  = ProcessTlib('fastq/Tlib190-15_S15_L001_R1_001.fastq', 'Vic11','D','NIS')
  Vic11WT_R1Rep1Dict = ProcessTlib('fastq/Tlib190-7_S7_L001_R1_001.fastq', 'Vic11','D','NIS')
  Vic11WT_R1Rep2Dict = ProcessTlib('fastq/Tlib190-8_S8_L001_R1_001.fastq', 'Vic11','D','NIS')
  Vic11WT_R1Rep3Dict = ProcessTlib('fastq/Tlib190-9_S9_L001_R1_001.fastq', 'Vic11','D','NIS')
  Output(Vic11WT_InputDict,Vic11WT_R1Rep1Dict,Vic11WT_R1Rep2Dict,Vic11WT_R1Rep3Dict, 'result/Vic11WT_Tlib.count','NIS') 
  Vic11D190E_InputDict  = ProcessTlib('fastq/Tlib190-16_S16_L001_R1_001.fastq', 'Vic11','E','NIS')
  Vic11D190E_R1Rep1Dict = ProcessTlib('fastq/Tlib190-10_S10_L001_R1_001.fastq', 'Vic11','E','NIS')
  Vic11D190E_R1Rep2Dict = ProcessTlib('fastq/Tlib190-11_S11_L001_R1_001.fastq', 'Vic11','E','NIS')
  Vic11D190E_R1Rep3Dict = ProcessTlib('fastq/Tlib190-12_S12_L001_R1_001.fastq', 'Vic11','E','NIS')
  Output(Vic11D190E_InputDict,Vic11D190E_R1Rep1Dict,Vic11D190E_R1Rep2Dict,Vic11D190E_R1Rep3Dict,'result/Vic11D190E_Tlib.count','NIS') 

if __name__ == "__main__":
  main()
