#!/usr/bin/python
import os
import sys
import string
import operator
import numpy as np
from Bio import SeqIO
from itertools import imap
from collections import Counter

def formatting(infile,outfile,inputcutoff,muts):
  print "writing: %s" % outfile
  infile  = open(infile,'r')
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','mutclass','InputCount','fitness'])+"\n")
  for line in infile.xreadlines(): 
    if 'mut' in line: continue
    else:
      line = line.rstrip().rsplit("\t")
      mut  = line[0]
      mutclass = int(line[1])
      iptcount = int(line[2])
      Rep1Fit  = float(line[9])
      Rep2Fit  = float(line[10])
      Rep3Fit  = float(line[11])
      if iptcount >= inputcutoff: 
        if mut=='WT' or mut.count('-')+1==sum([mut.count(m) for m in muts]):
          fit = np.mean([Rep1Fit,Rep2Fit,Rep3Fit])
          outfile.write("\t".join(map(str,[mut, mutclass, iptcount, fit]))+"\n")
  infile.close()
  outfile.close()

def main():
  inputcutoff = 7
  muts    = ['N145K','G186V','N189S','D190E','F193S','P194L','S219F','N225D','P227S']
  infile  = 'result/Bris07_MultiMutLib.tsv'
  outfile = 'result/Bris07_MutfitTable.tsv'
  formatting(infile, outfile, inputcutoff, muts)

if __name__ == "__main__":
  main()
