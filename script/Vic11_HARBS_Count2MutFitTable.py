#!/usr/bin/python
import os
import sys
import string
import operator
import numpy as np
from Bio import SeqIO
from itertools import imap
from collections import Counter

def mut2ID(mut, WT):
  WT_resi = WT.rsplit('-')
  if mut == 'WT':
    return ''.join([resi[-1] for resi in WT_resi])
  for resi in WT_resi:
    pos = resi[0:-1]
    aa  = resi[-1]
    if pos not in mut:
      mut += '-'+resi
  mut = sorted(mut.rsplit('-'),key=lambda x:x[0:-1])
  ID  = ''.join([resi[-1] for resi in mut])
  return ID

def ID2mut(ID,WT): 
    mut   = []
    WT_resi = WT.rsplit('-')
    WT_pos  = [resi[0:-1] for resi in WT_resi]
    WT_ID   = ''.join([resi[-1] for resi in WT_resi])
    if ID == WT_ID: return 'WT'
    for pos, aa in zip(WT_pos,ID):
      resi = str(pos)+aa
      if resi in WT_resi: continue
      else: mut.append(resi)
    mut = '-'.join(mut)
    return mut

def CutoffMutFitTable(infile, outfile, WT, inputcutoff):
  infile    = open(infile,'r')
  outfile   = open(outfile,'w')
  outfile.write('Mut'+"\t"+'Fitness'+"\n")
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    ID  = line[0]
    inputcount = int(line[2])
    Fitness  = str(np.mean([float(line[9]),float(line[10]),float(line[11])]))
    if '-' in ID: continue
    if inputcount < inputcutoff: continue
    mut = ID2mut(ID, WT)
    #outfile.write(mut+"\t"+ID+"\t"+Fitness+"\n") 
    outfile.write(mut+"\t"+Fitness+"\n") 
  infile.close()
  outfile.close()

def hashin(infile,WT):
  infile = open(infile,'r')
  fithash = {}
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    fithash[mut2ID(line[0],WT)] = line[1]
  infile.close()
  return fithash

def main():
  inputcutoff=20
  CutoffMutFitTable('result/HK68WT_Tlib.count', 'result/HK68WT_MutFitTable.tsv','225G-226L-228S', inputcutoff)
  CutoffMutFitTable('result/HK68E190D_Tlib.count', 'result/HK68E190D_MutFitTable.tsv','225G-226L-228S', inputcutoff)
  CutoffMutFitTable('result/Vic11WT_Tlib.count', 'result/Vic11WT_MutFitTable.tsv','225N-226I-228S', inputcutoff)
  CutoffMutFitTable('result/Vic11D190E_Tlib.count', 'result/Vic11D190E_MutFitTable.tsv','225N-226I-228S', inputcutoff)
  
if __name__ == "__main__":
  main()
