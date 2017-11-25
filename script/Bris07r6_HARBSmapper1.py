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

def ProcessMultilib(R1file):
  print "Reading %s" % R1file
  R2file = R1file.replace('_R1_','_R2_')
  Primerlength = 23
  R1frameshift = 1
  R2frameshift = 2
  roilength    = 255-R1frameshift-R2frameshift
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  variants = [] 
  for R1record in R1records:
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi = R1seq[Primerlength+R1frameshift:Primerlength+R1frameshift+roilength]
    R2roi = R2seq[Primerlength+R2frameshift:Primerlength+R2frameshift+roilength]
    if 'N' in R1roi or 'N' in R2roi: continue
    R1pep = translation(R1roi)
    R2pep = translation(rc(R2roi))
    if R1pep == R2pep:
      variants.append(R1pep)
  return Counter(variants)

def Bris07mut2ID(mut,refseq):
  shift = 144
  haplo = []
  assert(len(mut)==len(refseq))
  for n in range(len(mut)):
    pos = n+shift
    if refseq[n]!=mut[n]:
       haplo.append(refseq[n]+str(pos)+mut[n])
  return '-'.join(haplo)

def Output(InputDict, R1Rep1Dict, R1Rep2Dict, R1Rep3Dict, outfile, refseq):
  print "Compiling results into %s" % outfile
  outfile = open(outfile,'w')
  muts = list(set(InputDict.keys()+R1Rep1Dict.keys()+R1Rep2Dict.keys()+R1Rep3Dict.keys()))
  WT_R1Rep1EnrichRatio = float(R1Rep1Dict[refseq]+1)/float(InputDict[refseq]+1)
  WT_R1Rep2EnrichRatio = float(R1Rep2Dict[refseq]+1)/float(InputDict[refseq]+1)
  WT_R1Rep3EnrichRatio = float(R1Rep3Dict[refseq]+1)/float(InputDict[refseq]+1)
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
    mutclass    = 0 if refseq == mut else hamming(refseq,mut)
    ID          = 'WT' if refseq == mut else Bris07mut2ID(mut,refseq)
    outfile.write("\t".join(map(str,[ID,mutclass,InputDict[mut],
                                     R1Rep1Dict[mut],R1Rep2Dict[mut],R1Rep3Dict[mut],
                                     R1Rep1EnrichRatio,R1Rep2EnrichRatio,R1Rep3EnrichRatio,
                                     R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness]))+"\n")
  outfile.close()

def main():
  refseq  = 'NNSFFSRLNWLTHLKFKYPALNVTMPNNEKFDKLYIWGVHHPSTDNEQISLYAQASGRITVSTKRSQQTVIPNIGSRPRVRGIS'
  outfile = 'result/Bris07r6_MultiMutLib.tsv'
  Bris07_InputDict   = ProcessMultilib('fastq/rev6mut9-Lib1_S8_L001_R1_001.fastq')
  Bris07_R1Rep1Dict  = ProcessMultilib('fastq/rev6mut9-Lib2_S9_L001_R1_001.fastq')
  Bris07_R1Rep2Dict  = ProcessMultilib('fastq/rev6mut9-Lib3_S10_L001_R1_001.fastq')
  Bris07_R1Rep3Dict  = ProcessMultilib('fastq/rev6mut9-Lib4_S11_L001_R1_001.fastq')
  Output(Bris07_InputDict,Bris07_R1Rep1Dict,Bris07_R1Rep2Dict,Bris07_R1Rep3Dict,outfile,refseq)

if __name__ == "__main__":
  main()
