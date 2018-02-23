## ANALYSIS FOR H3N2 HA RBS DEEP MUTATIONAL SCANNING EXPERIMENTS
A Glu to Asp substitution at hemagglutinin (HA) residue 190 (E190D), which locates in the receptor-binding site and is part of the antigenic site B, emerged during the natural evolution of human H3N2 influenza virus in early 1990s. This study discovers that reverting Asp190 to Glu190 is lethal to recent human H3N2 strains. This study aims to examine the mechanistic basis of this evolutionary entrenchment. Three sets of deep mutational scanning experiment were performed based on A/Hong Kong/1/1968 (HK68), A/Moscow/10/1999 (Mos99), A/Wyoming/3/2003 (Wy03), A/Brisbane/10/2007 (Bris07), A/Victoria/361/2011 (Vic11):
  * Combinatorial mutant library of reversions from Bris07 to Wy03
  * Combinatorial mutant library of reversions from Bris07 rev6 mutant to Mos99
  * Triple mutant libraries (residues 225, 226, and 228):
    * HK68 WT
    * HK68 E190D
    * Vic11 WT
    * Vic11 D190E

### INPUT FILE
* All sequencing raw reads, which can be downloaded from NIH SRA database [PRJNA377321](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA377321), should be placed in fastq/ folder:
  * Combinatorial mutant library of reversions from Bris07 to Wy03
    * Input library: fastq/NW-4\_S4\_L001\_R1\_001.fastq and fastq/NW-4\_S4\_L001\_R2\_001.fastq
    * Passaged library (Replicate 1): fastq/NW-5\_S5\_L001\_R1\_001.fastq and fastq/NW-5\_S5\_L001\_R2\_001.fastq
    * Passaged library (Replicate 2): fastq/NW-6\_S6\_L001\_R1\_001.fastq and fastq/NW-6\_S6\_L001\_R2\_001.fastq
    * Passaged library (Replicate 3): fastq/NW-7\_S7\_L001\_R1\_001.fastq and fastq/NW-7\_S7\_L001\_R2\_001.fastq
  * Combinatorial mutant library of reversions from Bris07 rev6 mutant to Mos99
    * Input library: fastq/rev6mut9-Lib1\_S8\_L001\_R1\_001.fastq and fastq/rev6mut9-Lib1\_S8\_L001\_R2\_001.fastq
    * Passaged library (Replicate 1): fastq/rev6mut9-Lib2\_S9\_L001\_R1\_001.fastq and fastq/rev6mut9-Lib2\_S9\_L001\_R2\_001.fastq
    * Passaged library (Replicate 2): fastq/rev6mut9-Lib2\_S10\_L001\_R1\_001.fastq and fastq/rev6mut9-Lib2\_S10\_L001\_R2\_001.fastq
    * Passaged library (Replicate 3): fastq/rev6mut9-Lib2\_S11\_L001\_R1\_001.fastq and fastq/rev6mut9-Lib2\_S11\_L001\_R2\_001.fastq
  * Triple mutant libraries (residues 225, 226, and 228):
    * HK68 WT Input library: fastq/Tlib190-13\_S13\_L001\_R1\_001.fastq and fastq/Tlib190-13\_S13\_L001\_R2\_001.fastq
    * HK68 WT Passaged library (Replicate 1): fastq/Tlib190-1\_S1\_L001\_R1\_001.fastq and fastq/Tlib190-1\_S1\_L001\_R2\_001.fastq
    * HK68 WT Passaged library (Replicate 2): fastq/Tlib190-2\_S2\_L001\_R1\_001.fastq and fastq/Tlib190-2\_S2\_L001\_R2\_001.fastq
    * HK68 WT Passaged library (Replicate 3): fastq/Tlib190-3\_S3\_L001\_R1\_001.fastq and fastq/Tlib190-3\_S3\_L001\_R2\_001.fastq
    * HK68 E190D Input library: fastq/Tlib190-14\_S14\_L001\_R1\_001.fastq and fastq/Tlib190-14\_S14\_L001\_R2\_001.fastq
    * HK68 E190D Passaged library (Replicate 1): fastq/Tlib190-4\_S4\_L001\_R1\_001.fastq and fastq/Tlib190-4\_S4\_L001\_R2\_001.fastq
    * HK68 E190D Passaged library (Replicate 2): fastq/Tlib190-5\_S5\_L001\_R1\_001.fastq and fastq/Tlib190-5\_S5\_L001\_R2\_001.fastq
    * HK68 E190D Passaged library (Replicate 3): fastq/Tlib190-6\_S6\_L001\_R1\_001.fastq and fastq/Tlib190-6\_S6\_L001\_R2\_001.fastq
    * Vic11 WT Input library: fastq/Tlib190-15\_S15\_L001\_R1\_001.fastq and fastq/Tlib190-15\_S15\_L001\_R2\_001.fastq
    * Vic11 WT Passaged library (Replicate 1): fastq/Tlib190-7\_S7\_L001\_R1\_001.fastq and fastq/Tlib190-7\_S7\_L001\_R2\_001.fastq
    * Vic11 WT Passaged library (Replicate 2): fastq/Tlib190-8\_S8\_L001\_R1\_001.fastq and fastq/Tlib190-8\_S8\_L001\_R2\_001.fastq
    * Vic11 WT Passaged library (Replicate 3): fastq/Tlib190-9\_S9\_L001\_R1\_001.fastq and fastq/Tlib190-9\_S9\_L001\_R2\_001.fastq
    * Vic11 D190E Input library: fastq/Tlib190-16\_S16\_L001\_R1\_001.fastq and fastq/Tlib190-16\_S16\_L001\_R2\_001.fastq'
    * Vic11 D190E Passaged library (Replicate 1): fastq/Tlib190-10\_S10\_L001\_R1\_001.fastq and fastq/Tlib190-10\_S10\_L001\_R2\_001.fastq
    * Vic11 D190E Passaged library (Replicate 2): fastq/Tlib190-11\_S11\_L001\_R1\_001.fastq and fastq/Tlib190-11\_S11\_L001\_R2\_001.fastq
    * Vic11 D190E Passaged library (Replicate 3): fastq/Tlib190-12\_S12\_L001\_R1\_001.fastq and fastq/Tlib190-12\_S12\_L001\_R2\_001.fastq

### ANALYSIS PEPLINE
#### For Bris07 combinatorial mutant library
1. script/Bris07\_HARBSmapper1.py
  * Input file: fastq/NW-\*.fastq
  * Output file: result/Bris07\_MultiMutLib
2. script/Bris07\_HARBSmapper2.py
  * Input file: result/Bris07\_MultiMutLib.tsv
  * Output file: result/Bris07\_MutfitTable.tsv
3. script/Bris07\_PlotCompareRep.R
  * Input file: result/Bris07\_MultiMutLib.tsv
  * Output file: graph/Bris07\_RepCompare.png
4. script/Bris07\_PlotEnrichMut.R 
  * Input file: result/Bris07\_MutfitTable.tsv
  * Output file: graph/Bris07\_MutEnich\_D190E.png
#### For Bris07-rev6 combinatorial mutant library
1. script/Bris07r6\_HARBSmapper1.py
  * Input file: fastq/rev6mut9-\*.fastq
  * Output file: result/Bris07r6\_MultiMutLib.tsv
2. script/Bris07r6\_HARBSmapper2.py  
  * Input file: result/Bris07r6\_MultiMutLib.tsv
  * Output file: result/Bris07r6\_MutfitTable.tsv
3. script/Bris07r6\_PlotCompareRep.R
  * Input file: result/Bris07r6\_MultiMutLib.tsv
  * Output file: graph/Bris07r6\_RepCompare.png
4. script/Bris07r6\_PlotEnrichMut.R 
  * Input file: result/Bris07r6\_MutfitTable.tsv
  * Output file: graph/Bris07r6\_MutEnich.png
#### For HK68/Vic11 triple mutant libraries
1. script/Vic11\_HARBSmapper.py
  * Input file: fastq/Tlib190-\*.fastq
  * Output file: result/\*\_Tlib.count
2. script/Vic11\_HARBS\_Count2MutFitTable.py
  * Input file: result/\*\_Tlib.count
  * Output file: 
    * result/HK68WT\_MutFitTable.tsv
    * result/HK68E190D\_MutFitTable.tsv
    * result/Vic11WT\_MutFitTable.tsv
    * result/Vic11D190E\_MutFitTable.tsv
3. script/Vic11\_PlotCompareRep.R
  * Input file: result/\*\_Tlib.count
  * Output file: graph/Tlib\_RepCompare.png
4. script/Vic11\_PlotCompareLib.R 
  * Input file: result/\*\_Tlib.count
  * Output file: graph/Tlib\_RFcompare\_\*.png
