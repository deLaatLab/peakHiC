#!/usr/bin/env bash
#full_run - a wrapper script that calls peakHiC; change file to run

#full_run.sh: change parameters within script

##optional: set up environment##
#conda env create -f peakhic_environmnent.yml
#R> install.packages("isotone")
#optional: tmux new -s peakhic

##configure and rename conf_hg38_GM12878_example_data.yml file to create appropriate peakHiCObj##
#R> source("../R/peakHiC_functions.R")
#R> peakHiCConf("../conf_X.yml")

conda activate peakHiC
cd /path/to/peakHiCFolder/RUN #~/localdev/peakHiC/RUN

##run peakHiC##
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
PEAKOBJ=/path/to/peakHiCObj #~/data/peakHiCObj/peakHiCObj.rds

for CHR in $CHROMS; do 
echo $PEAKOBJ
echo $CHR
Rscript createPartitionV4CsByChr.R -chr $CHR -peakHiCObj $PEAKOBJ
Rscript callPartitionPeaksbyChr.R -chr $CHR -peakHiCObj $PEAKOBJ
Rscript processLoops.R -chr $CHR -peakHiCObj $PEAKOBJ
done


