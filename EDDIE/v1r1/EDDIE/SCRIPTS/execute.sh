#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1/datasets

#---algorithms, Rclustering.R
#alg     <- vector(length=16)
#alg[1]  <- "louvain"
#alg[2]  <- "infomap"
#alg[3]  <- "fc"
#alg[4]  <- "lec"
#alg[5]  <- "sgG1"
#alg[6]  <- "sgG2"
#alg[7]  <- "sgG5"
#alg[8]  <- "wt"
#alg[9]  <- "spectral"
#alg[10] <- "louvain2"
#alg[11] <- "infomap2"
#alg[12] <- "fc2"
#alg[13] <- "lec2"
#alg[14] <- "sgG12"
#alg[15] <- "wt2"
#alg[16] <- "Spectral2"

ALG=10
TYPE=2
MASK=20
CNMIN=-1
CNMAX=-1

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

cp -r $SCRIPTDIR/Rclustering.R .
cp -r $DATADIR/PPI_Presynaptic.gml .
#cp -r $DATADIR/PPI_PSP.gml .
#cp -r $DATADIR/PPI_PSP_consensus.gml .

#cp -r $DATADIR/example1.gml .
#cp -r $DATADIR/example2.gml .
#cp -r $DATADIR/karate.gml .
#cp -r $DATADIR/PPI_Network.gml .

#Database paper (28/01/2020)
#cp -r $DATADIR/full_FilteredPSD_ppi.gml .

# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
module load R 
#module load igmm/apps/R/3.2.2 #this should work for igraph & poweRlaw

time Rscript Rclustering.R $ALG $TYPE $MASK $CNMIN $CNMAX

#if [ $ALG == 9 ]; then
#
# mkdir OUT
# EXE=run
# 
# cp -r $EXECDIR/$EXE .
# cp -r $SCRIPTDIR/Rformat.R .
# time ./$EXE -file "edgelist.csv" -cols 2 -a 3
#
# Rscript Rformat.R
#
#fi

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.txt $OUTDIR

echo "$0 done!"

exit 0
