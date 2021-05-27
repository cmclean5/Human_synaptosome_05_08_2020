#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

#ARRAYID=$SGE_TASK_ID
#SUBDIR=$JOB_ID.$ARRAYID
SUBDIR=$JOB_ID

echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/powerlaw
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/powerlaw/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/powerlaw

DATA=datasets.tar.gz
#SEED=$(($JOB_ID+$ARRAYID))
SEED=$JOB_ID
NSIM=1000
NTHREADS=4

NSTUDY=2

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

cp -r $DATADIR/$DATA .
cp -r $SCRIPTDIR/Centrality.R .

# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
#module load R #issue with this??? 11/8/17
module load igmm/apps/R/3.2.2 #this should work for igraph & poweRlaw
#module load igmm/apps/R/3.6.0

tar -xzvf $DATA

Rscript Centrality.R $SEED $NSIM $NTHREADS $NSTUDY

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.csv        $OUTDIR

echo "$0 done!"

exit 0
