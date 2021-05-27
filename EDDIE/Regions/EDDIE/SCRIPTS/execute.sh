#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Regions
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Regions/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Regions/datasets
GRAPHDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Regions/Graphs

SEED=$SUBDIR
ALG=1
PERMS=100

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

#--- load R scripts
cp -r $SCRIPTDIR/createRandomStudySet.R .

#--- load graphs
#cp -r $GRAPHDIR/full_FilteredPSD_ppi.gml .
cp -r $GRAPHDIR/PPI_comp.gml .

#--- load datasets
cp -r $DATADIR/*.csv .

# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
module load R 

Rscript createRandomStudySet.R  $SEED $ALG $PERMS

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v results_*.csv $OUTDIR

echo "$0 done!"

exit 0
