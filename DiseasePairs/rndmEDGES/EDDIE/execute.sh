#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/localization/datasets

SEED=$SUBDIR
PERMS=10

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

#--- load R scripts
cp -r $SCRIPTDIR/DiseaseLoc.R .
cp -r $SCRIPTDIR/annotationTYPES.R .

#--- load graphs
#cp -r $DATADIR/PPI_Presynaptic.gml .
#cp -r $DATADIR/PPI_PSP.gml .
cp -r $DATADIR/PPI_PSP_consensus.gml .
#cp -r $DATADIR/PPI_PSP_consensus2.gml .


# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
module load R 
#module load igmm/apps/R/3.2.2

Rscript DiseaseLoc.R $SEED $PERMS

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.csv $OUTDIR

echo "$0 done!"

exit 0
