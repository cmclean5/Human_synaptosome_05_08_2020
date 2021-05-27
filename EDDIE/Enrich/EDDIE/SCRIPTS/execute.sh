#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Enrich

SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Enrich/EDDIE/SCRIPTS

ANNODIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Enrich/Annotations

CLUSTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Enrich/PPI_PSP

PARMDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/Enrich/parameterFiles

RUNS=500
#RUNS=100

#FILTER=0
FILTER=1

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

anno=Annotations
Dir=PPI_PSP
parm=parameterFiles

if [ ! -d $anno ]; then
    mkdir $anno
fi

if [ ! -d $Dir ]; then
    mkdir $Dir
fi

if [ ! -d $parm ]; then
    mkdir $parm
fi

cp -v $ANNODIR/*.csv  $anno
cp -v $CLUSTDIR/*.csv $Dir
cp -v $PARMDIR/*.csv  $parm

EXE=run
SCRIPT1=submitClustEnrch.sh

cp -r $EXECDIR/$EXE .
cp -r $EXECDIR/$SCRIPT1 .
cp -r $SCRIPTDIR/runEnrichment.R .

# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
module load R 

chmod +x $EXE
chmod +x $SCRIPT1

time Rscript runEnrichment.R $RUNS $SUBDIR $FILTER

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

##cp -v enrichment.RData $OUTDIR
##cp -v enrichment.rds $OUTDIR
cp -v random_studies*.csv.gz $OUTDIR

echo "$0 done!"

exit 0
