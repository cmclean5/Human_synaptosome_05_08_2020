#!/bin/sh

echo "Running on Eddie..."

WORKINGDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/powerlaw/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute.sh

START=1
END=5 #1500

name="poweRlaw_20"

chmod +x $EXE

for i in `seq $START $END`
    do
     qsub -N $name -l h_rt=00:30:00 -pe sharedmem 4 -l h_vmem=8G $EXE
     #qsub -N $name -l h_rt=00:30:00 $EXE
  done

echo "$0 done!"

exit 0
    
