#!/bin/sh

echo "Running on Eddie..."

WORKINGDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/powerlaw/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute.sh

name="poweRlaw_17"

chmod +x $EXE

qsub -t 1-2:1 -N $name -l h_rt=01:00:00 -pe sharedmem 4 $EXE

echo "$0 done!"

exit 0
    
