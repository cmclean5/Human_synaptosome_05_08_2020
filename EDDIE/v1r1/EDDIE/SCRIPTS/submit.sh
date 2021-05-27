#!/bin/sh

echo "Running on Eddie..."

#load module R
#module load R

WORKINGDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/v1r1/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute.sh

START=1
END=550

name="Human_Synaptosome_20"

chmod +x $EXE

for i in `seq $START $END`
    do
     qsub -N $name -l h_rt=00:30:00 -l h_vmem=8G $EXE
     #qsub -N $name -l h_rt=01:00:00 -l h_vmem=8G $EXE
  done

echo "$0 done!"

exit 0
    
