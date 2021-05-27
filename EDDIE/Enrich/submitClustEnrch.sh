#!/bin/sh

echo "Run Cluster Enrichment Studies"

#Read in offset for random number seed 
seed=$1

#echo "add to seed ((${seed}))"

#Working Dir
WORKINGDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Main Package location
#MAINDIR=$(dirname "$WORKINGDIR")

#Graph Location
NTWRKDIR=$WORKINGDIR/Graphs

#Community Location
COMDIR=$WORKINGDIR/PPI_PSP

#Annotation Location
ANNODIR=$WORKINGDIR/Annotations

#---Read in Parameters
#---Graphs
graphFILE=$WORKINGDIR/parameterFiles/graphs.csv
declare -a graphON SUBDIR

IFS=$'\n\t' ARRAY=($(<$graphFILE))

len=${#ARRAY[@]}
k=0
for ((i=0; i<$len; i=i+2)); do
    graphON[$k]=${ARRAY[$i]}
    #echo ${graphON[$k]}
    k=$(($k+1))
done    

k=0
for ((i=1; i<$len; i=i+2)); do
    SUBDIR[$k]=${ARRAY[$i]}
    #echo ${SUBDIR[$k]}
    k=$(($k+1))
done    
#---

#---Annotations
annoFILE=$WORKINGDIR/parameterFiles/annotations.csv
declare -a annoON annoFILES annoTITLES

IFS=$'\n\t' ARRAY=($(<$annoFILE))

len=${#ARRAY[@]}
k=0
for ((i=0; i<$len; i=i+3)); do
    annoON[$k]=${ARRAY[$i]}
    #echo ${annoON[$k]}
    k=$(($k+1))
done  

k=0
for ((i=1; i<$len; i=i+3)); do
    annoFILES[$k]=${ARRAY[$i]}
    #echo ${annoFILES[$k]}
    k=$(($k+1))
done    

k=0
for ((i=2; i<$len; i=i+3)); do
    annoTITLES[$k]=${ARRAY[$i]}
    #echo ${annoTITLES[$k]}
    k=$(($k+1))
done    
#---

#---Clustering Alg.
comFILE=$WORKINGDIR/parameterFiles/clusteringAlg.csv
declare -a comON comFILES comTITLES

IFS=$'\n\t' ARRAY=($(<$comFILE))

#len=${#ARRAY[@]}    
#for ((i=0; i<$len; i++)); do    
#    echo "$i \t" ${ARRAY[$i]}
#done    

len=${#ARRAY[@]}
k=0
for ((i=0; i<$len; i=i+3)); do
    comON[$k]=${ARRAY[$i]}
    #echo ${comON[$k]}
    k=$(($k+1))
done  


k=0
for ((i=1; i<$len; i=i+3)); do
    comFILES[$k]=${ARRAY[$i]}
    #echo ${comFILES[$k]}
    k=$(($k+1))
done    

k=0
for ((i=2; i<$len; i=i+3)); do
    comTITLES[$k]=${ARRAY[$i]}
    #echo ${comTITLES[$k]}
    k=$(($k+1))
done    
#---
#--- END READ-IN

#---Check OUT and RESULTS directories exist
#--- RESULTS DIR
RESDIR=$WORKINGDIR/RESULTS
if [ ! -d $RESDIR ]; then
    mkdir $RESDIR
fi

#--- OUT DIR
tmpDIR=OUT
TEMPDIR=$WORKINGDIR/$tmpDIR
if [ ! -d $TEMPDIR ]; then
    mkdir $TEMPDIR
fi
#---

EXE=run
chmod +x $EXE

#---SET LOOP COUNTERS
gSTART=0
gEND=${#SUBDIR[@]}
gEND=$(($gEND-1))

START=0
END=${#annoTITLES[@]}
END=$(($END-1))

kSTART=0
kEND=${#comTITLES[@]}
kEND=$(($kEND-1))
#--- END

MINov1=2

#---loop over graphs
for g in `seq $gSTART $gEND`
do

    if [ ${graphON[$g]} -eq "1" ]; then 
    
	GRAPHENRICHOUT=$WORKINGDIR/RESULTS/${SUBDIR[$g]}
	if [ ! -d $GRAPHENRICHOUT ]; then
	    mkdir $GRAPHENRICHOUT
	fi
    
	#---loop over clustering algorithms
	for k in `seq $kSTART $kEND`
	do
	    #echo "Running over clustering algorithm ${comTITLES[$k]}:"

	    #---Create output directory for cluster enrichment results 
	    OUTDIR=$WORKINGDIR/RESULTS/${SUBDIR[$g]}/${comTITLES[$k]}
	    if [ ! -d $OUTDIR ]; then
		mkdir $OUTDIR
	    fi
	    
	    #---loop over annotation set 
	    for a in `seq $START $END`
	    do    
		#echo "Annotation set: ${annoTITLES[$a]}"
		

		if [[ ${comON[$k]} -eq "1" && ${annoON[$a]} -ge "1" ]]; then
		    
		    echo "Running over clustering algorithm ${comTITLES[$k]}:"
		    echo "Annotation set: ${annoTITLES[$a]}"

		    #--- Clean temp dir, i.e. "OUT/" 
		    #rm $WORKINGDIR/$tmpDIR/*   
		    rm $TEMPDIR/*.csv 
		    
		    sleep 1

		    #---run cluster enrichment for annotation set
		    time ./$EXE -Comfile $COMDIR/${comFILES[$k]} -Annofile $ANNODIR/${annoFILES[$a]} -ext ${annoTITLES[$a]} -o $tmpDIR -opt 5 -setFDR BY -setSEED ${seed} -printID -noPerm -onesided
		    
		    #---Copy cluster enrichment results to output directory
		    cp -v $TEMPDIR/*.csv $OUTDIR

		    #--- Clean temp dir, i.e. "OUT/" 
		    rm $TEMPDIR/*.csv 
		    
		fi
		
	    done
	done
    fi
done

echo "$0 done!"

exit 0
    
