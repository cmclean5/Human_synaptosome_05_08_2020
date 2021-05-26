#!/bin/bash

echo "Run Overlap Cluster/Regions Enrichment Studies"

#Working Dir
WORKINGDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#Main Package location
MAINDIR=$(dirname "$WORKINGDIR")

#Bridging Region Location
REGDIR=$MAINDIR/Consensus/

#Community Location
COMDIR=$MAINDIR/Clustering

#Annotation Location
ANNODIR=$MAINDIR/Annotations

#---Read in Parameters
#---Graphs
graphFILE=$MAINDIR/parameterFiles/graphs.csv
declare -a graphON SUBDIR

IFS=$'\n\t' ARRAY=($(<$graphFILE))

len=${#ARRAY[@]}
k=0
for ((i=0; i<$len; i=i+2)); do
    graphON[$k]=${ARRAY[$i]}
    echo ${graphON[$k]}
    k=$(($k+1))
done    

k=0
for ((i=1; i<$len; i=i+2)); do
    SUBDIR[$k]=${ARRAY[$i]}
    echo ${SUBDIR[$k]}
    k=$(($k+1))
done    
#---


#---Annotations
annoFILE=$MAINDIR/parameterFiles/annotations.csv
declare -a annoON annoFILES annoTITLES

IFS=$'\n\t' ARRAY=($(<$annoFILE))

len=${#ARRAY[@]}
k=0
KEYind=0
for ((i=0; i<$len; i=i+3)); do
    annoON[$k]=${ARRAY[$i]}
    echo ${annoON[$k]}
    if [[ ${annoON[$k]} -eq "2" ]]; then
	KEYind=$k
    fi    
    k=$(($k+1))
done  

k=0
for ((i=1; i<$len; i=i+3)); do
    annoFILES[$k]=${ARRAY[$i]}
    echo ${annoFILES[$k]}
    k=$(($k+1))
done    

k=0
for ((i=2; i<$len; i=i+3)); do
    annoTITLES[$k]=${ARRAY[$i]}
    echo ${annoTITLES[$k]}
    k=$(($k+1))
done    
#---

#---Clustering Alg.
comFILE=$MAINDIR/parameterFiles/clusteringAlg.csv
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
    echo ${comON[$k]}
    k=$(($k+1))
done  


k=0
for ((i=1; i<$len; i=i+3)); do
    comFILES[$k]=${ARRAY[$i]}
    echo ${comFILES[$k]}
    k=$(($k+1))
done    

k=0
for ((i=2; i<$len; i=i+3)); do
    comTITLES[$k]=${ARRAY[$i]}
    echo ${comTITLES[$k]}
    k=$(($k+1))
done    
#---

#---Regions
regFILE=$MAINDIR/parameterFiles/RegionsAlg.csv
declare -a regON regFILES regTITLES

IFS=$'\n\t' ARRAY=($(<$regFILE))

len=${#ARRAY[@]}
k=0
for ((i=0; i<$len; i=i+3)); do
    regON[$k]=${ARRAY[$i]}
    echo ${regON[$k]}
    k=$(($k+1))
done 

k=0
for ((i=1; i<$len; i=i+3)); do
    regFILES[$k]=${ARRAY[$i]}
    echo ${regFILES[$k]}
    k=$(($k+1))
done    

k=0
for ((i=2; i<$len; i=i+3)); do
    regTITLES[$k]=${ARRAY[$i]}
    echo ${regTITLES[$k]}
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
OUTDIR=$WORKINGDIR/$tmpDIR
if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
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

rSTART=0
rEND=${#regTITLES[@]}
rEND=$(($rEND-1))
#--- END

#--- SET MIN OVERLAP and geneset size
MINOV1=5
MINOV2=5
MINOV3=5

#--- ADD OFFSET test enrichment of non-Bridging proteins
OFFSET="y"
#OFFSET="n"

#---loop over graphs
for g in `seq $gSTART $gEND`
do

    if [ ${graphON[$g]} -eq "1" ]; then

	GRAPHENRICHOUT=$WORKINGDIR/RESULTS/${SUBDIR[$g]}
	if [ ! -d $GRAPHENRICHOUT ]; then
	    mkdir $GRAPHENRICHOUT
	fi

	REGIONENRICHOUT=$WORKINGDIR/RESULTS/${SUBDIR[$g]}/REGIONS
	if [ ! -d $REGIONENRICHOUT ]; then
	    mkdir $REGIONENRICHOUT
	fi
    
    
	#---loop over clustering algorithms
	for k in `seq $kSTART $kEND`
	do
	    echo "Running over clustering algorithm ${comTITLES[$k]}:"

	    N=$(awk 'END {print NR}' $COMDIR/${SUBDIR[$g]}/${comFILES[$k]})
	    N=$(($N-1))

	    echo "Network Size N: $N"
	
	    #---loop over annotation set 
	    for a in `seq $START $END`
	    do    
		echo "Annotation set: ${annoTITLES[$a]}"

		#---Create output directory for cluster enrichment results 
		OUTDIR=$GRAPHENRICHOUT/${comTITLES[$k]}
		if [ ! -d $OUTDIR ]; then
		    mkdir $OUTDIR
		fi


		if [[ ${comON[$k]} -eq "1" && ${annoON[$a]} -ge "1" ]]; then
		
		
		    #---clear tmp dir
		    rm $tmpDIR/*.csv
		    
		    #---run overlap between two annotation files in community relative to network
		    time ./$EXE -Comfile $COMDIR/${SUBDIR[$g]}/${comFILES[$k]} -Annofile $ANNODIR/${annoFILES[$KEYind]} -Annofile $ANNODIR/${annoFILES[$a]} -ext "${annoTITLES[$KEYind]}_${annoTITLES[$a]}" -o $tmpDIR -opt 2 -minOV $MINOV1 -minOV $MINOV2 -minOV $MINOV3 -Chi2 -twosided
			
		    #---Copy cluster enrichment results to output directory
		    cp $WORKINGDIR/$tmpDIR/* $OUTDIR
		
		    #---clear OUT/ dir
		    rm $tmpDIR/*.csv
	    		  
		    fi
		    
	    done
	done
	

	#---loop over Regions
	for k in `seq $rSTART $rEND`
	do
	    echo "Running Regions over clustering algorithm ${regTITLES[$k]}:"

	    N=$(awk 'END {print NR}' $REGDIR/${SUBDIR[$g]}/REGIONS/${regFILES[$k]})
	    N=$(($N-1))
	    
	    echo "Network Size N: $N"
	
	    #---loop over annotation set 
	    for a in `seq $START $END`
	    do    
		echo "Annotation set: ${annoTITLES[$a]}"

		#---Create output directory for cluster enrichment results 
		OUTDIR=$REGIONENRICHOUT/${regTITLES[$k]}
		if [ ! -d $OUTDIR ]; then
		    mkdir $OUTDIR
		fi
		

		if [[ ${regON[$k]} -eq "1" && ${annoON[$a]} -ge "1" ]]; then
		
		    #---clear tmp dir
		    rm $tmpDIR/*.csv

		    #---overlap between two annotation files in Bridging Region relative to network
		    time ./$EXE -Comfile $REGDIR/${SUBDIR[$g]}/REGIONS/${regFILES[$k]} -Annofile $ANNODIR/${annoFILES[$KEYind]} -Annofile $ANNODIR/${annoFILES[$a]} -ext "${annoTITLES[$KEYind]}_${annoTITLES[$a]}" -o $tmpDIR -opt 2 -minOV $MINOV1 -minOV $MINOV2 -minOV $MINOV3 -offset $OFFSET -Chi2 -twosided

		    #---regional enrichment for annotation set 1 
		    time ./$EXE -Comfile $REGDIR/${SUBDIR[$g]}/REGIONS/${regFILES[$k]} -Annofile $ANNODIR/${annoFILES[$KEYind]} -ext ${annoTITLES[$KEYind]} -o $tmpDIR -opt 1 -setFDR BY -noPerm -printAn -onesided

		    #---regional enrichment for annotation set 2 
		    time ./$EXE -Comfile $REGDIR/${SUBDIR[$g]}/REGIONS/${regFILES[$k]} -Annofile $ANNODIR/${annoFILES[$a]} -ext ${annoTITLES[$a]} -o $tmpDIR -opt 1 -setFDR BY -noPerm -printAn -onesided
		    
		    #---Copy cluster enrichment results to output directory
		    cp $WORKINGDIR/$tmpDIR/* $OUTDIR

		    #---clear OUT/ dir
		    rm $tmpDIR/*.csv
		    

		    #if [[ "$KEYind" -ne "$a" ]]; then

		    #---run overlap between annotation files at network level 
		    time ./run -Comfile $REGDIR/${SUBDIR[$g]}/REGIONS/${regFILES[$k]} -Annofile $ANNODIR/${annoFILES[$KEYind]} -Annofile $ANNODIR/${annoFILES[$a]} -ext "${annoTITLES[$KEYind]}_${annoTITLES[$a]}" -o $tmpDIR -opt 3 -minOV $MINOV1 -minOV $MINOV2 -minOV $MINOV3 -setN "${annoTITLES[$KEYind]}" -setN "${annoTITLES[$a]}" -offset $OFFSET 

		    
	    	    #---Copy cluster enrichment results to output directory
		    cp $WORKINGDIR/$tmpDIR/* $OUTDIR

		    #---clear OUT/ dir
		    rm $tmpDIR/*.csv

		    #fi
		    
		fi
		    
	    done
	done
	
    fi
     
done
    
echo "$0 done!"

exit 0
    
