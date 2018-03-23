#!/bin/bash

### transfer_gc_from_nci_mdss.bash ###

scgDir="/share/ScratchGeneral/jamtor/"

inDir="$scgDir/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/logs/"
shortDir="/short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"
clusterDir="$scgDir/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"

# get each file uID:
files=("$inDir/rfPT6" "$inDir/rfPT7" "$inDir/rfPT8" "$inDir/rfPT9")
#files=( $(ls $inDir/* | grep ^.*[A-Z][A-Z][0-9].* | sed "s/\\://") )
#files=( $(ls $inDir/* | grep FT | sed "s/\\://") )

echo -e
echo "The files are:"
echo ${files[@]}

for file in ${files[@]}; do 
	uID=`basename $file`
	#uID="prPT11"
	#uID=`basename $file | sed "s/\\.tar.gz//"`
	echo -e
	echo "The uID is: $uID"

	# transfer out of mass storage:
	#echo -e
	#echo "Transferring jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"$uID".tar.gz out of NCI mass storage..."
	#mdssLine="mdss get jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"$uID".tar.gz /short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"
	#ssh jt3341@raijin.nci.org.au "$mdssLine"

	# transfer to cluster:
	echo -e
	echo "Transferring jt3341@raijin.nci.org.au:$shortDir/$uID.tar.gz to $clusterDir..."
	rsync -avPS jt3341@raijin.nci.org.au:$shortDir/$uID.tar.gz $clusterDir

	if [ -f "$clusterDir/$uID.tar.gz" ]; then
		echo -e
		echo "Removing $uID.tar.gz from NCI short"
		ssh jt3341@raijin.nci.org.au "rm $shortDir/$uID.tar.gz"
		echo -e
		echo "Untarring $clusterDir/$uID.tar.gz"	
		tar -zvxf $clusterDir/$uID.tar.gz -C $clusterDir
		echo "Removing $clusterDir/$uID.tar.gz"
		rm $clusterDir/$uID.tar.gz
	fi;

done;