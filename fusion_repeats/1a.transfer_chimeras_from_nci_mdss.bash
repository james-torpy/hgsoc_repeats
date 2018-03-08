#!/bin/bash

### 1a.transfer_chimeras_from_nci_mdss.bash ###

scgDir="/share/ScratchGeneral/jamtor/"

inDir="$scgDir/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/logs/"
shortDir="/short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"
clusterDir="$scgDir/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"
chimeraDir="$scgDir/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/chimeras/"
scriptDir="/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/fusion_repeats/"
logDir="$scgDir/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/qsub_logs/"

# get each file uID:
#files=( $(ls $inDir/* | grep ^.*[A-Z][A-Z][0-9].* | sed "s/\\://") )
files=( $(ls $inDir/* | grep FT | sed "s/\\://") )

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
#	mdssLine="mdss get jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"$uID".tar.gz /short/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC/exp9/"
#	ssh jt3341@raijin.nci.org.au "$mdssLine"

	# transfer to cluster:
#	rsync -avPS jt3341@raijin.nci.org.au:$shortDir/$uID.tar.gz $clusterDir

	if [ -f "$clusterDir/$uID.tar.gz" ]; then
#		echo -e
#		echo "Removing $uID.tar.gz from NCI short"
#		ssh jt3341@raijin.nci.org.au "rm $shortDir/$uID.tar.gz"

		# untar file, make chimera directory, and copy Chimeric.out.junction to chimera directory:
		qsub -N getch$uID -b y -wd $logDir -j y -R y -pe smp 1 -V \
		"$scriptDir/1b.retrieve_chimeras.bash" $clusterDir $chimeraDir $uID
			
	fi;

done;