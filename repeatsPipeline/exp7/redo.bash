# make directory hierachy:
projectname="hgsoc_repeats"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts"
logDir="$scriptsPath/repeatsPipeline/exp7/logs/"
oldDir="$logDir/old/"

mkdir -p $oldDir

source $scriptsPath/repeatsPipeline/exp7/test.bash



# fetch logs with errors, and identify the samples to be redone:
for f in $logDir/*; do
	err=`grep -E 'error|ERROR|Error|halted|unexpected' $f | grep -v ERROR: | grep -v error:`
	echo ${err[@]} >> /home/jamtor/test.txt
	



	if [ ${#err} -gt 1 ]; then
		id=`basename $f | sed s'/^.*_//' | sed 's/.o[0-9].*//'`
		echo -e
		echo $id >> /home/jamtor/test.txt

		if [ -e "$rawDir/fullsamples/bowtell_primary/record.txt" ]; then
			g=( $(grep $id "$rawDir/fullsamples/bowtell_primary/record.txt") )
			len=`echo ${g[@]} | wc -w`
			echo -e
			echo $id >> "$rawDir/fullsamples/bowtell_primary/files.txt"
			echo $id >> "$rawDir/fullsamples/bowtell_primary/record.txt"
			echo $id >> /home/jamtor/test.txt



			export largecores=$((4*$len+8))
			export countcores=$((4*$len+14))

			echo "$id is being re-done with $largecores largecores and $countcores countcores, number of redos = $len" >> /home/jamtor/test.txt
			echo "$id is being re-done with $largecores largecores and $countcores countcores, number of redos = $len"



		else echo -e
			 echo "$id is being re-done with 8 largecores and 14 countcores"
			 echo $rawDir/fullsamples/bowtell_primary/$id.bam >> "$rawDir/fullsamples/bowtell_primary/files.txt"
			 echo $id >> "$rawDir/fullsamples/bowtell_primary/record.txt"

			echo $rawDir/fullsamples/bowtell_primary/$id.bam >> /home/jamtor/test.txt

			echo "$id is being re-done with 8 largecores and 14 countcores" >> /home/jamtor/test.txt

			 export largecores=8
			 export countcores=14

			 echo "$id is being re-done with $largecores largecores and $countcores countcores" >> /home/jamtor/test.txt

		fi;

	# move logs to oldDir:
	mv $logDir/*$id* $oldDir
	fi;
done;

# call repeats pipeline for samples with errors:
source $scriptsPath/repeatsPipeline/exp7/1.bam_to_fastq.bash


