# make directory hierachy:
projectname="hgsoc_repeats"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts"
logDir="$scriptsPath/repeatsPipeline/exp8/logs/"
oldDir="$logDir/old/"

mkdir -p $oldDir


# fetch logs with errors, and identify the samples to be redone:
for f in $logDir/*; do
	err=`grep -E 'error|ERROR|Error|halted|unexpected' $f | grep -v ERROR: | grep -v error:`
	if [ ${#err} -gt 1 ]; then
		id=`basename $f | sed s'/^.*_//' | sed 's/.o[0-9].*//'`
		if [ -e "$rawDir/fullsamples/bowtell_primary/record.txt" ]; then
			g=( $(grep $id "$rawDir/fullsamples/bowtell_primary/record.txt") )
			len=`echo ${g[@]} | wc -w`
			echo -e
			echo $id >> "$rawDir/fullsamples/bowtell_primary/files.txt"
			echo $id >> "$rawDir/fullsamples/bowtell_primary/record.txt"

			export fqcores=$((2*$len+4))
			export ribocores=$((2*$len+12))
			export gccores=$((2*$len+8))
			export sortcores=$((2*$len+9))
			export countcores=$((2*$len+14))
			export mem=$((4*$len+24))G

			echo "$id is being re-done with $fqcores fqcores,  $ribocores ribocores, \
			$gccores gccores, $sortcores sortcores, $countcores countcores and $mem memory, \
			number of redos = $len"

		else echo -e
			 echo "$id is being re-done with 8 largecores and 14 countcores"
			 echo $rawDir/fullsamples/bowtell_primary/$id.bam >> "$rawDir/fullsamples/bowtell_primary/files.txt"
			 echo $id >> "$rawDir/fullsamples/bowtell_primary/record.txt"

			export fqcores=4
			export ribocores=12
			export gccores=8
			export sortcores=9
			export countcores=14
			export mem=28G

			 echo "$id is being re-done with $fqcores fqcores,  $ribocores ribocores, \
			$gccores gccores, $sortcores sortcores, $countcores countcores and $mem memory"

		fi;

	# move logs to oldDir:
	mv $logDir/*$id* $oldDir
	fi;
done;

# call repeats pipeline for samples with errors:
source $scriptsPath/repeatsPipeline/exp8/1.bamtofastq.bash


