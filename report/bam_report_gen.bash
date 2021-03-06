#!/bin/bash

### bam_report_gen.bash ###

# For each sample subset counts alignments in mapped bam,primary alignments in bam,
# primary alignments reported by STAR and multimappers reported by STAR and
# generates a report to be inputted into method_report_gen.R
# Also sorts each bam for R script.

# define variables:
projectName="hgsoc_repeats"
methodName="method2"
bamExt=".bam"
numcores=6

# define directories:
homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectName"
resultsDir="$projectDir/RNA-seq/results"
bamPath="$resultsDir/star/$methodName"
reportDir="$resultsDir/method_reports/$methodName"
scriptDir="$projectDir/RNA-seq/scripts"
varDir="$scriptDir/R/vars"
logDir="$scriptDir/logs"

mkdir -p $reportDir
mkdir -p $varDir
mkdir -p $logDir

RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"

echo -e
echo "This is the bamPath:"
echo $bamPath
echo -e
echo "This is the reportDir:"
echo $reportDir

# fetch directory names for files:
i=0
directory_names=`ls $bamPath | grep -v subset | grep -v FT[2-3] | grep -v PR[3-4] | grep -v grant_endo[2-3] | grep -v grant_endosis[2-3] | grep -v gtx_ft[2-3] | grep -v kur_v[2-3]`
for directory in ${directory_names[@]}; do
	echo -e
	echo The directory used is: $directory
	directoriesTotal[i]=$directory
	let i++
done;

echo -e
echo The total number of directories is: ${#directoriesTotal[@]}

# fetch inFiles:
j=0
while [ $j -lt ${#directoriesTotal[@]} ]; do
	uniqueDir="$bamPath/${directoriesTotal[$j]}"
	uID=`basename $uniqueDir`
	inFile=$uniqueDir/Aligned.sortedByCoord.out.bam

	echo -e
	echo "The file used is:"
	echo $inFile

	# fetch number of input reads reported by STAR (6th line, field 2 of Log.final.out),
	# and add this to the report:
	input=$(sed '6q;d' $uniqueDir/Log.final.out | cut -f2)
	echo -e >> "$reportDir/"$uID\_report.txt
	echo -e "No. input reads reported by STAR:\t$input" > "$reportDir/"$uID\_report.txt
	
	# fetch number of uniquely mapped reads reported by STAR (9th line, field 2 of Log.final.out),
	# and add this to the report:
	unique=$(sed '9q;d' $uniqueDir/Log.final.out | cut -f2)
	echo -e "No. uniquely mapped reads reported by STAR:\t$unique" >> "$reportDir/"$uID\_report.txt
	
	# fetch number of multimapping reads reported by STAR (24th line, field 2 of Log.final.out),
	# and add this to the report:
	multi=$(sed '24q;d' $uniqueDir/Log.final.out | cut -f2)
	echo -e "No. multimapping reads reported by STAR:\t$multi" >> "$reportDir/"$uID\_report.txt

	# fetch mismatch rate per base reported by STAR (18th line, field 2 of Log.final.out),
	# and add this to the report:
	multi=$(sed '18q;d' $uniqueDir/Log.final.out | cut -f2)
	echo -e "Mismatch rate per base, % reported by STAR:\t$multi" >> "$reportDir/"$uID\_report.txt

	# convert inFile to .sam and count lines, adding this to the report generated by
	# fastq_report_gen.bash:
	count=$(samtools view $inFile | wc -l)
	echo -e >> "$reportDir/"$uID\_report.txt
	echo -e "Number of alignments in bam:\t$count" >> "$reportDir/"$uID\_report.txt
	
	# convert inFile to .sam and count primary alignment lines, adding this to the report:
	primary=$(samtools view $inFile | grep -v 339 | grep -v 355 | grep -v 403 | grep -v 419 | wc -l)
	echo -e "Number of primary alignments in bam:\t$primary" >> "$reportDir/"$uID\_report.txt

	# define report file variable to be passed to method_reports.R:
	rFile="$reportDir/$uID"_report.txt
	# define script to be called:
	outScript="$scriptDir/R/mapping_report.R"
	# define submitLine:
	#submitLine="--vanilla < $outScript "$varDir/$uID"_vars.R"
	submitLine="--vanilla --args $methodName $rFile $uID < $outScript"

	echo -e
	echo "The submitLine is:"
	echo $submitLine

	# submit method_reports.R to the cluster for each sample, with associated variables:
	qsub -b y -j y -N rep$uID -wd $logDir -pe smp $numcores -V $RDir $submitLine

	j=$(($j+1))
done;





