#!/bin/bash

### call_compBarplot.bash ###

# Call loadBams.R to load bams and save as RDS files

# define variables:
projectName="hgsoc_repeats"
methodName="method1"
sampleType=""
numcores=3

# define directories:
homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectName/"
resultsDir="$projectDir/RNA-seq/results/"
RobjectDir="$projectDir/RNA-seq/Robjects/"
scriptDir="$projectDir/RNA-seq/scripts/R/"
logDir="$scriptDir/logs/"

mkdir -p $logDir

RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"
outScript="$scriptDir/countGencodeRepeatsRibo.R"

echo -e
echo "This is the RobjectDir:"
echo $RobjectDir
echo -e
echo "This is the logDir:"
echo $logDir

# fetch inFiles:
inFiles=( `ls $RobjectDir | grep $sampleType\_bamGR.rds | grep -v bowtell_FT2_bamGR | grep -v subset | grep -v CG | grep -v lib` )
echo -e
echo "The inFiles are:"
echo ${inFiles[@]}

if [ sampleType="" ]; then
	sampleType="full"
fi

for file in ${inFiles[@]}; do
	uID=`basename $file | sed 's/_bamGR.rds//'`
	# create command to send bam file to loadBams.R for loading:
	submitLine="--vanilla --args $file $uID $sampleType < $outScript"

	echo -e
	echo "The submitLine is:"
	echo $submitLine

	# submit method_reports.R to the cluster for each sample, with associated variables:
	qsub -b y -j y -N compCount_$sampleType -wd $logDir -pe smp $numcores -V $RDir $submitLine

done;


