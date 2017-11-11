#!/bin/bash

### 4.bam_report_gen.bash ###

# For each sample subset counts alignments in mapped bam,primary alignments in bam,
# primary alignments reported by STAR and multimappers reported by STAR and
# generates a report to be inputted into method_report_gen.R
# Also sorts each bam for R script.

echo -e
echo "### 4.bam_report_gen.bash ###"
echo -e

# make directory hierachy:
projectname="hgsoc_repeats"
sampleType="fullsamples/bowtell_primary"
expName="exp5"

homeDir="/share/ClusterShare/thingamajigs/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
rawDir="$projectDir/raw_files/"
resultsDir="$projectDir/results/"

# scripts/logs directories:
scriptsPath="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/"

# set in/outPaths:
bamPath="$resultsDir/star/GC/$expName/$uID"
reportDir="$projectDir/reports/$expName/"

mkdir -p $reportDir

# fetch uID:
uID=$1
fastqPath=$2

# remove fastq files and non-sorted bam as they are no longer needed:
#echo -e
#echo "Removing fastq files..."
#rm $fastqPath/$uID.1.fastq.gz
#rm $fastqPath/$uID.2.fastq.gz
#rm $bamPath/Aligned.out.bam

echo -e
echo "This is the uID:"
echo $uID
echo -e
echo "This is the reportDir:"
echo $reportDir
echo -e

mkdir -p $reportDir

echo -e
echo "This is the bamPath:"
echo $bamPath
echo -e
echo "This is the reportDir:"
echo $reportDir

uniqueDir="$bamPath/$uID"
inFile="$uniqueDir//Aligned.sortedByName.out.bam"

echo -e
echo "This is the inFile:"
echo $inFile

# fetch number of input reads reported by STAR (6th line, field 2 of Log.final.out),
# and add this to the report:
input=$(sed '6q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e No. input reads reported by STAR:\\t$input
echo -e No. input reads reported by STAR:\\t$input >> $reportDir/$uID.report.txt

# fetch number of uniquely mapped reads reported by STAR (9th line, field 2 of Log.final.out),
# and add this to the report:
unique=$(sed '9q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e No. uniquely mapped reads reported by STAR:\\t$unique
echo -e No. uniquely mapped reads reported by STAR:\\t$unique >> $reportDir/$uID.report.txt

# fetch percentage of uniquely mapped reads reported by STAR (10th line, field 2 of Log.final.out),
# and add this to the report:
uniqueP=$(sed '10q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e Percentage uniquely mapped reads reported by STAR:\\t$uniqueP
echo -e Percentage uniquely mapped reads reported by STAR:\\t$uniqueP >> $reportDir/$uID.report.txt

# fetch number of multimapping reads reported by STAR (24th line, field 2 of Log.final.out),
# and add this to the report:
multi=$(sed '24q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e No. multimapping reads reported by STAR:\\t$multi
echo -e No. multimapping reads reported by STAR:\\t$multi >> $reportDir/$uID.report.txt

# fetch percentage of multimapping reads reported by STAR (25th line, field 2 of Log.final.out),
# and add this to the report:
multiP=$(sed '25q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e Percentage multimapping reads reported by STAR:\\t$multiP
echo -e Percentage multimapping reads reported by STAR:\\t$multiP >> $reportDir/$uID.report.txt

# fetch percentage of unmapped reads reported by STAR (30th line, field 2 of Log.final.out),
# and add this to the report:
nonMap=$(sed '30q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e Percentage unmapped reads reported by STAR:\\t$nonMap
echo -e Percentage unmapped reads reported by STAR:\\t$nonMap >> $reportDir/$uID.report.txt

# fetch mismatch rate per base reported by STAR (18th line, field 2 of Log.final.out),
# and add this to the report:
misM=$(sed '18q;d' $uniqueDir/Log.final.out | cut -f2)
echo -e Mismatch rate per base, % reported by STAR:\\t$misM
echo -e Mismatch rate per base, % reported by STAR:\\t$misM >> $reportDir/$uID.report.txt




