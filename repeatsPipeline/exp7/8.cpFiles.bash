#!/bin/bash

expName="exp7"

homeDir1="/share/ClusterShare/thingamajigs/jamtor/"
projectDir1="$homeDir1/projects/hgsoc_repeats/RNA-seq/"
reportsDir1="$projectDir1/reports/$expName"
htseqDir1="$projectDir1/results/htseq/$expName"

homeDir2="/share/ScratchGeneral/jamtor/"
projectDir2="$homeDir2/projects/hgsoc_repeats/RNA-seq/"
reportsDir2="$projectDir2/reports/$expName"
htseqDir2="$projectDir2/results/htseq/$expName"

mkdir -p $reportsDir2
cp -r $reportsDir1/* $reportsDir2

mkdir -p $htseqDir2
cp -r $htseqDir1/* $htseqDir2