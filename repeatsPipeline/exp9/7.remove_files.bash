#!/bin/bash

echo -e
echo "### 7.remove_files.bash ###"
echo -e
echo -e

rmFile=$1
uID=$2
scriptDir=$3

rm $rmFile

echo $uID >> $scriptDir/done_samples.txt