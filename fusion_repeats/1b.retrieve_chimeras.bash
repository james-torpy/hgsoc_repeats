#!/bin/bash

### 1b.retrieve_chimeras.bash ###

# Takes output from 1a.transfer_chimeras_from_nci_mdss.bash,
# untars it and retrieves Chimeric.out.junction

# fetch input values:
clusterDir=$1
chimeraDir=$2
uID=$3

echo $clusterDir
echo $chimeraDir
echo $uID

# untar file:
tar -zxvf "$clusterDir/$uID.tar.gz" -C $clusterDir

# make chimera directory:
mkdir -p "$chimeraDir/$uID/"

# copy Chimeric.out.junction to chimera directory:
cp $clusterDir/$clusterDir/$uID/Chimeric.out.junction $chimeraDir/$uID/
if [ -f "$chimeraDir/$uID/Chimeric.out.junction" ]; then
	rm -r "$clusterDir/$clusterDir/$uID/"
fi;

rm "$clusterDir/$uID.tar.gz"