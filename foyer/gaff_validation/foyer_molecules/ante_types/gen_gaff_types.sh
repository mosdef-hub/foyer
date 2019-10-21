#!/bin/bash

ls ../mol2/*mol2| while read filen
do
    cname=`echo ${filen} | sed 's/\.mol2//' | sed 's/\.\.\/mol2\///'`
    echo "STARTING COMPOUND: ${cname}"
    antechamber -i ../mol2/${cname}.mol2 -fi mol2 -o ${cname}.mol2 -fo mol2 -at gaff -s 2
done &> antechamber_log.txt


