#!/bin/bash

mdb_path="/afs/crc.nd.edu/user/r/rdefever/repos/MiniDrugBank/split/molecules"

for filen in {1..371}
do
    cname="mdb_${filen}"
    echo "STARTING COMPOUND: ${cname}"
    antechamber -i ${mdb_path}/${cname}.mol2 -fi mol2 -o ${cname}.mol2 -fo mol2 -at gaff -s 2
done &> antechamber_log.txt


