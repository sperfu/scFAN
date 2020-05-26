#!/bin/bash

## change directory
cd Data
mkdir new_TF_for_GM12878
cd new_TF_for_GM12878
## download TF TBX21 BATF BMI1 CBFB SMAD1 in cell line GM12878

## TBX21
wget https://www.encodeproject.org/files/ENCFF902EUL/@@download/ENCFF902EUL.bed.gz
## CBFB
wget https://www.encodeproject.org/files/ENCFF838ABE/@@download/ENCFF838ABE.bed.gz
## SMAD1
wget https://www.encodeproject.org/files/ENCFF916SIZ/@@download/ENCFF916SIZ.bed.gz
## BMI1
wget https://www.encodeproject.org/files/ENCFF159QMJ/@@download/ENCFF159QMJ.bed.gz
## BATF
wget https://www.encodeproject.org/files/ENCFF887VZZ/@@download/ENCFF887VZZ.bed.gz

## rename files
mv ENCFF902EUL.bed.gz GM12878_TBX21.bed.gz
mv ENCFF838ABE.bed.gz GM12878_CBFB.bed.gz
mv ENCFF916SIZ.bed.gz GM12878_SMAD1.bed.gz
mv ENCFF159QMJ.bed.gz GM12878_BMI1.bed.gz
mv ENCFF887VZZ.bed.gz GM12878_BATF.bed.gz

## generate chip.txt
ls *bed.gz | awk -F '\_|\.' '{OFS="\t"}{print $0,$2,$0}' > chip.txt

## copy other files here
cp ../GM12878/bigwig.txt ./
ln -s ../GM12878/SRR891269.forward.1x2.bw SRR891269.forward.1x2.bw
ln -s ../wgEncodeDukeMapabilityUniqueness35bp.bigWig wgEncodeDukeMapabilityUniqueness35bp.bigWig

## ready for training
cd ../../

python scFAN_train.py -i Data/new_TF_for_GM12878 -e 5 -oc new_TF_GM12878


