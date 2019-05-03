#!/bin/bash

folder=$1
sample=$2
if [ -d output_PuritySim_Set${folder} ];then 
	rm -r output_PuritySim_Set${folder}
	mkdir output_PuritySim_Set${folder}
	mkdir output_PuritySim_Set${folder}/our_method
	mkdir output_PuritySim_Set${folder}/ABSOLUTE
fi
python script/code_simulation/purity_simulation.py $folder $sample

cd "PuritySim_Set"${folder}"/our_method"
files=$(ls|awk -F ".xls" '{print $1}')
cd ..
for i in $files; do bash ../script/code_simulation/generate_ABSOLUTE_input.sh $i;done
cd ABSOLUTE
cnv=$(ls *cnv*)
for i in $cnv; do ../../script/code_simulation/AddNeutralCNSegments.py $i;done
