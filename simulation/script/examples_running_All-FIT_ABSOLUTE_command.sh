#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -t 0:30:00
#SBATCH -p perceval
#SBATCH --export=ALL
#SBATCH -D /home/jl1444/project/automerging_project/data_for_paper/simulation/All-FIT/simulation/log

script_directory=/home/jl1444
cd ${script_directory}/All-FIT/simulation
echo [Directory]`pwd`
starting=$(date)
echo "[Start] "$starting

input_file=$1
input_folder=$2
#python ../All-FIT.py -i ${script_directory}/All-FIT/simulation/PuritySim_Set${input_folder}/our_method/${input_file}".xls" -d ${script_directory}/All-FIT/simulation/output_PuritySim_Set${input_folder}/our_method -o ${input_file}

python ../All-FIT.py -i ${script_directory}/All-FIT/simulation/PuritySim_Set${input_folder}/our_method-removed-germline/${input_file}".xls" -d ${script_directory}/All-FIT/simulation/output_PuritySim_Set${input_folder}/our_method-removed-germline -o ${input_file} -t somatic

#Rscript script/ABSOLUTE_command.r ${script_directory}/All-FIT/simulation/PuritySim_Set${input_folder}/ABSOLUTE/ ${script_directory}/All-FIT/simulation/output_PuritySim_Set${input_folder}/ABSOLUTE/ ${input_file}

end=$(date)
echo "[End] "$end

