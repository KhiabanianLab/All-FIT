#!/bin/bash

input=$1
python ../script/code_simulation/part1_generate_ABSOLUTE_input.py $input
bash ../script/code_simulation/part2_sort_ABSOLUTE_input.sh $input
