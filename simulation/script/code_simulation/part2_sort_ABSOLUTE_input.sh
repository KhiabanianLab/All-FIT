#!/bin/bash

loc=$(pwd)
cnv="tmp."$1".cnv.txt"
snv="tmp."$1".maf.txt"
out_cnv=$(echo ${cnv}|awk -F "tmp." '{print $2}')
out_snv=$(echo ${snv}|awk -F "tmp." '{print $2}')
cat ${loc}"/ABSOLUTE/"${snv}|(read -r; printf "%s\n" "$REPLY"; sort -t $'\t' -k7,7V -k4,4n) > ${loc}"/ABSOLUTE/"${out_snv}
cat ${loc}"/ABSOLUTE/"${cnv}|(read -r; printf "%s\n" "$REPLY"; sort -t $'\t' -k1,1V -k2,2n) > ${loc}"/ABSOLUTE/"${out_cnv}
rm ${loc}"/ABSOLUTE/"${snv}
rm ${loc}"/ABSOLUTE/"${cnv}
