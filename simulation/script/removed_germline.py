#!/usr/bin/python

import subprocess as sp
import sys

loc=sp.Popen(["pwd"],stdout=sp.PIPE).communicate()[0].strip()

input_file = sys.argv[1]
fh = open(loc+"/"+input_file,"r")
data = fh.readlines()[1:]
fh.close()

out = open(loc+"-removed-germline/"+input_file,"w")
out.write("ID\tAllele_Freq\tDepth\tPloidy\tModel\n")
for row in data:
	ele = row.strip().split("\t")
	if ele[4].split(",")[0] != "germline":
		out.write(row)
out.close()
