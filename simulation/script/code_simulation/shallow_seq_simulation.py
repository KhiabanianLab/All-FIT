#!/usr/bin/python

import subprocess as sp
from scipy.stats import binom
import numpy as np, sys

loc = sp.Popen(['pwd'],stdout=sp.PIPE).communicate()[0].decode().strip()+"/"
alpha = float(sys.argv[2]) #0.1 - taking 10% of total depth
trial = sys.argv[3]

input_file = sys.argv[1]
fh = open(loc+"PuritySim_Set3/our_method/"+input_file,"r")
content = fh.readlines()
header = content[0]
data = content[1:]
fh.close()

out = open(loc+"PuritySim_Set9/trial"+trial+"/alpha"+str(alpha)+"/"+input_file,"w")
out.write(header)
for row in data:
	row = row.strip().split("\t")
	n_depth = int(round(alpha * int(row[2])))
	AF = float(row[1])/100.0
	n_vardepth = float(int(np.random.binomial(n_depth,AF)))
	n_AF = round(n_vardepth/n_depth,2)
#	if row[0] == "Var2":	#specifically for SIM_DATA_3.7310.xls which has AF=0
#		continue
	if n_AF == 0:
		while True:
			n_vardepth = float(int(np.random.binomial(n_depth,AF)))
			n_AF = round(n_vardepth/n_depth,2)
			if n_AF > 0:
				break
	out.write(row[0]+"\t"+str(round(n_AF*100.0,2))+"\t"+str(n_depth)+"\t"+"\t".join(row[3:])+"\n")
out.close()
