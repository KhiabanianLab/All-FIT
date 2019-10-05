#!/usr/bin/python

import subprocess as sp
import numpy as np, sys
import random

loc = sp.Popen(['pwd'],stdout=sp.PIPE).communicate()[0].decode().strip()+"/"

input_file = sys.argv[1]

fh = open(loc+"PuritySim_Set3/our_method/"+input_file,"r")
content = fh.readlines()
header = content[0]
data = content[1:]
fh.close()

for i in range(100):	#percentage of mutation being kept
	for j in range(100):	#number of replicates
		out = open(loc+"PuritySim_Set10/"+input_file+"/"+input_file+".%i.%i"%(i,j),"w")
		out.write(header)
		index_list = random.sample(range(0,len(data)-1),round(i*0.01*len(data)))
		for ind in index_list:
			out.write(data[ind])
		out.write(data[-1])
		out.close()
