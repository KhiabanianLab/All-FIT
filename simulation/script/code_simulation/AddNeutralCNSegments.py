#!/usr/bin/python
# Assumes input segments are sorted in each chr.

import subprocess as sp
import sys

script_directory = "/home/jl1444/"
loc = sp.Popen(["pwd"],stdout=sp.PIPE).communicate()[0].strip()+"/"
input_file = sys.argv[1]

chr_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]

fk = open(loc+input_file,"r")
data = fk.readlines()
header = data[0]
content = data[1:]
fk.close()

cnv_dic = {}	# chr as key, [[start,end,num_probes,seg_mean]]
for row in content:
	row = row.strip().split("\t")
	if row[0] not in cnv_dic.keys():
		cnv_dic[row[0]] = [[row[1],row[2],row[3],row[4]]]
	else:
		cnv_dic[row[0]].append([row[1],row[2],row[3],row[4]])

tmp_name = input_file.split(".")
out = open(loc+".".join(tmp_name[0:-2])+".addNeutral."+".".join(tmp_name[-2:]),"w")
out.write(header)
for i in chr_list:
	fh = open(script_directory+"/All-FIT/simulation/SureSelect_v5/"+i+"/seg_probes","r")
	ref = fh.readlines()
	fh.close()

	if i not in cnv_dic.keys():
		start = ref[0].strip().split("\t")[0]
		end = ref[-1].strip().split("\t")[1]
		num_probe = len(ref)
		out.write(i+"\t"+start+"\t"+end+"\t"+str(num_probe)+"\t0\n")
	else:
		ref_pointer = 0
		for cnv_ind in range(len(cnv_dic[i])):
			num_probe = 0
			before = int(cnv_dic[i][cnv_ind][0])
			if cnv_ind == 0:
				start = ref[0].strip().split("\t")[0]
				for line in range(len(ref)):
					row = ref[line].strip().split("\t")
					if int(row[1]) < before:
						num_probe += 1
					else:
						ref_pointer = line
						break
				if num_probe > 0:
					out.write(i+"\t"+start+"\t"+str(before-1)+"\t"+str(num_probe)+"\t0\n")
			else:
				after = int(cnv_dic[i][cnv_ind-1][1])				
				for line in range(ref_pointer,len(ref),1):
					row = ref[line].strip().split("\t")
					if int(row[0]) > after and int(row[1]) < before:
						num_probe += 1
					elif int(row[1]) > before:
						end = ref[line-1].strip().split("\t")[1]
						ref_pointer = line
						break
				if num_probe > 0:
					out.write(i+"\t"+str(after+1)+"\t"+str(before-1)+"\t"+str(num_probe)+"\t0\n")
			out.write(i+"\t"+"\t".join(cnv_dic[i][cnv_ind])+"\n")
		ending = int(cnv_dic[i][-1][1])
		last_ref = ref[-1].strip().split("\t")[1]
		num_probe = 0
		for line in range(ref_pointer,len(ref),1):
			row = ref[line].strip().split("\t")
			if int(row[0]) > ending:
				num_probe += 1
		out.write(i+"\t"+str(ending+1)+"\t"+last_ref+"\t"+str(num_probe)+"\t0\n")
out.close()		
