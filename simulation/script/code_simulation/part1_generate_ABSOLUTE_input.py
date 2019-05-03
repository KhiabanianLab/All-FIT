#!/usr/bin/python

import random,sys,math
import subprocess as sp

loc = sp.Popen(["pwd"],stdout=sp.PIPE).communicate()[0].strip()+"/"

#Reading in seg_probes files to randomly chose coding position for writing into ABSOLUTE input
chr_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"]

chr_pos_dic = {}     #chr as key, [[start,end]] as value
for i in chr_list:
	fh = open("/home/jl1444/project/automerging_project/data_for_paper/SureSelect_v5/"+i+"/seg_probes","r")

	data = fh.readlines()
	fh.close()

	chr_pos_dic[i] = []
	for row in data:
		row = row.strip().split("\t")
		chr_pos_dic[i].append([int(row[0]),int(row[1])])

input_file = sys.argv[1]+".xls"
fk = open(loc+"our_method/"+input_file,"r")
ref = fk.readlines()[1:]
fk.close()
freq = []	#in decimal form
depth = []
ploidy = []
for row in ref:
	row = row.strip().split("\t")
	if row[3].split(",")[0] == "somatic":
		freq.append(float(row[0])/100)
		depth.append(int(row[1]))
		ploidy.append(int(row[2]))

#generating ABSOLUTE cnv input
tmp_name = input_file.split(".xls")[0]
file_cnv = loc+"ABSOLUTE/tmp."+tmp_name+".cnv.txt"
out_cnv = open(file_cnv,"w")
out_cnv.write("Chromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n")
corresponding_snv_list = []	#stored in format of snv file
for ind in range(len(ploidy)):
	if ploidy[ind] != 2:
		alt_count = int(round(freq[ind]*depth[ind]))
		ref_count = depth[ind]-alt_count
		rand_chr = random.choice(chr_list)
		rand_exon_index = random.randrange(len(chr_pos_dic[rand_chr]))
		start = chr_pos_dic[rand_chr][rand_exon_index][0]
		end = chr_pos_dic[rand_chr][rand_exon_index][1]
		rand_position = random.randint(start,end)
		out_cnv.write(rand_chr+"\t"+str(rand_position)+"\t"+str(rand_position+99)+"\t1\t"+str(math.log(ploidy[ind]*1.0/2,2))+"\n")
		corresponding_snv_list.append(str(ref_count)+"\t"+str(alt_count)+"\t\t"+str(rand_position)+"\t\t\t"+rand_chr)
		del chr_pos_dic[rand_chr][rand_exon_index]
		chr_pos_dic[rand_chr].append([start,rand_position-1])
		if rand_position+100 <= end:
			chr_pos_dic[rand_chr].append([rand_position+100,end])
out_cnv.close()

#generating ABSOLUTE snv input
file_snv = loc+"ABSOLUTE/tmp."+tmp_name+".maf.txt"
out_snv = open(file_snv,"w")
out_snv.write("t_ref_count\tt_alt_count\tdbSNP_Val_Status\tStart_position\tTumor_Sample_Barcode\tHugo_Symbol\tChromosome\n")
for ind in range(len(ploidy)):
	if ploidy[ind] == 2:
		alt_count = int(round(freq[ind]*depth[ind]))
		ref_count = depth[ind]-alt_count
		rand_chr = random.choice(chr_list)
		rand_exon_index = random.randrange(len(chr_pos_dic[rand_chr]))
		start = chr_pos_dic[rand_chr][rand_exon_index][0]
		end = chr_pos_dic[rand_chr][rand_exon_index][1]
		rand_position = random.randint(start,end)
		out_snv.write(str(ref_count)+"\t"+str(alt_count)+"\t\t"+str(rand_position)+"\t\t\t"+rand_chr+"\n")
		del chr_pos_dic[rand_chr][rand_exon_index]
		chr_pos_dic[rand_chr].append([start,rand_position-1])
		if rand_position+1 <= end:
			chr_pos_dic[rand_chr].append([rand_position+1,end])
if len(corresponding_snv_list) > 0:
	out_snv.write("\n".join(corresponding_snv_list)+"\n")
out_snv.close()
