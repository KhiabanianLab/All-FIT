#!/usr/bin/python

'''
This file creates multiple fake mutations and gives them an allele frequency,
depth and purity. The allele frequency (VAF) and depth are pushed into one
file and all the correct purity's are placed into another file so that they
may be referenced later.

**This version gives every mutation an equal amount of being chosen and adds in
a few subclonal mutations
'''

from scipy.stats import binom
import random, math
import os.path
import sys, numpy as np

def mut_calc(mut_type,purity,out_sam):
	'''Calculates and write the VAF,depth,ploidy,model to the sample data files '''
	purity = round((purity*0.01), 2)
	count_mut = 1
	for m_type in range(len(mut_type)):
		num_mut = mut_type[m_type]
		while num_mut > 0:
			if m_type == 0:
				Y = 2
				cnmut = 1
				model = "somatic, CNmut = 1"
				#Randomly chooses between clonal and subcloanl mutaions only for this option
				#There is currently a bit less than 1/4 chance that this mutation turns subclonal.
				#This can be adjusted by changing the random int and while loop below
				sub = random.randint(1,8)
				if sub == 1 and purity > 0.2:
					subpurity = purity - (random.randint(5,10) * 0.01)
					model += ", subclonal"
					f = subpurity/2
				elif sub == 2 :
					subpurity = purity/2
					model += ", subclonal"
					f = subpurity/2
				else:
					#Calculates the expected allele frequency
					f = purity/2	#somatic, CNmut = 1
			elif m_type == 1:
				Y = 1
				cnmut = 1
				model = "somatic, LOH CNmut = 1"
				f = purity/(2-purity)    #somatic, LOH CNmut = 1
			elif m_type == 2:
				Y = 2
				cnmut = 2
				model = "somatic, CNmut = 2"
				f = purity       #somatic, CNmut = 2
			elif m_type == 3:
				Y = 2
				cnmut = 1
				model = "germline, CNmut = 1"
				f = 0.5  #germline, CNmut = 1
			elif m_type == 4:
				Y = 1
				cnmut = 1
				model = "germline, LOH CNmut = 1"
				f = 1/(2-purity) #germline, LOH CNmut = 1
			elif m_type == 5:
				Y = 2
				cnmut = 2
				model = "germline, CNmut = 2"
				f = (1 + purity)/2       #germline, CNmut = 2
			elif m_type == 6:
				Y = random.randint(3,8)
				cnmut = random.randint(1,Y)
				model = "somatic, CNmut=%i"%cnmut
				f = (cnmut * purity)/((2*(1-purity)) + (Y * purity))     #somatic
			elif m_type == 7:
				Y = random.randint(3,8)
				cnmut = random.randint(1,Y)
				model = "germline, CNmut=%i"%cnmut
				f = ((1-purity) + (cnmut * purity))/((2 * (1 - purity)) + (Y * purity))  #germline
			else:	#Assuming at least one somatic heterozygous mutation
				Y = 2
				cnmut = 1
				model = "somatic, CNmut = 1"
				f = purity/2	#somatic, CNmut = 1

			#The depth is just a random number from 300-1000
			depth = random.randint(300, 1000)
			VAF = round((float(np.random.binomial(depth,f))/depth),2)
			if model == "somatic, CNmut = 1, subclonal":
				if binom.cdf(round(depth*VAF),depth,purity/2) >= 0.01:
#					print sam_count, "subclonal is not right"
					while True:
						VAF = round((float(np.random.binomial(depth,f))/depth),2)
						if binom.cdf(round(depth*VAF),depth,purity/2) < 0.01:
							break
			out_sam.write("Var"+str(count_mut)+"\t"+str(VAF*100)+"\t"+str(depth)+"\t"+str(Y)+"\t"+model+"\n")
			num_mut -= 1
			count_mut += 1
	out_sam.close()

#Only change this number if you want an entire new set of data
set_count = int(sys.argv[1])
#Changes the amount of samples you want per answer file
samples = int(sys.argv[2])

curr_directory = os.getcwd()
folder = "PuritySim_Set"+str(set_count)
final_directory = os.path.join(curr_directory,folder)
if not os.path.exists(final_directory):
	os.makedirs(final_directory)
	os.makedirs(os.path.join(final_directory,"our_method"))
	os.makedirs(os.path.join(final_directory,"ABSOLUTE"))

out_ans = open(final_directory+"/SIM_DATA_ANS_"+str(set_count)+".xls","w")
out_ans.write("File_num\tPurity\n")

for sam_count in range(1,samples+1,1):
	#Selects a random purity from 10-90
	purity = random.randint(10,90)
	#Adds the correct purity to the answer file
	out_ans.write(str(set_count)+"."+str(sam_count)+"\t"+str(purity)+"\n")

	#Creates the sample data files
	out_sam = open(final_directory+"/our_method/SIM_DATA_"+str(set_count)+"."+str(sam_count)+".xls","w")
	out_sam.write("ID\tAllele_Freq\tDepth\tPloidy\tModel\n")	

	#Change the number of mutations per sample here:
	mut_count=random.randint(20,100)
#	mut_count=random.randint(5,100)	#for PuritySim_Set7 to test the effect on num_mutation_per_sample on All-FIT
	mut_type = [0 for x in range(9)]
	#Assume we have at least one somatic, CNmut = 1 (somatic heterozygous mutation)
	mut_count -= 1
	mut_type[-1] += 1
	for count in range(mut_count):
		mutation = random.randint(0,7)
		mut_type[mutation] += 1
	mut_calc(mut_type,purity,out_sam)

out_ans.close()
