#!/usr/bin/python

import argparse
import subprocess as sp
import numpy as np
import os
from scipy import stats
from scipy.stats import binom
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap
import seaborn as sns

__author__ = "Jui Wan Loh"
__date__ = "Date: 02-2019"

def binom_interval(n_success, total, conf_int):
	quantile = (1 - conf_int) / 2.
	lower = beta.ppf(quantile, n_success, total - n_success + 1)
	upper = beta.ppf(1 - quantile, n_success + 1, total - n_success)
	if n_success == total:
		upper = 1.0
	if n_success == 0:
		lower = 0.0
	return (lower, upper)

def gray_box_accuracy(vaf_list,depth_list,purity,ploidy_list,varType):
	vaf_CI = 0.01
	f_incre = 0.005
	CCF = [{} for j in range(len(vaf_list))]
	weight = [{} for j in range(len(vaf_list))]
	alpha = 1 - vaf_CI
	for j in range(len(vaf_list)):
		d = depth_list[j] 
		ploidy = ploidy_list[j]
		if ploidy == 2:
			if varType == "somatic":
				name_likelihood = ["Somatic, LOH CNmut=1"]
			else:
		        	name_likelihood = ["Somatic, LOH CNmut=1","Germline, LOH CNmut=1"]
		else:
                        name_likelihood = []
		som_list = []
		germ_list = []
		for i in range(ploidy):
			som_list.append("Somatic, CNmut=%i"%(i+1))
			germ_list.append("Germline, CNmut=%i"%(i+1))
		name_likelihood[1:1] = som_list
		if varType != "somatic":
			name_likelihood.extend(germ_list)

		for name in name_likelihood:
			weight[j][name] = 0.0
		if d == 0:
			continue
		if ploidy == 2:
			CCF[j]["Somatic, LOH CNmut=1"] = vaf_list[j]/(purity/(2-purity)) #somatic LOH
			if varType != "somatic":
				CCF[j]["Germline, LOH CNmut=1"] = ((vaf_list[j]*(2-purity))-1+purity)/purity	#germline LOH
		for i in range(ploidy):
			CCF[j]["Somatic, CNmut=%i"%(i+1)] = vaf_list[j]/(((i+1)*purity)/(2*(1-purity)+ploidy*purity))	#somatic LOH high CN
			if varType != "somatic":
				CCF[j]["Germline, CNmut=%i"%(i+1)] = ((vaf_list[j]*(2-2*purity+ploidy*purity))-1+purity)/((i+1)*purity)	#germline LOH high CN
				
		freqLB,freqUB = binom_interval(round(d*vaf_list[j]),d,alpha)
		f_range = np.arange(freqLB,freqUB+f_incre,f_incre)
		f_range[-1] = freqUB
		aic = [{} for f in range(len(f_range))]
		for f in range(len(f_range)):
			freq = f_range[f]
			p = purity
			if ploidy == 2:		#dbinom-exact probability
				t_ans1 = binom.pmf(round(d*freq),d,p/(2-p)) + np.finfo(np.double).tiny     #somatic LOH
				aic[f]["Somatic, LOH CNmut=1"] = 2-2*np.log(t_ans1)
				if varType != "somatic":
					t_ans2 = binom.pmf(round(d*freq),d,1/(2-p)) + np.finfo(np.double).tiny     #germline LOH
					aic[f]["Germline, LOH CNmut=1"] = 2-2*np.log(t_ans2)
			for i in range(ploidy):
				t_ans3 = binom.pmf(round(d*freq),d,((i+1)*p)/(2*(1-p)+ploidy*p)) + np.finfo(np.double).tiny        #somatic LOH high CN
				aic[f]["Somatic, CNmut=%i"%(i+1)] = 2-2*np.log(t_ans3)
				if varType != "somatic":
					t_ans4 = binom.pmf(round(d*freq),d,(1-p+(i+1)*p)/(2*(1-p)+ploidy*p)) + np.finfo(np.double).tiny    #germline LOH high CN
					aic[f]["Germline, CNmut=%i"%(i+1)] = 2-2*np.log(t_ans4)
		smallest = aic[0]["Somatic, CNmut=1"]
		for f in range(len(f_range)):
			for name in name_likelihood:
				if smallest > aic[f][name]:
					smallest = aic[f][name]

		D = 0.0
		for f in range(len(f_range)):
			for name in name_likelihood:
				D += np.exp(-0.5*(aic[f][name]-smallest))
		for f in range(len(f_range)):
			for name in name_likelihood:
				weight[j][name] += (np.exp(-0.5*(aic[f][name]-smallest))/D)
	return CCF, weight

def conf_interval(x,x2,num_mut,std_dev):
	x_err = ((np.array(x2)/num_mut) - (np.array(x)/num_mut)**2)**0.5 
	return x-std_dev*np.array(x_err)

def predicted_purity_from_CFF(vaf_list,depth_list,SNV_list,out_name,ploidy_list,LOH_thres,std_dev,W_thres,varType):
	all_purity = np.arange(0.01,1,0.01)
	CCF = [[] for x in range(len(all_purity))]
	weight = [[] for x in range(len(all_purity))]
	out = open(out_name+".txt","w")

	for pur in range(len(all_purity)):
		CCF[pur],weight[pur] = gray_box_accuracy(vaf_list,depth_list,all_purity[pur],ploidy_list,varType)

	gene_CCF_weight = {}
	for var in range(len(SNV_list)):
		gene_CCF_weight[SNV_list[var]] = []
	sum_CCF_weight = [0.0 for x in range(len(all_purity))]	#after removing germline no LOH and subclonal
	sum_CCF_weight_sq = [0.0 for x in range(len(all_purity))]
	Bef1_sum_CCF_weight = [0.0 for x in range(len(all_purity))]	#before removing any mutation
	Bef1_sum_CCF_weight_sq = [0.0 for x in range(len(all_purity))]
	Bef2_sum_CCF_weight = [0.0 for x in range(len(all_purity))]	#after removing germline no LOH
	Bef2_sum_CCF_weight_sq = [0.0 for x in range(len(all_purity))]

	for var in range(len(SNV_list)):
		name_likelihood = weight[0][var].keys()
		for pur in range(len(all_purity)):
			temp = 0.0
			for name in name_likelihood:
				temp += ((CCF[pur][var][name]-1)**2)*weight[pur][var][name]
			gene_CCF_weight[SNV_list[var]].append(np.log10(temp))
			Bef1_sum_CCF_weight[pur]+=temp
			Bef1_sum_CCF_weight_sq[pur] += temp**2

	num_color = len(SNV_list)
	cm = plt.get_cmap('gist_ncar')
	dic_color = [cm(1.*i/num_color) for i in range(0,num_color,1)]
	fig = plt.figure()
	lines = []
	labels = []
	s1 = plot(fig,1,all_purity,[],[],SNV_list,gene_CCF_weight,dic_color,SNV_list,"Variants likelihood distribution",[],[],[])
	lines.extend(s1[0])
	labels.extend(s1[1])

	##remove germline mutation no LOH
	index_purity_bef1 = np.argmin(Bef1_sum_CCF_weight)
	pur1 = all_purity[index_purity_bef1]
	if varType == "somatic":
		n_SNV, n_vaf, n_depth, n_ploidy = SNV_list, vaf_list, depth_list, ploidy_list
		Bef2_sum_CCF_weight = Bef1_sum_CCF_weight
		Bef2_sum_CCF_weight_sq = Bef1_sum_CCF_weight_sq
	else:
		germ_mut_index = []
		for mut in range(len(SNV_list)):
			germ_weight = 0.0
			flag_LOH = False
			if ploidy_list[mut] == 2:
			        germ_list = ["Germline, LOH CNmut=1"]
			else:
	                        germ_list = []
			for i in range(ploidy_list[mut]):
				germ_list.append("Germline, CNmut=%i"%(i+1))
			for model in germ_list:
				germ_weight += weight[index_purity_bef1][mut][model]
				if model == "Germline, LOH CNmut=1" or model == "Germline, CNmut=%i"%ploidy_list[mut]:
					if weight[index_purity_bef1][mut][model] > LOH_thres:
						flag_LOH = True
			if not flag_LOH and germ_weight > W_thres:
				germ_mut_index.append(mut)

		n_SNV, n_vaf, n_depth, n_ploidy = remove_mut(vaf_list,depth_list,SNV_list,ploidy_list,germ_mut_index)
		for new_var in n_SNV:
			unlog = np.power(10,gene_CCF_weight[new_var])
			for pur in range(len(all_purity)):
				Bef2_sum_CCF_weight[pur] += unlog[pur]
				Bef2_sum_CCF_weight_sq[pur] += unlog[pur]**2
		
	###remove subclonal mutation
	index_purity_bef2 = np.argmin(Bef2_sum_CCF_weight)
	pur2 = all_purity[index_purity_bef2]
	subclonal_mut_index = []
	for mut in range(len(n_SNV)):
		f = n_vaf[mut]
		d = n_depth[mut]
		Y = n_ploidy[mut]
		if binom.cdf(round(d*f),d,pur2/(2*(1-pur2)+Y*pur2)) < 0.01:      #cumulative probability
			subclonal_mut_index.append(mut)
	n2_SNV, n2_vaf, n2_depth, n2_ploidy = remove_mut(n_vaf,n_depth,n_SNV,n_ploidy,subclonal_mut_index)
	for new_var in n2_SNV:
		unlog = np.power(10,gene_CCF_weight[new_var])
		for pur in range(len(all_purity)):
			sum_CCF_weight[pur] += unlog[pur]
			sum_CCF_weight_sq[pur] += unlog[pur]**2

	pur3 = all_purity[np.argmin(sum_CCF_weight)]
	#expected vaf and expected models
	exp1 = expected_vaf_model(SNV_list,all_purity,index_purity_bef1,weight,ploidy_list,vaf_list)
	exp2 = expected_vaf_model(SNV_list,all_purity,index_purity_bef2,weight,ploidy_list,vaf_list)
	exp3 = expected_vaf_model(SNV_list,all_purity,np.argmin(sum_CCF_weight),weight,ploidy_list,vaf_list)

	all_mod = list(set(exp1[1])|set(exp2[1])|set(exp3[1]))
	oth_color = [cm(1.*i/len(all_mod)) for i in range(0,len(all_mod),1)]

	Bef1_sum_CCF_weight_LB = conf_interval(Bef1_sum_CCF_weight,Bef1_sum_CCF_weight_sq,len(SNV_list),std_dev)
	CI_p1 = all_purity[Bef1_sum_CCF_weight_LB < np.amin(Bef1_sum_CCF_weight)].tolist()
	Bef2_sum_CCF_weight_LB = conf_interval(Bef2_sum_CCF_weight,Bef2_sum_CCF_weight_sq,len(n_SNV),std_dev)
	CI_p2 = all_purity[Bef2_sum_CCF_weight_LB < np.amin(Bef2_sum_CCF_weight)].tolist()
	sum_CCF_weight_LB = conf_interval(sum_CCF_weight,sum_CCF_weight_sq,len(n2_SNV),std_dev)
	CI_p3 = all_purity[sum_CCF_weight_LB < np.amin(sum_CCF_weight)].tolist()
	if (len(CI_p1) > 0 and check_CI(CI_p1)) or (len(CI_p2) > 0 and check_CI(CI_p2)):
		out.write("There is discontinuous confidence interval around purity estimates before reaching the last purity.\n")
	if len(CI_p3) > 0 and check_CI(CI_p3):
		out.write("Warning! There is discontinuous confidence interval around final purity estimate; we may not have enough power to break ambiguity between different purity estimates.\n")

	out.write("len_SNV_list:%i\n"%len(SNV_list))
	out.write("len_SNV_list_after_germline_removed:%i\n"%len(n_SNV))
	out.write("len_SNV_list_after_germline_and_subclone_removed:%i\n"%len(n2_SNV))
	out.write("All_variants\tPresence(T)/Absence(F)_after_germline_removed\tPresence(T)/Absence(F)_after_germline&subclone_removed\n")
	for var in gene_CCF_weight.keys():
		out.write(var+"\t"+str(var in n_SNV)+"\t"+str(var in n2_SNV)+"\n")
	out.write("purity:\t"+str(round(pur1,2))+"\t"+",".join([str(round(x,2)) for x in CI_p1])+"\n")
	out.write("purity_after_removing_germline_no_LOH:\t"+str(round(pur2,2))+"\t"+",".join([str(round(x,2)) for x in CI_p2])+"\n")
	out.write("purity_after_removing_germline_no_LOH_and_subclonal:\t"+str(round(pur3,2))+"\t"+",".join([str(round(x,2)) for x in CI_p3])+"\n")
	out.write(out_name.split("/")[-1]+"\t"+str(round(pur3,2))+"\t"+",".join([str(round(x,2)) for x in CI_p3])+"\n")
	out.close()

	s2 = plot(fig,2,all_purity,Bef1_sum_CCF_weight,Bef1_sum_CCF_weight_LB,[],[],[],[],"Purity estimated from all variants",[],[],[])
	lines.extend(s2[0])
	labels.extend(s2[1])
	s3 = plot(fig,3,[],[],[],[],[],oth_color,[],"VAF distribution at "+r'$p$'+" = %.2f"%(round(pur1,2)),vaf_list,exp1,all_mod)
	lines.extend(s3[0])
	labels.extend(s3[1])
	s4 = plot(fig,4,[],[],[],[],[],[],[],"CCF distribution at "+r'$p$'+" = %.2f"%(round(pur1,2)),exp1[2],[],[])
	s5 = plot(fig,5,all_purity,[],[],n_SNV,gene_CCF_weight,dic_color,SNV_list,"Variants likelihood distribution after removing germline heterozygous mutations",[],[],[])
	s6 = plot(fig,6,all_purity,Bef2_sum_CCF_weight,Bef2_sum_CCF_weight_LB,[],[],[],[],"Purity estimated after excluding germline heterozygous mutations",[],[],[])
	s7 = plot(fig,7,[],[],[],[],[],oth_color,[],"VAF distribution at "+r'$p$'+" = %.2f"%(round(pur2,2)),vaf_list,exp2,all_mod)
	for l in range(len(s7[1])):
		if s7[1][l] not in labels:
			labels.append(s7[1][l])
			lines.append(s7[0][l])
	s8 = plot(fig,8,[],[],[],[],[],[],[],"CCF distribution at "+r'$p$'+" = %.2f"%(round(pur2,2)),exp2[2],[],[])
	s9 = plot(fig,9,all_purity,[],[],n2_SNV,gene_CCF_weight,dic_color,SNV_list,"Variants likelihood distribution after removing germline heterozygous and sub-clonal mutations",[],[],[])
	s10 = plot(fig,10,all_purity,sum_CCF_weight,sum_CCF_weight_LB,[],[],[],[],"Purity estimated after excluding germline heterozygous and sub-clonal mutations",[],[],[])
	s11 = plot(fig,11,[],[],[],[],[],oth_color,[],"VAF distribution at "+r'$p$'+" = %.2f"%(round(pur3,2)),vaf_list,exp3,all_mod)
	for l in range(len(s11[1])):
		if s11[1][l] not in labels:
			labels.append(s11[1][l])
			lines.append(s11[0][l])
	s12 = plot(fig,12,[],[],[],[],[],[],[],"CCF distribution at "+r'$p$'+" = %.2f"%(round(pur3,2)),exp3[2],[],[])
	fig.set_size_inches(20,12)
	plt.tight_layout()
	fig.savefig(out_name+"_likelihood.png")
#	fig.savefig(out_name+"_likelihood.eps",format="eps",dpi=350)
	plt.close(fig)

	##Plotting figure legend
	fig_leg = plt.figure()
	fig_leg.legend(lines, labels,'center',ncol=4)
	fig_leg.savefig(out_name+"_likelihood_legend.png",bbox_inches='tight')
#	fig_leg.savefig(out_name+"_likelihood_legend.eps",bbox_inches='tight',format="eps",dpi=350)
	plt.close(fig_leg)

	return str(all_purity[np.argmin(sum_CCF_weight)]),",".join([str(x) for x in CI_p3])

def expected_vaf_model(snv_list,all_purity,index,weight,ploidy_list,obs_vaf_list):
	vaf_list = []
	ccf_list = []
	model_list = []
	purity = all_purity[index]
	for snv in range(len(snv_list)):
		largest_w = 0.0
		best_model = ""
		name_likelihood = weight[index][snv].keys()
		for name in name_likelihood:
			if weight[index][snv][name] > largest_w:
				largest_w = weight[index][snv][name]
				best_model = name
		cmut = int(best_model.split("=")[1])
		tmp_Y = ploidy_list[snv]
		if len(best_model.split("LOH")) == 2:
			tmp_Y = 1
		if best_model.split(",")[0] == "Somatic":
			vaf = cmut*purity/((2*(1-purity))+(tmp_Y*purity))
			ccf = obs_vaf_list[snv]/vaf
		else:
			vaf = (1-purity+(cmut*purity))/((2*(1-purity))+(tmp_Y*purity))
			ccf = ((obs_vaf_list[snv]*((2*(1-purity))+(tmp_Y*purity)))-1+purity)/(cmut*purity)
		if vaf not in vaf_list:
			vaf_list.append(vaf)
			model_list.append(best_model+", ploidy="+str(tmp_Y))
		ccf_list.append(ccf)
	return [vaf_list,model_list,ccf_list]

def check_CI(CI_p):
	length = (CI_p[-1]*100)-(CI_p[0]*100)+1
	if len(CI_p) == length:
		return False
	else:
		return True     #indicate discontinuous confidence interval

def plot(fig,count,all_purity,sumCW,sumCW_LB,snv,gene_CCF_weight,dic_color,SNV_list,title_name,x_dist,exp_vaf_mod,all_model):
	lines = []
	labels = []
	legend_fontsize = 16
	panel_per_row = 4
	ax = fig.add_subplot(3,4,count)
	if count == 1:
		for each in range(len(snv)):
			lines.append(ax.plot(all_purity,gene_CCF_weight[snv[each]],label=snv[each],color=dic_color[SNV_list.index(snv[each])])[0])
			labels.append(snv[each])
	elif count == 2:
		lines.append(ax.plot(all_purity,np.log10(sumCW),label="sum_CCF_W",color='black')[0])
		lines.append(ax.plot(all_purity,np.log10(sumCW_LB),label='sum_CCF_W_LB',color='brown',linestyle='dashed')[0])
		lines.append(ax.axhline(y=np.log10(np.amin(sumCW)),label="min_sum_CCF_W",color='grey',linestyle='dotted'))
		labels.extend(["sum_CCF_W",'sum_CCF_W_LB',"min_sum_CCF_W"])
	elif count == 5 or count == 9:
		for each in range(len(snv)):
			ax.plot(all_purity,gene_CCF_weight[snv[each]],color=dic_color[SNV_list.index(snv[each])])
	elif count == 6 or count == 10:
		ax.plot(all_purity,np.log10(sumCW),color='black')
		ax.plot(all_purity,np.log10(sumCW_LB),color='brown',linestyle='dashed')
		ax.axhline(y=np.log10(np.amin(sumCW)),color='grey',linestyle='dotted')
	elif count%panel_per_row == 3 or count%panel_per_row == 0:
		sns.distplot(x_dist,color="lightgrey",kde=False,ax=ax)
		ax2 = ax.twinx()
		ax2.yaxis.set_ticks([])
		sns.kdeplot(x_dist,color="dimgrey",ax=ax2)
		if count%panel_per_row == 3:
			exp_vaf = exp_vaf_mod[0]
			exp_model = exp_vaf_mod[1]
			for v in range(len(exp_vaf)):
				lines.append(ax.axvline(exp_vaf[v],c=dic_color[all_model.index(exp_model[v])],linestyle='dashed',label=exp_model[v]))
				labels.append(exp_model[v])
			ax.set_xlabel("Observed variant allele frequencies",fontsize=legend_fontsize)
		else:
			ax.set_xlabel("Cancer cell fractions",fontsize=legend_fontsize)
		ax.set_ylabel("Counts",fontsize=legend_fontsize)
	ax.set_title("\n".join(textwrap.wrap(title_name, 35)),y=1.03,fontsize=18)
	if count%panel_per_row == 1:
		ax.set_title("\n".join(textwrap.wrap(title_name, 37)),y=1.03,fontsize=18)
		plt.xlabel("Purity",fontsize=legend_fontsize)
		plt.ylabel(r'$log_{10} L(p)$',fontsize=legend_fontsize)
	elif count%panel_per_row == 2:
		plt.xlabel("Purity",fontsize=legend_fontsize)
		plt.ylabel(r'$log_{10}\sum L(p)$',fontsize=legend_fontsize)
	return [lines,labels]

def remove_mut(vaf_list,depth_list,SNV_list,ploidy_list,index):
	n_SNV = []
	n_vaf = []
	n_depth = []
	n_ploidy = []
	for ind in range(len(SNV_list)):
		if ind not in index:
			n_SNV.append(SNV_list[ind])
			n_vaf.append(vaf_list[ind])
			n_depth.append(depth_list[ind])
			n_ploidy.append(ploidy_list[ind])
	return n_SNV, n_vaf, n_depth, n_ploidy

def main():
	parser = argparse.ArgumentParser(description="This is a program that estimates specimen purity from tumor-only sample sequenced with deep sequencing, called All-FIT (Allele Frequency based Imputation of Tumor Purity). We do not provide model for a variant that is being called on chrX (male) or chrY. Users may consider the option of removing all variants on chrX (male) and on chrY, if there is any in the input sample. Input sample should be a tab-delimited file with 4 columns and a header of ID\tAllele_Freq\tDepth\tPloidy")
	parser.add_argument('-i','--inputDirFile',help='Input file name with path',required=True)
	parser.add_argument('-d','--outputDir',help='Output directory',required=True)
	parser.add_argument('-o','--outputName',help='Output file name prefix',required=True)
	parser.add_argument('-s','--standardDeviation',default=2,help='How many standard deviation or confidence interval of estimated purity')
	parser.add_argument('-t','--typeOfInput',choices=['somatic','all'],default='all',help='Types of variants input whether germline variants are removed(somatic) or not(all)')
	args = parser.parse_args()

	varType = args.typeOfInput	
	std_dev = args.standardDeviation
	LOH_thres = 0.5
	W_thres = 0.7
	CCF_pred_purity = ""
	CI_pred_purity = ""
	ind_vaf = []
	ind_depth = []
	ind_SNV = []
	ind_ploidy = []
	out_name = args.outputName
	out_dir = args.outputDir
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	filename = args.inputDirFile
	fh = open(filename,"r")
	data = fh.readlines()[1:]
	fh.close()
	
	for row in data:
		row = row.strip().split("\t")
		ind_vaf.append(float(row[1])/100)
		ind_depth.append(int(row[2]))
		ind_SNV.append(row[0]+":%.2f"%ind_vaf[-1]+":"+row[2])
		ind_ploidy.append(int(row[3]))

	CCF_pred_purity,CI_pred_purity = predicted_purity_from_CFF(ind_vaf,ind_depth,ind_SNV,out_dir+"/"+out_name,ind_ploidy,LOH_thres,std_dev,W_thres,varType)
	print(out_name,"\t",CCF_pred_purity,"\t",CI_pred_purity)

if __name__ == "__main__":
	main()
