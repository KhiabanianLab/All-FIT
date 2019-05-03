#!/usr/bin/python

import subprocess as sp
import numpy as np, matplotlib
from scipy.stats import pearsonr
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy import stats

def conf_int(x,y):
	r = pearsonr(x,y)[0]
	r_z = np.arctanh(r)
	se = 1/np.sqrt(len(x)-3)
	alpha = 0.05
	z = stats.norm.ppf(1-alpha/2)
	lo_z, hi_z = r_z-z*se, r_z+z*se
	lo, hi = np.tanh((lo_z, hi_z))
	return r,lo,hi

loc = sp.Popen(["pwd"],stdout=sp.PIPE).communicate()[0].strip()+"/"
#print loc

fk = open(loc+"set6_output","r")
ref = fk.readlines()
fk.close()
step = 25
#step = 10
highP = [[] for x in range(100/step)]	#[0-25),[25-50),[50-75),[75-]
highIndex = [[] for x in range(100/step)]
for row in range(len(ref)):
	num_type = np.array([int(x.strip()) for x in ref[row].split("[")[1].split("]")[0].split(",")])
	high_count = num_type[[6,7]].sum()
        num_mut = num_type.sum()
	percent = high_count*100.0/num_mut
	category = int(percent/step)
	highIndex[category].append(row)
	highP[category].append(percent)


fh = open(loc+"our_Abs_real_purity.txt","r")
#fh = open(loc+"our_method-removed-germline_real_purity.txt","r")
data = fh.readlines()[1:]
fh.close()

our_corr = [[],[],[]]   #[r,lower_confidence, upper_confidence]
Abs_corr = [[],[],[]]
inter = []
fsize = 18
for interval in range(len(highP)):
	if len(highIndex[interval]) == 0:
		continue
	#variables for making plot
	x_value = []
	y_value = []
	y_max = []
	y_min = []
	z_value = []
	
	our_count = 0
	our_count_less = 0
	our_count_more = 0
	our_count_exact = 0
	Abs_count_less = 0
	Abs_count_more = 0
	Abs_count = 0
	real_v = []
	est_v = []
	Abs_v = []
	for line in highIndex[interval]:
		row = data[line]
		row = row.strip().split("\t")
		real_p = round(float(row[-1]),2)
		real_v.append(real_p)
	        est_p = round(float(row[1]),2)
		est_v.append(est_p)
		Abs_p = round(float(row[3]),2)
		Abs_v.append(Abs_p)
		x_value.append(int(round(real_p*100)))
		y_value.append(int(round(est_p*100)))
		z_value.append(int(round(Abs_p*100)))
		if real_p == Abs_p:
			Abs_count += 1
	#		print row[0]	#uncomment when printing list of files with Absolute_good_pattern
		elif Abs_p > real_p:
			Abs_count_more += 1
		else:
			Abs_count_less += 1
		if row[2].split(",")[0] == '' :
			print row[0],"does not have purity range"
			if real_p == est_p:
				our_count += 1
				our_count_exact += 1
			elif est_p > real_p:
				our_count_more += 1
			else:
				our_count_less += 1
			y_max.append(int(round(est_p*100)))
			y_min.append(int(round(est_p*100)))
			continue
	        p_range = [round(float(x),2) for x in row[2].split(",")]
		if real_p in p_range:
			our_count += 1
	#		print row[0]	#uncomment when printing list of files with our_method_good_pattern
	        else:
	                if est_p > real_p:
	                        our_count_more += 1
	                elif est_p < real_p:
	                        our_count_less += 1
	        if real_p == est_p:
	                our_count_exact += 1
		y_max.append(int(round(p_range[-1]*100)))
		y_min.append(int(round(p_range[0]*100)))

	if step*0.01*(interval+1) < 1.1:
		our_r_conf = conf_int(real_v,est_v)
		Abs_r_conf = conf_int(real_v,Abs_v)
		inter.append(step*0.01*(interval+1))
		our_corr[0].append(our_r_conf[0])
		our_corr[1].append(our_r_conf[0]-our_r_conf[1])
		our_corr[2].append(our_r_conf[2]-our_r_conf[0])
		Abs_corr[0].append(Abs_r_conf[0])
		Abs_corr[1].append(Abs_r_conf[0]-Abs_r_conf[1])
		Abs_corr[2].append(Abs_r_conf[2]-Abs_r_conf[0])

#	print "INTERVAL",interval	
#	print "scipy correlation between real purity and our method purity with p-value:",pearsonr(real_v,est_v)
#	threshold = 1e-323
#	if np.argmax([pearsonr(real_v,est_v)[1],threshold]) == 0:
#		print "p-value > threshold (",threshold,")"
#	else:
#		print "p-value < threshold (",threshold,")"
#	print "scipy correlation between real purity and Absolute purity with p-value:",pearsonr(real_v,Abs_v)
#	print "scipy correlation between our method purity and Absolute purity with p-value:",pearsonr(est_v,Abs_v)
#	print "for our method:"
#	print "accuracy within confidence region:",our_count*100.0/len(highIndex[interval])
#	print "count our method purity same as real purity:",our_count_exact
#	print "count our method purity higher than real purity and not in confidence region:",our_count_more
#	print "count our method purity lower than real purity and not in confidence region:",our_count_less
#	print "for Absolute:"
#	print "count Absolute purity same as real purity and accuracy:",Abs_count,Abs_count*100.0/len(data)
#	print "count Absolute purity higher than real purity:",Abs_count_more
#	print "count Absolute purity lower than real purity:",Abs_count_less
#	
#
	##Making plot
	big_x = [x_value,x_value]
	big_y = [y_value,z_value]
	x_label = ["Simulated purity","Simulated purity"]
	y_label = ["Estimated purity (All-FIT)","Estimated purity (ABSOLUTE)"]
	
	fig = plt.figure()
	cm = plt.cm.get_cmap('jet')
	count = 0
	for k in range(0,2,1):
		count += 1
		ax = fig.add_subplot(2,1,count)
		x = np.array(big_x[k])
		y = np.array(big_y[k])
		# Calculate the point density
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		
		# Sort the points by density, so that the densest points are plotted last
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]
	
		if k == 0:	
			for i in range(len(y)):
			        ax.vlines(big_x[k][i],y_min[i],y_max[i],linewidth=0.1)
			cax = ax.scatter(x, y, c=z, s=10, edgecolor='',cmap=cm, vmin=0, vmax=0.016)
		else:
			cax = ax.scatter(x, y, c=z, s=10, edgecolor='',cmap=cm, vmin=0, vmax=0.0012)
		cb=fig.colorbar(cax)
		cb.ax.tick_params(labelsize=fsize)
		plt.xlabel(x_label[k],fontsize=fsize)
		plt.ylabel(y_label[k],fontsize=fsize)
		ax.set_aspect('equal')
		ax.set_xlim(0,105)
		ax.set_ylim(0,105)
		plt.setp(ax.get_xticklabels(), fontsize=fsize)
		plt.setp(ax.get_yticklabels(), fontsize=fsize)
	fig.set_size_inches(7,11)
	fig.savefig(loc+"combined_accuracy_purity_interval%i.png"%interval)
	fig.savefig(loc+"combined_accuracy_purity_interval%i.eps"%interval,format="eps",dpi=350)
	plt.close(fig)

##plot correlation graph
#xticklabel = ["[0.0-%.1f)"%inter[0]]
#for i in range(0,len(inter)-1,1):
#	xticklabel.append("[%.1f-%.1f)"%(inter[i],inter[i+1]))
#
#fig = plt.figure()
#ax = plt.axes()
#ax.set_ylim(bottom=0.0,top=1.1)
#ax.errorbar(inter,our_corr[0],yerr = [our_corr[1],our_corr[2]],color='b',label="All-FIT",capsize=3)
#ax.errorbar(inter,Abs_corr[0],yerr = [Abs_corr[1],Abs_corr[2]],color='r',label="ABSOLUTE",capsize=3)
## get handles
#handles, labels = ax.get_legend_handles_labels()
## remove the errorbars
#handles = [h[0] for h in handles]
## use them in the legend
#ax.legend(handles, labels, loc='lower left',numpoints=1,fontsize=fsize)
#plt.xticks(inter,xticklabel)
#plt.setp(ax.get_xticklabels(), rotation=20,fontsize=fsize)
#plt.setp(ax.get_yticklabels(), fontsize=fsize)
#plt.xlabel("Fraction of mutations with ploidy > 2",fontsize=fsize)
#plt.ylabel("Pearson's "+r'$r$'+" between estimated and simulated purity",fontsize=fsize)
#fig.set_size_inches(9,9)
#fig.savefig("correlation.png")
#fig.savefig("correlation.eps",format='eps',dpi=350)
#plt.close(fig)

