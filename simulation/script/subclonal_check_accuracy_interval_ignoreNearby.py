#!/usr/bin/python

import subprocess as sp, sys
import numpy as np, matplotlib
from scipy.stats import pearsonr
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy import stats
import textwrap

def conf_int(x,y):
	r = pearsonr(x,y)[0]
	r_z = np.arctanh(r)
	se = 1/np.sqrt(len(x)-3)
	alpha = 0.05
	z = stats.norm.ppf(1-alpha/2)
	lo_z, hi_z = r_z-z*se, r_z+z*se
	lo, hi = np.tanh((lo_z, hi_z))
	return r,lo,hi

loc = sp.Popen(["pwd"],stdout=sp.PIPE).communicate()[0].decode().strip()+"/"
method = sys.argv[1]

fk = open(loc+"all_pattern","r")
ref = fk.readlines()
fk.close()
step = 25
#step = 10
subclonalP = [[] for x in range(int(100/step))]	#[0-25),[25-50),[50-75),[75-]
subclonalIndex = [[] for x in range(int(100/step))]
for row in range(len(ref)):
	sub_count = int(ref[row].strip().split("\t")[1])
	num_type = np.array([int(x.strip()) for x in ref[row].split("[")[1].split("]")[0].split(",")])
#	if method == "our_method-removed-germline":
#		num_mut = num_type[[0,1,2,6,8]].sum()   #counting somatic mutations only
#	else:
#	        num_mut = num_type.sum()
	num_mut = num_type.sum()
	percent = sub_count*100.0/num_mut
	category = int(percent/step)
	subclonalIndex[category].append(row)
	subclonalP[category].append(percent)


if method == "our_method-removed-germline":
	fh = open(loc+"our_method-removed-germline_Abs_real_purity.txt","r")	#for somatic mutations only
else:
	fh = open(loc+"our_Abs_real_purity.txt","r")
data = fh.readlines()[1:]
fh.close()
out_dir = loc+method+"_subclonal_analysis/"

our_corr = [[],[],[]]	#[r,lower_confidence, upper_confidence]
Abs_corr = [[],[],[]]
inter = []
fsize = 14
for interval in range(len(subclonalP)):
	if len(subclonalIndex[interval]) == 0:
		continue

	if interval == 3:
		multiply = 10.0
	else:
		multiply = 1000.0
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
	for line in subclonalIndex[interval]:
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
			print(row[0],"does not have purity range")
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

#	if step*0.01*(interval+1) < 0.8:
	our_r_conf = conf_int(real_v,est_v)
	Abs_r_conf = conf_int(real_v,Abs_v)
	inter.append(step*0.01*(interval+1))
	our_corr[0].append(our_r_conf[0])
	our_corr[1].append(our_r_conf[0]-our_r_conf[1])
	our_corr[2].append(our_r_conf[2]-our_r_conf[0])
	Abs_corr[0].append(Abs_r_conf[0])
	Abs_corr[1].append(Abs_r_conf[0]-Abs_r_conf[1])
	Abs_corr[2].append(Abs_r_conf[2]-Abs_r_conf[0])

	print("INTERVAL",interval)
	print("scipy correlation between real purity and our method purity with p-value:",pearsonr(real_v,est_v))
	threshold = 1e-323
	if np.argmax([pearsonr(real_v,est_v)[1],threshold]) == 0:
		print("p-value > threshold (",threshold,")")
	else:
		print("p-value < threshold (",threshold,")")
	print("scipy correlation between real purity and Absolute purity with p-value:",pearsonr(real_v,Abs_v))
	print("scipy correlation between our method purity and Absolute purity with p-value:",pearsonr(est_v,Abs_v))
	print("for our method:")
	print("accuracy within confidence region:",our_count*100.0/len(subclonalIndex[interval]))
	print("count our method purity same as real purity:",our_count_exact)
	print("count our method purity higher than real purity and not in confidence region:",our_count_more)
	print("count our method purity lower than real purity and not in confidence region:",our_count_less)
	print("for Absolute:")
	print("count Absolute purity same as real purity and accuracy:",Abs_count,Abs_count*100.0/len(subclonalIndex[interval]))
	print("count Absolute purity higher than real purity:",Abs_count_more)
	print("count Absolute purity lower than real purity:",Abs_count_less)
      
      
	#Making plot
	big_x = [x_value,x_value]
	big_y = [y_value,z_value]
	x_label = ["Simulated purity","Simulated purity"]
	y_label = ["Estimated purity (All-FIT)","Estimated purity (ABSOLUTE)"]
	
	fig = plt.figure()
	cm = plt.cm.get_cmap('jet')
	count = 0
	for k in range(0,2,1):
#	for k in range(0,1,1):
		count += 1
		ax = fig.add_subplot(2,1,count)
#		ax = plt.axes()
		x = np.array(big_x[k])
		y = np.array(big_y[k])

		# Calculate the point density
		xy_count = np.zeros(shape=(101,101))
		xy_index = [[[] for k in range(101)] for k in range(101)]
		for ind in range(len(big_x[k])):
#			print(k,ind,big_y[k][ind],big_x[k][ind])
			xy_count[big_y[k][ind]][big_x[k][ind]] += 1
			xy_index[big_y[k][ind]][big_x[k][ind]].append(ind)

		sorted_xy_idx = np.dstack(np.unravel_index(np.argsort(xy_count.ravel()), (101, 101)))[0]
		idx = []
		z = [0 for k in range(len(big_x[k]))]
		for each in sorted_xy_idx:
			x_idx = each[1]
			y_idx = each[0]
			idx.extend(xy_index[y_idx][x_idx])
			for ind in xy_index[y_idx][x_idx]:
				z[ind] = xy_count[y_idx][x_idx]*multiply/len(subclonalIndex[interval])
		idx = np.array(idx)
		z = np.array(z)

		# Sort the points by density, so that the densest points are plotted last
		x, y, z = x[idx], y[idx], z[idx]
	
		if k == 0:	
			for i in range(len(y)):
			        ax.vlines(big_x[k][i],y_min[i],y_max[i],linewidth=0.1)
			cax = ax.scatter(x, y, c=z, s=10, edgecolor='',cmap=cm, vmin=0, vmax=22.5)
		else:
			cax = ax.scatter(x, y, c=z, s=10, edgecolor='',cmap=cm, vmin=0, vmax=22.5)
		cb=fig.colorbar(cax)
		cb.ax.tick_params(labelsize=fsize)
		cb.ax.set_xlabel("\n".join(textwrap.wrap('N per coord per %i'%multiply, 12)), fontsize = fsize)
		cb.ax.xaxis.set_label_coords(2.50, -0.032)
		plt.xlabel(x_label[k],fontsize=fsize)
		plt.ylabel(y_label[k],fontsize=fsize)
		ax.set_aspect('equal')
		ax.set_xlim(0,105)
		ax.set_ylim(0,105)
		plt.setp(ax.get_xticklabels(), fontsize=fsize)
		plt.setp(ax.get_yticklabels(), fontsize=fsize)
	fig.set_size_inches(6,10)
	fig.savefig(out_dir+"combined_accuracy_purity_interval%i.png"%interval)
	fig.savefig(out_dir+"combined_accuracy_purity_interval%i.eps"%interval,format="eps",dpi=350)
	plt.close(fig)


#plot correlation graph
xticklabel = ["[0.0-%.2f)\nN = %i"%(inter[0],len(subclonalIndex[0]))]
for i in range(0,len(inter)-1,1):
	xticklabel.append("[%.2f-%.2f)\nN = %i"%(inter[i],inter[i+1],len(subclonalIndex[i+1])))

fig = plt.figure()
ax = plt.axes()
ax.set_ylim(bottom=0.0,top=1.1)
ax.errorbar(inter,our_corr[0],yerr = [our_corr[1],our_corr[2]],color='b',label="All-FIT",capsize=3)
ax.errorbar(inter,Abs_corr[0],yerr = [Abs_corr[1],Abs_corr[2]],color='r',label="ABSOLUTE",capsize=3)
# get handles
handles, labels = ax.get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
plt.rcParams["legend.labelspacing"] = 0.25
leg = ax.legend(handles, labels, loc='lower left',numpoints=1,fontsize=fsize)
leg.get_frame().set_linewidth(0.0)
leg.get_frame().set_facecolor('none')
plt.xticks(inter,xticklabel)
plt.setp(ax.get_xticklabels(), fontsize=fsize)
plt.setp(ax.get_yticklabels(), fontsize=fsize)
plt.xlabel("Percentage of sub-clonal mutations",fontsize=fsize)
plt.ylabel("Pearson's correlation between estimated and simulated purity",fontsize=fsize)
fig.set_size_inches(9,8)
fig.savefig(out_dir+"correlation.png")
fig.savefig(out_dir+"correlation.eps",format='eps',dpi=350)
plt.close(fig)
