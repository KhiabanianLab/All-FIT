#!/usr/bin/python

import subprocess as sp
import numpy as np, matplotlib
from scipy.stats import pearsonr
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

loc = sp.Popen(["pwd"],stdout=sp.PIPE).communicate()[0].strip()+"/"
#print loc

fh = open(loc+"our_Abs_real_purity.txt","r")
data = fh.readlines()[1:]
fh.close()

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
for row in data:
	row = row.strip().split("\t")
	real_p = round(float(row[4]),2)
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
#			print row[0]
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

print "scipy correlation between real purity and our method purity with p-value:",pearsonr(real_v,est_v)
threshold = 1e-323
if np.argmax([pearsonr(real_v,est_v)[1],threshold]) == 0:
	print "p-value > threshold (",threshold,")"
else:
	print "p-value < threshold (",threshold,")"
print "scipy correlation between real purity and Absolute purity with p-value:",pearsonr(real_v,Abs_v)
print "scipy correlation between our method purity and Absolute purity with p-value:",pearsonr(est_v,Abs_v)
print "for our method:"
print "accuracy within confidence region:",our_count*100.0/len(data)
print "count our method purity same as real purity:",our_count_exact
print "count our method purity higher than real purity and not in confidence region:",our_count_more
print "count our method purity lower than real purity and not in confidence region:",our_count_less
print "for Absolute:"
print "count Absolute purity same as real purity and accuracy:",Abs_count,Abs_count*100.0/len(data)
print "count Absolute purity higher than real purity:",Abs_count_more
print "count Absolute purity lower than real purity:",Abs_count_less


###Making plot
big_x = [x_value,x_value,z_value]
big_y = [y_value,z_value,y_value]
x_label = ["Simulated purity","Simulated purity"]
y_label = ["Estimated purity (All-FIT)","Estimated purity (ABSOLUTE)"]
legend_fontsize = 12

fig = plt.figure()
count = 0
cm = plt.cm.get_cmap('jet')
for k in range(0,2,1):
	count += 1
	ax = fig.add_subplot(1,2,count)
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
		#       ax.hlines(y_value[i],y_min[i],y_max[i])
		        ax.vlines(big_x[k][i],y_min[i],y_max[i],linewidth=0.1)
	cax = ax.scatter(x, y, s=10, c=z, edgecolor='', cmap=cm,vmin=0, vmax=0.016)
	fig.colorbar(cax)
	plt.xlabel(x_label[k],fontsize=legend_fontsize)
	plt.ylabel(y_label[k],fontsize=legend_fontsize)
	ax.set_aspect('equal')
	ax.set_xlim(0,105)
	ax.set_ylim(0,105)
fig.set_size_inches(12,5)
fig.savefig(loc+"accuracy_purity.png")
fig.savefig(loc+"accuracy_purity.eps",format="eps",dpi=350)
plt.close(fig)
