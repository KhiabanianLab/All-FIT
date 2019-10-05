#!/usr/bin/python

import subprocess as sp, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

loc=sp.Popen(['pwd'],stdout=sp.PIPE).communicate()[0].decode().strip()+"/"

num_mut = []
fh = open(loc+"num_mut.txt","r")
ref = fh.readlines()
fh.close()
for row in ref:
	row = row.strip().split("\t")
	num_mut.append(int(row[1]))

p1 = sp.Popen(['ls',loc],stdout=sp.PIPE)
input_files = sp.Popen(['grep','purity.txt'],stdin=p1.stdout,stdout=sp.PIPE).communicate()[0].decode().strip().split("\n")

x_label = ["2-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95-99"]
y_group = [{} for x in range(len(x_label))]
y_data = [[] for x in range(len(x_label))]

for inp in input_files:
	fk = open(loc+inp,"r")
	data = fk.readlines()[101:]
	fk.close()

	diff_real_est = []
	for each in range(len(x_label)):
		y_group[each][inp] = []
	for row in data:
		row = row.strip().split("\t")
		real_p = round(float(row[3]),2)
		est_p = round(float(row[1]),2)
		err_p = round(np.absolute(est_p-real_p),2)
		diff_real_est.append(err_p)
	
	####for violin plot#############################
	for each in range(len(num_mut)):
		index = int(num_mut[each]/5)
		y_group[index][inp].append(diff_real_est[each])

for each in range(len(x_label)):
	for inp in input_files:
		y_data[each].append(np.median(y_group[each][inp]))

fsize = 16
fig = plt.figure()
ax = plt.axes()
ax = sns.violinplot(data=y_data,inner=None,scale="count",color=".8")
ax = sns.stripplot(data=y_data,size=2,jitter=True)
plt.xticks(np.arange(len(x_label)),x_label,rotation=45,horizontalalignment='right')
plt.xlabel("Number of mutations",fontsize=fsize)
plt.ylabel("Difference between estimated purity and simulated purity",fontsize=fsize)
plt.setp(ax.get_xticklabels(), fontsize=fsize)
plt.setp(ax.get_yticklabels(), fontsize=fsize)
#plt.yticks(np.arange(0, 100, 10))
fig.set_size_inches(18,13)
fig.savefig(loc+"err_dist_num_mut.png")
fig.savefig(loc+"err_dist_num_mut.eps",format="eps",dpi=350)
plt.close(fig)
