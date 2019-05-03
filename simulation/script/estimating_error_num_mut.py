#!/usr/bin/python

import subprocess as sp, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

loc=sp.Popen(['pwd'],stdout=sp.PIPE).communicate()[0].decode().strip()+"/"

num_mut = []
diff_real_est = []

fh = open(loc+"num_mut.txt","r")
ref = fh.readlines()
fh.close()
for row in ref:
	row = row.strip().split("\t")
	num_mut.append(int(row[1]))

fk = open(loc+"our_real_purity.txt","r")
data = fk.readlines()[1:]
fk.close()
for row in data:
	row = row.strip().split("\t")
	real_p = round(float(row[3]),2)
	est_p = round(float(row[1]),2)
	err_p = round(est_p-real_p,2)
	diff_real_est.append(err_p)

x_label = ["5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95-99","100"]

####for violin plot#############################
y_group = [[] for x in range(len(x_label))]
for each in range(len(num_mut)):
	index = int(num_mut[each]/5)-1
	y_group[index].append(diff_real_est[each])
for each in range(len(x_label)):
	print(x_label[each],"\t",len(y_group[each]),"\t",np.std(y_group[each]))
###############################################

fig = plt.figure()
ax = plt.axes()
ax = sns.violinplot(data=y_group,inner=None,scale="count",color=".8")
ax = sns.stripplot(data=y_group,size=2,jitter=True)
plt.xticks(np.arange(len(x_label)),x_label,rotation=45,horizontalalignment='right')
plt.xlabel("Intervals of number of mutations")
plt.ylabel("Difference between estimated purity and simulated purity")
fig.set_size_inches(14,8)
fig.savefig(loc+"err_dist_num_mut.png")
fig.savefig(loc+"err_dist_num_mut.eps",format="eps",dpi=350)
plt.close(fig)

