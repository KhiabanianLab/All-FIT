#!/usr/bin/python

import subprocess as sp, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import pandas as pd

loc=sp.Popen(['pwd'],stdout=sp.PIPE).communicate()[0].decode().strip()+"/"
input_file = sys.argv[1]
fz = open(input_file,"r")
ref = fz.readlines()[1:]
fz.close()
input_list = [[] for x in range(81)]	#storing file_num for each purity - SIM_DATA_ANS_3.xls
for row in ref:
	row = row.strip().split("\t")
	input_list[int(row[1])-10].append(row[0])
	
alpha = ["alpha0.%i"%i for i in range(1,10,1)]
alpha.append("alpha1.0")
trial = ["trial%i"%i for i in range(1,101,1)]
fsize = 16
x_label = ["10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-90"]
y_data = [[[] for y in range(len(alpha))] for x in range(len(x_label))]
for a in range(len(alpha)):
	print(alpha[a])
	y_group = [{} for x in range(len(input_list))]
	for i in range(len(input_list)):
		for name in input_list[i]:
			y_group[i][name] = []
	for t in range(len(trial)):
		fk = open(loc+"%s_%s_our_real_purity.txt"%(trial[t],alpha[a]),"r")
		data = fk.readlines()[1:]
		fk.close()
		for row in data:
			row = row.strip().split("\t")
			real_p = round(float(row[3]),2)
			est_p = round(float(row[1]),2)
			err_p = round(np.absolute(est_p-real_p),2)
			y_group[round(real_p*100)-10][row[0]].append(err_p)
	for i in range(len(input_list)):
		for name in input_list[i]:
			if i == 80:
				index = 7
			else:
				index = int(i/10)
			y_data[index][a].append(np.median(y_group[i][name]))
y_plot = []
for x in range(len(x_label)):
	for a in range(len(alpha)):
		for diff in y_data[x][a]:
			y_plot.append([diff,x_label[x],alpha[a]])

y_df = pd.DataFrame(np.array(y_plot),columns = ['diff','p','alpha'])
y_df["diff"]=y_df["diff"].astype(float)

#print(y_df.query('alpha == "alpha0.4" and p == "10-19"')["diff"].value_counts())

fig = plt.figure()
ax = plt.axes()
my_pal={"alpha0.1":"lightblue", "alpha0.2":"lightpink","alpha0.3":"darkorange","alpha0.4":"maroon","alpha0.5":"blue","alpha0.6":"darkmagenta","alpha0.7":"green","alpha0.8":"grey","alpha0.9":"mediumorchid","alpha1.0":"springgreen"}
ax = sns.violinplot(x="p",y="diff",data=y_df,hue="alpha",inner=None,scale="count",palette=my_pal)
#ax = sns.stripplot(x="p",y="diff",data=y_df,hue="alpha",size=2,jitter=0.05,color="black",dodge=True)
x_ticklabel = []
plt.xticks(np.arange(len(x_label)),x_label)
plt.xlabel("Simulated purity",fontsize=fsize)
plt.ylabel("Difference between estimated purity and simulated purity",fontsize=fsize)
ax.legend(["alpha = 0.1 (30x-100x)","alpha = 0.2 (60x-200x)","alpha = 0.3 (90x-300x)","alpha = 0.4 (120x-400x)","alpha = 0.5 (150x-500x)","alpha = 0.6 (180x-600x)","alpha = 0.7 (210x-700x)","alpha = 0.8 (240x-800x)","alpha = 0.9 (270x-900x)","alpha = 1.0 (300x-1000x)"],fontsize=fsize)
plt.setp(ax.get_xticklabels(), fontsize=fsize)
plt.setp(ax.get_yticklabels(), fontsize=fsize)
fig.set_size_inches(17,13)
fig.savefig(loc+"no_dot_err_dist_real_p.png")
fig.savefig(loc+"no_dot_err_dist_real_p.eps",format="eps",dpi=350)
plt.close(fig)
