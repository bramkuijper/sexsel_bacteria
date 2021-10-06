#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

def get_data_start(filename, startcolname):

    f = open(filename)
    fl = f.readlines()

    for (idx, line in enumerate(fl)):
        if re.search("^" + startcolname, line) != None:
            return(idx)

    return(None)

# get the filename from the command line
filename = sys.argv[1]

data = pd.read_csv(
        filename,
        sep=";",
        skiprows=parline)

# old data instead of new
# we need to rename columns
if "q0" not in data.columns.values:
    data = rename_cols(data)

#row_each = math.floor(float(nrows)/1000)
#
## make dataset shorter so that we don't plot megabytes
#data = data.iloc[range(0,nrows,int(row_each)),:]

b1 = "#549CFF" 
r1 = "#FF5500"

# initialize and specify size 
fig = plt.figure(figsize=(10,13))

num_rows = 6

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting % type 1 offspring
ax1 = plt.subplot(num_rows,1,1)

ax1.plot(
        data["generation"],data["q0"],'b',
        data["generation"],data["q1"],'darkgreen',
        data["generation"],data["q2"],'black',
        )

ax1.set_ylabel(r'Prob. offspring is $z_{1}$')
ax1.legend((r'$q_{0}$',
                r'$q_{\mathrm{one signals}}$',
                r'$q_{\mathrm{both signal}}$'
                ))
ax1.set_ylim(-0.05,1.05)


# signals
plt.subplot(num_rows,1,2)

plt.plot(data["generation"],data["sm0"],'y',
        data["generation"],data["sm1"],'g',
        data["generation"],data["sf0"],'magenta',
        data["generation"],data["sf1"],'b',
        linewidth=1)
plt.legend((
            r'$s_{\mathrm{pat},e_{1}}$',
            r'$s_{\mathrm{pat},e_{2}}$',
            r'$s_{\mathrm{mat},e_{1}}$',
            r'$s_{\mathrm{mat},e_{2}}$'))

plt.ylabel(r'Class freqs')
plt.ylim(-0.05,1.05)

# class freqs
plt.subplot(num_rows,1,3)

plt.plot(data["generation"],data["u0"],'y',
        data["generation"],data["u1"],'g',
        data["generation"],data["u2"],'magenta',
        data["generation"],data["u3"],'b',
        linewidth=1)
plt.legend((r'$u_{1}$',r'$u_{2}$',r'$u_{3}$',r'$u_{4}$'))

plt.ylabel(r'Class freqs')
plt.ylim(-0.05,1.05)


# add 3rd subplot depicting relatedness
plt.subplot(num_rows,1,4)

plt.plot(data["generation"],data["v0"],'y',
        data["generation"],data["v1"],'g',
        data["generation"],data["v2"],'magenta',
        data["generation"],data["v3"],'b',
        data["generation"],data["ev"],'k',
        linewidth=1)
plt.legend((r'$v_{1}$',r'$v_{2}$',r'$v_{3}$',r'$v_{4}$',r'$\lambda$'))
plt.ylabel(r'Reproductive values')

# relatedness
plt.subplot(num_rows,1,5)

plt.plot(data["generation"],data["Qmm0"],'y',
        data["generation"],data["Qmm1"],'g',
        data["generation"],data["Qfm0"],'b',
        data["generation"],data["Qfm1"],'r',
        data["generation"],data["Qff0"],'magenta',
        data["generation"],data["Qff1"],'k',
        linewidth=1)
plt.legend(
        (
            r'$Q_{\mathrm{mm},1}$',
            r'$Q_{\mathrm{mm},2}$',
            r'$Q_{\mathrm{fm},1}$',
            r'$Q_{\mathrm{fm},2}$',
            r'$Q_{\mathrm{ff},1}$',
            r'$Q_{\mathrm{ff},2}$'))

plt.ylabel(r'Relatedness')
plt.ylim(-0.05,1.05)

# prob z1|ei
plt.subplot(num_rows,1,6)

data["z1_e1"] = data["q0"] * (1-data["sm0"]) * (1-data["sf0"]) + data["q1"] * ((1-data["sf0"]) * data["sm0"] + (1-data["sm0"]) *data["sf0"]) + data["q2"] * data["sf0"] * data["sm0"]

data["z1_e2"] = data["q0"] * (1-data["sm1"]) * (1-data["sf1"]) + data["q1"] * ((1-data["sf1"]) * data["sm1"] + (1-data["sm1"]) *data["sf1"]) + data["q2"] * data["sf1"] * data["sm1"]


plt.plot(data["generation"],data["z1_e1"],'b',
        data["generation"],data["z1_e2"],'r',
        linewidth=1)

plt.legend((r'$\mathrm{Pr}\left(z1\mid e1\right)$',
            r'$\mathrm{Pr}\left(z1\mid e2\right)$'))

plt.ylabel(r'Class freqs')
plt.ylim(-0.05,1.05)

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
