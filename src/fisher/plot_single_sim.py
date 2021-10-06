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

    for idx, line in enumerate(fl):
        if re.search("^" + startcolname, line) != None:
            return(idx)

    return(None)

# get the filename from the command line
filename = sys.argv[1]

parline = get_data_start(filename,"time")

data = pd.read_csv(
        filename,
        sep=";",
        skiprows=parline)

# initialize and specify size 
fig = plt.figure(figsize=(10,13))

num_rows = 6

# add first subplot depicting % type 1 offspring
ax1 = plt.subplot(num_rows,1,1)

# plot the susceptibles
for t_idx in range(1,3):
    for p_idx in range(1,3):
        ax1.plot(
                data["time"]
                ,data[f"St{t_idx}p{p_idx}"]
                ,label="$S_{t_{" + str(t_idx) + "}p_{" + str(p_idx) + "}}$"
                )

ax1.set_ylabel(r'Susceptible')
ax1.legend()
ax1.set_ylim(-0.05,1.05)

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
