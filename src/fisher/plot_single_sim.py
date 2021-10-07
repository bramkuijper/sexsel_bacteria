#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
fig = plt.figure(figsize=(10,15))

num_rows = 8

print(data.head())

row_ctr = 0

# subplot for susceptibles
gs = gridspec.GridSpec(
            nrows = 8
            ,ncols = 1)

ax1 = plt.subplot(gs[row_ctr, 0])

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


row_ctr += 1

for row_idx, plasmid_genotype in enumerate(["t1p1","t2p1","t1p2","t2p2"]):

    # subplot for infected
    ax1 = plt.subplot(gs[row_ctr, 0])

    for t_idx in range(1,3):
        for p_idx in range(1,3):
            ax1.plot(
                    data["time"]
                    ,data[f"It{t_idx}p{p_idx}{plasmid_genotype}"]
                    ,label="$I_{t_{" + str(t_idx) + "}p_{" + str(p_idx) + "}" +\
                            plasmid_genotype + "}$"
                    )

    row_ctr += 1
    ax1.set_ylabel(r'Infected')
    ax1.legend()


ax1 = plt.subplot(gs[row_ctr, 0])
ax1.plot(
    data["time"]
    ,data["S"]
    ,label="S")

ax1.plot(
    data["time"]
    ,data["I"]
    ,label="I")

ax1.plot(
    data["time"]
    ,data["N"]
    ,label="N")

ax1.set_ylabel(r'Totals')
ax1.legend()

row_ctr += 1

ax1 = plt.subplot(gs[row_ctr, 0])
ax1.plot(
    data["time"]
    ,data["p2"]
    ,label="$p_{2}$")

ax1.plot(
    data["time"]
    ,data["t2"]
    ,label="$t_{2}$")

ax1.legend()

row_ctr += 1

ax1 = plt.subplot(gs[row_ctr, 0])
ax1.plot(
    data["time"]
    ,data["p2_chr"]
    ,label="$p_{2,chr}$")

ax1.plot(
    data["time"]
    ,data["p2_plm"]
    ,label="$p_{2,plm}$")

ax1.plot(
    data["time"]
    ,data["t2_chr"]
    ,label="$t_{2,chr}$")
ax1.plot(
    data["time"]
    ,data["t2_plm"]
    ,label="$t_{2,plm}$")

ax1.legend()

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
