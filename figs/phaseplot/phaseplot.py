#!/usr/bin/env python3

import pandas as pd
import numpy as np
import math
import re

from pythontools import multipanel
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator


fontpath = "/System/Library/Fonts/Supplemental/" 

rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

import matplotlib as mpl
mpl.use("pgf")

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
        r"\usepackage{units}",         # load additional packages
        r"\usepackage{mathspec}",         # load additional packages
        r"\setmainfont[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
        r"\setmathsfont(Digits,Latin,Greek)[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
        r"\setmathrm[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
         ])
}

mpl.rcParams.update(pgf_with_custom_preamble)

# find the initial line of a file
def find_init_line(filename):
    f = open(filename)
    fl = f.readlines()
    f.close()

    for i, line_i in enumerate(fl):
        if re.search("^time",line_i) != None:
            return(i)

    return(0)





# read in the data and select subset
the_data = pd.read_csv(
        filepath_or_buffer="summary_environmental_signaling.csv"
        ,sep=";"
        ,dtype={"file":object})

the_data = the_data.astype(dtype={"file":"string"})

subset = the_data[
        (the_data["cp"] > 0.0)
        & (the_data["mu_t2"] > 0.0)
        & (the_data["p2_t0"].isin([0.02,0.98])
            | the_data["t2_t0"].isin([0.02,0.98]))
        ]

print(subset["p2"])

nrows_subset = subset.shape[0]
assert(nrows_subset > 0)

# start the figure
the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1]
        ,filename="phase_plot_single.pdf"
        ,hspace=0.3
        ,width=6
        ,height=6
        )

# initialize the panel axes
the_axis = the_fig.start_block(
    row=0
    ,col=0)

for row_i in range(0,nrows_subset):
    the_file = str(subset["file"].values[row_i])
    print(the_file)

    dat = pd.read_csv(
            filepath_or_buffer=the_file
            ,skiprows=find_init_line(the_file)
            ,sep=";"
            )

    print(dat.columns.values)

    the_axis.plot(dat["p2"]
            ,dat["t2"]
            ,color="blue")

# end the figure
the_fig.end_block(
        the_axis
#        ,xlim=xlim
#        ,ylim=ylim
#        ,y_ticks_minor = 5
#        ,x_ticks_minor = 1
#        ,x_ticks_major_multiple =1
#        ,y_ticks_major_multiple =0.25 
#        ,xticks=row == the_fig.rows - 1
#        ,yticks=col==0
#        ,title=title
        ,xlabel=""
        ,ylabel=""
#        ,loc_title=True
#        ,loc_title_pos=[-0.05,1.05]
        )

the_fig.close(tight=True)
