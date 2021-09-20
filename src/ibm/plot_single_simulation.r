#!/usr/bin/env Rscript


# this script should be run from the command line rather than from within Rstudio
# to run it, go to the folder where it is saved
# and then try
# ./plot_single_simulation.r output_file.csv

# get rid of everything there
rm(list=ls())
# ggplot2 graphics library
library("ggplot2")

# some extra functionality to 
library("gridExtra")

# get command line arguments
#args = commandArgs(trailingOnly=TRUE)

# give an error message if you do not provide it with a simulation file name
#if (length(args) < 1)
#{
#    print("Usage: /path/plot_single_simulation.r output_file_blabla.csv")
#    stop()
#}

# list of simulation output files
files <- list.files(
        path=".",
        pattern="^sim_bact_sexsel.*?\\.csv$")

print(files)
# output files of simulation files
# consist of two parts
# 1. the actual data per timestep
# 2. parameters

# or (when order is reversed):
# 1. parameters
# 2. the actual data per timestep


# now some functions to find out where parameters are located in the file

# in case parameter listing is at the start of the file
# find out where the data starts
find_out_data_start <- function(filename) 
{
    # read in file line by line
    f <- readLines(filename)

    # loop through each line and detect a line starting with
    # the words 'time' or 'generation', as this is the header
    # after which the simulation output for each timestep is
    # printed
    for (line_i in seq(1,length(f)))
    {
        if (length(grep("^(time|generation)",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
} # end find_out_data_start()

# in case the 
# find out where the parameter listing starts
# so that we can read in the data part of the file 
# without having it messed up by the subsequent parameter listing
find_out_param_line <- function(filename) 
{
    # read in file line by line
    f <- readLines(filename)

    # make a reverse sequence of all the line numbers
    # so that we can go through all the lines in reverse
    # order (starting at the end, rather than at the beginning
    # allows us to find the end of the data
    # quicker.
    seqq <- seq(length(f),1,-1)

    # go through each line in the data file and find first line
    # where data is printed (i.e., a line which starts with a number)
    for (line_i in seqq)
    {
        # find the first line starting with a number
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }

    return(NA)
} # end find_out_param_line

data.first <- F
the.data <-NULL
the.base.name <- NULL

for (f_idx in 1:length(files)) 
{
	if (data.first)
	{
	    parameter_row <- find_out_param_line(files[f_idx])

	    if (is.na(parameter_row))
	    {
		print("cannot find data...")
		stop()
	    }

	    the.data <- read.table(files[f_idx], header=T, nrow=parameter_row - 1, sep=";")
	} else
	{
	    # read in data frame of corresponding simulation
	    parameter_row <- find_out_data_start(files[f_idx])
            the.base.name <- basename(files[f_idx])
            the.data <- read.table(files[f_idx], header=T, skip=parameter_row - 1, sep=";")
	    if(max(the.data$time) > 1e+6)
		{
		the.data <- the.data[seq(1,nrow(the.data),100), ]
		}
	}
} # end for file loop

mytheme <- theme_classic() + 
	   theme(axis.text = element_text(size=14), 
	   axis.title=element_text(size=14), 
	   legend.title = element_text(size=14), 
	   legend.text = element_text(size=14))  

# now use ggplot2 to plot stuff
# If no previous exposure to ggplot2, find one's own way of plotting
# or read the manual at https://ggplot2.tidyverse.org/ 

# first plot: mean_preference
p1 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = Ns, colour="Susceptible")) +
	    geom_line(aes(y = Ni, colour = "Infected")) + 
	    mytheme + 
	    xlab("Generation") + 
	    #ylim(c(0,10)) +
	    ylab("Pop size")

p2 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = mean_freq_p2_total, colour="Overall frequency p2")) +
	    geom_line(aes(y = mean_freq_t2_total, colour="Overall frequency t2")) +
	    mytheme + 
	    xlab("Generation") + 
	    ylab("Overall allele frequencies") +
	    ylim(0,1)
	    
p3 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = mean_freq_p2_susceptible, colour="Freq p2 Susceptible")) +
	    geom_line(aes(y = mean_freq_t2_susceptible, colour="Freq t2 susceptible")) +
	    xlab("Generation") + 
	    ylab("Allele frequencies in Susceptibles") +
	    ylim(0,1) +
	    mytheme  

p4 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = mean_freq_p2_infected, colour="Freq p2 infected")) +
	    geom_line(aes(y = mean_freq_t2_infected, colour="Freq t2 infected")) +
	    xlab("Generation") + 
	    ylab("Allele frequencies in Infected") +
	    ylim(0,1) +
	    mytheme  

p5 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = mean_freq_p2_plasmid, colour="Freq p2 plasmid")) +
	    geom_line(aes(y = mean_freq_t2_plasmid, colour="Freq t2 plasmid")) +
	    xlab("Generation") + 
	    ylab("Allele frequencies in plasmid") +
	    ylim(0,1) +
	    mytheme  

p6 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = mean_freq_p2_chr, colour="Freq p2 chromosome")) +
	    geom_line(aes(y = mean_freq_t2_chr, colour="Freq t2 chromosome")) +
	    xlab("Generation") + 
	    ylab("Allele frequencies in chromosome") +
	    ylim(0,1)+ 
	    mytheme  
p7 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = St1p1, group=rep, colour="St1p1")) +
	    geom_line(aes(y = St2p1, group=rep, colour="St2p1")) +
	    geom_line(aes(y = St1p2, group=rep, colour="St1p2")) +
	    geom_line(aes(y = St2p2, group=rep, colour="St2p2")) +
	    xlab("Generation") + 
	    ylab("Abs frequencies Susceptible genotypes") +
	    mytheme  

p8 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = It1p1t1p1, group=rep, colour="It1p1t1p1")) +
	    geom_line(aes(y = It2p1t1p1, colour="It2p1t1p1")) +
	    geom_line(aes(y = It1p2t1p1, group=rep,colour="It1p2t1p1")) +
	    geom_line(aes(y = It2p2t1p1, group=rep,colour="It2p2t1p1")) +
	    xlab("Generation") + 
	    ylab("Abs frequencies Susceptible genotypes") +
	    mytheme  

p9 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = It1p1t2p1, colour="It1p1t2p1")) +
	    geom_line(aes(y = It2p1t2p1, colour="It2p1t2p1")) +
	    geom_line(aes(y = It1p2t2p1, group=rep,colour="It1p2t2p1")) +
	    geom_line(aes(y = It2p2t2p1, colour="It2p2t2p1")) +
	    xlab("Generation") + 
	    ylab("Abs frequencies Susceptible genotypes") +
	    mytheme  

p10 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = It1p1t1p2, colour="It1p1t1p2")) +
	    geom_line(aes(y = It2p1t1p2, colour="It2p1t1p2")) +
	    geom_line(aes(y = It1p2t1p2, colour="It1p2t1p2")) +
	    geom_line(aes(y = It2p2t1p2, colour="It2p2t1p2")) +
	    xlab("Generation") + 
	    ylab("Abs frequencies Susceptible genotypes") +
	    mytheme  

p11 <- ggplot(data=the.data
	,aes(x=time)) +
	    geom_line(aes(y = It1p1t2p2, colour="It1p1t2p2")) +
	    geom_line(aes(y = It2p1t2p2, colour="It2p1t2p2")) +
	    geom_line(aes(y = It1p2t2p2, colour="It1p2t2p2")) +
	    geom_line(aes(y = It2p2t2p2, colour="It2p2t2p2")) +
	    xlab("Generation") + 
	    ylab("Abs frequencies Susceptible genotypes") +
	    mytheme  

big_plot <- arrangeGrob(p1, 
						p2, 
						p3,
						p4, p5,	
						p6,
						p7, p8,
						p9, p10, p11,
						nrow=3,ncol=4)


output_file_name <- paste(
	"graph_"
	,the.base.name
	,".pdf"
	,sep="")

ggsave(output_file_name, big_plot, width=25,height = 15)

