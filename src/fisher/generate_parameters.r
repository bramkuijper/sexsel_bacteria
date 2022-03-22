#!/usr/bin/env Rscript

# R script to generate batch files with which
# many instances of the simulation can be run

# number of replicates
# for each unique parameter combinations
nrep <- 1

# density-dependence parameter
kappa <- 0.0001

# max fecundity of the host
bmax <- 50

# cost of preference 
cp <- c(0.00,0.01)

# clearance rate good plasmid
loss_gamma <- 0.01

# cost of trait
ct <- c(0.01)

# cost of plasmid
delta <- 0.8

# 1- conjugation rate 
pival <- c(0.1)

# death rates
d <- 0.01 

gdb_on <- F

# recombination rate
r <- 1e-5 


mu_p1 <- 0
mu_p2 <- 0

mu_t1 <- 0
mu_t2 <- c(0,0.01)

# initial frequencies of preference and trait
# vary them along a square 
p2_t0 <- c(0.02,seq(0.1,0.9,0.1),0.98)
t2_t0 <- p2_t0

#
#
#init_vals_tp <- rbind(init_vals_tp, cbind(t0_val,rep(0.98,times=length(t0_val))))
#
#init_vals_tp <- rbind(init_vals_tp,cbind(rep(0.02,times=length(p0_val)),p0_val))
#
#init_vals_tp <- rbind(init_vals_tp,cbind(rep(0.98,times=length(p0_val)),p0_val))

#t2_t0 <- init_vals_tp[,1]
#p2_t0 <- init_vals_tp[,2]

D_t0 <- 0.15

max_time <- 5e08

N <- 10000

frac_infected <- 0.2

# preference factor
a <- 1

# whether t is crispr like or environmental
t_crispr <- 1



# preference dominance coefficient
#h <- c(0,0.5,1)
hp <- 0.5

# trait dominance coefficient
#l <- c(0,0.5,1)
ht <- 0.5

# get the directory name of this script
# we use this so that we can
script.dir <- getwd()

# executable
exe <- "./solve_fisher.exe"

if (gdb_on)
{
    exe <- paste("gdb --args",exe)
}

# save current time point to make
# filenames that contain a time stamp
curr.time <- Sys.time()
# put the time in a "dd_mm_YYYY_HHMMSS" format
curr.time.str <- format(x=curr.time, format="%d_%m_%Y_%H%M%S")

# basename for the simulation output file
# i.e., each simulation will generate an output file
# which starts with 'sim_bact_sexsel_06_07_2020_120359'
basename <- paste("sim_bact_sexsel_",curr.time.str,sep="")

# the name of the batch file that contains the commmands to 
# run the simulations
batch_file_name <- file.path(script.dir,paste("batch_file_",curr.time.str,".sh",sep=""))

print(paste("printing batch file ",batch_file_name))

# a file that summarizes the parameters that were varied
# just for the sake of bookkeeping
record_file_name <- file.path(script.dir,
        paste("batch_record_",curr.time.str,".txt",sep=""))

# obtain all combinations by doing expand.grid()
# which calculates all combinations of factors
# and returns a table of it
# the order in which you put stuff into expand.grid()
# will be the order in which it will be put into 
# the batch file (and hence read in by the simulation)
combinations <- as.data.frame(
        expand.grid(
                 exe=exe
                        ,kappa=kappa
                        ,gamma=loss_gamma
                        ,delta=delta
                        ,max_time=max_time
                        ,d=d
                        ,N=N
                        ,frac_infected=frac_infected
                        ,p2_t0=p2_t0
                        ,t2_t0=t2_t0
                        ,D_t0=D_t0
                        ,cp=cp
                        ,ct=ct
                        ,pival=pival
                        ,r=r
                        ,a=a
                        ,hp=hp
                        ,ht=ht
                        ,mu_t1=mu_t1
                        ,mu_t2=mu_t2
                        ,mu_p1=mu_p1
                        ,mu_p2=mu_p2
                        ,t_crispr=t_crispr
                        ,stringsAsFactors=F # prevent strings being recoded into numbers
                ))


nrows <- nrow(combinations)

# index number for output file of
# each simulation
file_idx <- 1


batch_file_contents <- ""

# loop through the replicates 
for (rep_idx in 1:nrep)
{
    # loop through the rows of the combinations
    # dataframe
    for (row_idx in 1:nrows)
    {
        # generate the output file name for the
        # current simulation
        filename <- paste(basename,"_",file_idx,sep="")

        # update count of the file index counter
        file_idx <- file_idx + 1

        # put spaces in between all the elements of the rows
        batch_file_contents <- paste(batch_file_contents,paste(
                combinations[row_idx,]
                ,collapse=" ")
                ,filename
                ,"\n")
    }
}

# write the batch file
write(x=batch_file_contents
        ,file=batch_file_name)

# write a file summarizing the parameters that are varied
summarize.params <- function(...)
{
    # get all the arguments that are supplied to this function
    n.args <- ...length()

    # get all the names of the things that were
    # provided as arguments
    # apparently the first item is always ""
    # so we will skip that
    arg.names <- names(as.list(match.call()))
    arg.names <- arg.names[2:length(arg.names)]

    stopifnot(n.args == length(arg.names))

    output_summary <- ""

    for (i in 1:n.args)
    {
        print(...elt(i))
        output_summary <- paste(
                output_summary
                ,arg.names[[i]]
                ,"="
                ,paste(...elt(i),collapse=" ")
                ,"\n"
                ,sep=""
                )
    }

    write(x=output_summary
            ,file=record_file_name)
}

summarize.params(
                exe=exe
                ,kappa=kappa
                ,gamma=loss_gamma
                ,delta=delta
                ,max_time=max_time
                ,d=d
                ,N=N
                ,p2_t0=p2_t0
                ,t2_t0=t2_t0
                ,D_t0=D_t0
                ,cp=cp
                ,ct=ct
                ,pival=pival
                ,a=a
                ,hp=hp
                ,ht=ht
                ,mu_t1=mu_t1
                ,mu_t2=mu_t2
                ,mu_p1=mu_p1
                ,mu_p2=mu_p2
                )
#
#
#

