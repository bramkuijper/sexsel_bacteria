#!/usr/bin/env Rscript

# (c) Ana Duarte 2021
# R script to generate batch files with which
# many instances of the simulation can be run

# number of replicates
# for each unique parameter combinations
nrep <- 10

# maximum number of times
max_time <- 10e8

# initial fraction of infected individuals
init_fraction_infected <-0.1

# p good plasmid init
p_good_init <- 0.0 

N <- 7000
# density-dependence parameter
kappa <- 1.0 / N

# max fecundity of the host
bmax <- 50

# initial resistance
init_x <- 0.0

# cost of resistance 
c <- seq(0,5,by=0.5) 

# clearance rate good plasmid
gamma_G <- 1
gamma_B <- 0

# cost of trait
#epsilon <- c(0.001,0.01, 0.1)
#epsilon <- 0.0

# cost of plasmid
#delta <- c(0.01,0.1)
#delta <- 0.01

# force of infection 
#psi <- c(0.01,0.05,0.1) 
psi_G <- 0.0 
psi_B <- 10 

# death rates
d <- 1 
dG <- 1 
dB <- 50 

# recombination rate
#r <- 1e-6

# probability of co-infection 
sigma <- 0.0

# mutation rate
mu_x <- 0.01
sdmu_x <- 0.01

# initial frequencies of preference and trait
#init_p2 <- c(0.01,0.1,0.8)
#init_t2 <- c(0.01,0.1,0.8)
#init_p2 <- 0.1
#init_t2 <- 0.1

# plasmid formation rate
#lambda <- 0.0

# preference factor
#alpha <- c(0.001,0.01,0.1)
#alpha <- 2

# preference dominance coefficient
#h <- c(0,0.5,1)
#h <- 0.5

# trait dominance coefficient
#l <- c(0,0.5,1)
#l <- 0.5

# get the directory name of this script
# we use this so that we can
script.dir <- getwd()

# executable
exe <- "./sexsel_bacteria.exe"

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
record_file_name <- file.path(script.dir,paste("batch_record_",curr.time.str,".txt",sep=""))

# obtain all combinations by doing expand.grid()
# which calculates all combinations of factors
# and returns a table of it
# the order in which you put stuff into expand.grid()
# will be the order in which it will be put into 
# the batch file (and hence read in by the simulation)
combinations <- as.data.frame(
        expand.grid(
                exe=exe
                ,max_time=max_time
                ,p_good_init=p_good_init
                ,kappa=kappa
                ,bmax=bmax
                ,c=c
                ,gamma_G=gamma_G
                ,gamma_B=gamma_B
                ,psi_G=psi_G
                ,psi_B=psi_B
                ,d=d
                ,dG=dG
                ,dB=dB
		,sigma=sigma
                ,mu_x=mu_x
                ,sdmu_x=sdmu_x
                ,init_x=init_x
		,init_fraction_infected=init_fraction_infected
		,N=N
                ,stringsAsFactors=F
                ))

nrows <- nrow(combinations)

print(combinations)
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
        output_summary <- paste(output_summary
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
                ,max_time=max_time
                ,p_good_init=p_good_init
                ,kappa=kappa
                ,bmax=bmax
                ,c=c
                ,gamma_G=gamma_G
                ,gamma_B=gamma_B
                ,psi_G=psi_G
                ,psi_B=psi_B
                ,d=d
                ,dG=dG
                ,dB=dB
		,sigma=sigma
                ,mu_x=mu_x
                ,sdmu_x=sdmu_x
                ,init_x=init_x
		,init_fraction_infected=init_fraction_infected
		,N=N
                )
#
#
#
