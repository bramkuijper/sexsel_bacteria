#!/usr/bin/env Rscript

# (c) Ana Duarte 2021
# R script to generate batch files with which
# many instances of the simulation can be run

# number of replicates
# for each unique parameter combinations
nrep <- 1 

# maximum number of times
max_time <- 1e9

# p noplasmid init
p_noplasmid_init <- 0.9 

N <- 7000

# density-dependence parameter
kappa <- 1.0 / N

# max fecundity of the host
bmax <- 50

# cost of preference 
c <- 0.01

# cost of trait
epsilon <- 0.0

# cost of plasmid
delta <- 0.0

# co-infection probablity
sigma <- 0.0

# resistance
pi <- 0.5 

# plasmid acceptance probability (once conjugation happened) 
lambda <- 1.0 

# death rates
d <- 1 
d_t2 <- 2 

# recombination rate
r <- 0.0 

# mutation rate
#mu_p <- c(0.001,0.01)
#mu_t <- c(0.001,0.01)

mu_p_chr <- 1e-5 
mu_t_chr <- 0.0 
mu_p_plasmid <- 0.0 
mu_t_plasmid <- 1e-5 
nu <- 1e-5 

# initial frequencies of preference and trait
init_p2_chr <- 0.1
init_t2_chr <- 0.0
init_p2_plasmid <- 0.0
init_t2_plasmid <- 0.1

# preference factor
#alpha <- c(0.001,0.01,0.1)
alpha <-0.0 

# preference dominance coefficient
h <- 1.0 

# trait dominance coefficient
l <- 1.0

#force infection t1 and t2 plasmids
# currently not necessary?
psi_t1 <- 0.0
psi_t2 <- 0.0

# fitness coefficient for t1 and t2 plasmids 
f_t1 <- 1.0
f_t2 <- 0.1

#coefficient for force of infection t1 and t2 plasmids
beta_t1 <- 1.0
beta_t2 <- 5.0

# clearance rate plasmid
gamma_loss_t1 <- 1
gamma_loss_t2 <- 0

# get the directory name of this script
# we use this so that we can
script.dir <- getwd()

# executable
exe <- "./fishersexsel_Gandon.exe"

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
                ,p_noplasmid_init=p_noplasmid_init
                ,kappa=kappa
                ,bmax=bmax
                ,c=c
		,epsilon=epsilon
		,delta=delta
		,sigma=sigma
                ,pi=pi
		,lambda=lambda
                ,d=d
                ,d_t2=d_t2
		,r=r
                ,mu_p_chr=mu_p_chr
                ,mu_t_chr=mu_t_chr
                ,mu_p_plasmid=mu_p_plasmid
                ,mu_t_plasmid=mu_t_plasmid
                ,nu=nu
                ,init_p2_chr=init_p2_chr
                ,init_p2_plasmid=init_p2_plasmid
                ,init_t2_chr=init_t2_chr
                ,init_t2_plasmid=init_t2_plasmid
		,alpha=alpha
		,h=h
		,l=l
                ,N=N
		,psi_t1=psi_t1
		,psi_t2=psi_t2
		,f_t1=f_t1
		,f_t2=f_t2
		,beta_t1=beta_t1
		,beta_t2=beta_t2
		,gamma_loss_t1=gamma_loss_t1
		,gamma_loss_t2=gamma_loss_t2
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
                ,p_noplasmid_init=p_noplasmid_init
                ,kappa=kappa
                ,bmax=bmax
                ,c=c
		,epsilon=epsilon
		,delta=delta
		,sigma=sigma
                ,pi=pi
		,lambda=lambda
                ,d=d
                ,d_t2=d_t2
		,r=r
                ,mu_p_chr=mu_p_chr
                ,mu_t_chr=mu_t_chr
                ,mu_p_plasmid=mu_p_plasmid
                ,mu_t_plasmid=mu_t_plasmid
                ,nu=nu
                ,init_p2_chr=init_p2_chr
                ,init_p2_plasmid=init_p2_plasmid
                ,init_t2_chr=init_t2_chr
                ,init_t2_plasmid=init_t2_plasmid
		,alpha=alpha
		,h=h
		,l=l
                ,N=N
		,psi_t1=psi_t1
		,psi_t2=psi_t2
		,f_t1=f_t1
		,f_t2=f_t2
		,beta_t1=beta_t1
		,beta_t2=beta_t2
		,gamma_loss_t1=gamma_loss_t1
		,gamma_loss_t2=gamma_loss_t2
                )
#
#
#
