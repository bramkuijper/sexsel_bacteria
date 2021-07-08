#!/usr/bin/env python3

# takes as 
import os
import os.path
import datetime
import time
import re
import json
import sys
import shutil
import string
import argparse

class JobArray:

    # pattern to match lines starting with './exename'
    linepattern = re.compile("^\.\/([a-zA-Z_\d]*.exe)\s+(.*)")

    # list with command line args
    linelist = []

    # list with executables
    executables = []

    email = "a.l.w.kuijper@exeter.ac.uk"

    parent_dir_prefix = "hpcbatch"
    jobfile_prefix = "hpcjob"
    core_folder_prefix = "core"
    runtime_mins = 300 
    maxprocesses_per_batch = 100

    bash_jobarray_varname = "jobarray"
    bash_sge_index = "SLURM_ARRAY_TASK_ID"

    ntasks_per_node = 1
    memory = "1gb"

    qsub_command_filename = "qsub_command.sh"

    # the directory in which we want to find
    # the execuble. By convention this should
    # be the directory in which we are now
    exedir = os.getcwd()

    # the location of the configuration file
    config_file = "jobarray_config.json"

    def __init__(self):

        parser = argparse.ArgumentParser(
                description="Generate batchfiles that " +\
                "can be submitted through slurm")

        # specify output file name
        parser.add_argument('-i', help="The file containing all the jobs that need to be put in batches")
        args = vars(parser.parse_args())

        self.job_input_file = args["i"]

        if not self.job_input_file:
            parser.print_help()
            sys.exit(1)

        if not os.path.exists(self.job_input_file):
            print("cannot find file " + self.job_input_file)
            sys.exit(1)
    
    def get_settings_from_file(self):
        """
        Read the settings file and check its contents

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        #try:
        #    json.load(settings_file)

        # check whether json config file exists
        if not os.path.exists(settings_file):
            raise Exception("Cannot find the json settings file that contains your email address, etc")


    # read in the file with jobs
    def readin(self):

        """ 
        Reads in the file with job commands and putch it into batches

        Parameters
        ----------
        self: the parent object

        filename: str
            the filename of the file which contains the job commands

        Returns
        -------
        Nothing

        """

        # open a batch file (generated with some other
        # script) which contains a list of commands
        f = open(self.job_input_file, "r")

        alllines = f.readlines()

        f.close()

        # loop through all the lines and see whether 
        # the current line is indeed a list of all commands
        for line in alllines:

            # see if this line is one with an executable and
            # command line arguments
            mo = re.search(self.linepattern,line)

            if mo is not None:
                
                self.linelist.append(line.strip())
                self.executables.append(mo.group(1))

        # just retain only the unique names of the exes
        self.executables = set(self.executables)

        if len(self.linelist) > 0:
            self.make_batches()
        else:
            print("No lines collected that matched the begin-pattern. Cowardly exiting.")

    # the main executing function
    def make_batches(self):

        # loop through the list of lines and group them in 
        # lists of lines of size maxprocesses_per_batch
        for starting_line in range(0
                ,len(self.linelist)
                ,self.maxprocesses_per_batch):

            # get the last line of the list of joblines
            # this is either just the last line or 
            # self.maxprocesses_per_batch-th line
            ending_line = starting_line + self.maxprocesses_per_batch
            if ending_line > len(self.linelist):
                endling_line = len(self.linelist)

            # wait so that timestamps of folders
            # are not identical
            time.sleep(3)

            # create main folder
            now = datetime.datetime.now()

            # set the time of creation of the parent directory
            # which will also be used as an id of the jobfiles
            self.idtime = now.strftime("%d_%m_%Y_%H%M%S")

            homedir = os.path.expanduser("~")

            # create the batch directory
            self.batch_dir = os.path.join(
                    homedir
                    ,self.parent_dir_prefix + "_" + self.idtime)

            os.mkdir(self.batch_dir)

            # get current wd
            olddir = os.getcwd()

            print("creating directory " + self.batch_dir)

            # go to that directory
            os.chdir(self.batch_dir)

            linelist_subset = self.linelist[starting_line:ending_line]

            assert(len(linelist_subset) > 0)

            # now make the jobfile
            self.make_jobfile(linelist_subset)

            # go back to previous dir
            os.chdir(olddir)

    # get all the command line commands, i.e.,
    # ./exe 0.1 0.3 0.4 0.5 0.8 0.9
    # and put them into a bash array so that they
    # can be selected in a certain job
    def joblines_to_jobfile_content(self, linelist):

        contents = ""
        for i, line in enumerate(linelist):
            contents += self.bash_jobarray_varname + "[" + str(i + 1) + "]=\"" + line + "\"\n"

        number_jobs = len(linelist)

        return((number_jobs, contents))


    # make the jobfile 
    # corresponding to each directory
    def make_jobfile(self, linelist):

        # the jobname should be the main folder name and the jobnumber
        jobname = os.path.basename(self.batch_dir)

        # copy each executable to the batch folder
        for exename in self.executables:

            exepath_from = os.path.join(self.exedir, exename.strip())
            exepath_to = os.path.join(self.batch_dir, exename.strip())

            # cannot find the executable in the 'from' directory
            if not os.path.exists(exepath_from):
                print("Could not find the executable file " 
                        + exename.strip() + " within the directory " 
                        + self.exedir + ". Either the name of the executable is wrong, or one should cd to the directory in which the executable file resides and run jobarray from there.")


            if not os.path.exists(self.batch_dir):
                print("Woah, something goes wrong. I cannot copy the executable " + exename.strip() + " to the new directory " + self.batch_dir + ", because that directory does not exist, despite the fact that I created it.")

            # copy 
            shutil.copy(
                    exepath_from
                    ,exepath_to
                    )
                 
        # we are not going to specify a walltime anymore
#        walltime = self.exes_per_core * self.runtime_mins

#        parentdirs = os.getcwd().split("/")

        number_tasks, the_content = self.joblines_to_jobfile_content(linelist)

        # make the jobfile
        clustercontent = "#!/bin/bash -l\n"  \
                        + "#SBATCH --job-name=" + jobname + "\n"  \
                        + "#SBATCH --mail-type=NONE\n" \
                        + "#SBATCH --mail-user=" + self.email + "\n" \
                        + "#SBATCH --ntasks=" + str(self.ntasks_per_node) + "\n" \
                        + "#SBATCH --mem=" + self.memory + "\n" \
                        + "#SBATCH --output=" + jobname + "_%A_%a.log\n" \
                        + "#SBATCH --array=1-" + str(number_tasks) + "\n" \
                        + "cd " + self.batch_dir + "\n"  \
                        + the_content + "\n\n" \
                        + "srun -n 1 --ntasks-per-node=" + str(self.ntasks_per_node) + " \"$(${" + self.bash_jobarray_varname + "[$" + self.bash_sge_index + "]})\""

        self.jobfilename = self.jobfile_prefix + "_" + self.idtime

        f = open(self.jobfilename, "w")
        f.write(clustercontent)
        f.close()

    def make_qsub_command_file(self):

        qsub_command = "#!/usr/bin/env bash\n\n"
        qsub_command += "qsub"

        assert(self.jobfilename is not None)

        # specify array min and max (i.e., total number of job ids, starting from 1)
        qsub_command += " -t 1-" + str(len(self.linelist)) + " -tc " + str(self.maxprocesses)
        qsub_command += " " + self.jobfilename + "\n"

        f =  open(self.qsub_command_filename,"w")
        f.write(qsub_command)
        f.close()

splot = JobArray()

splot.readin()

