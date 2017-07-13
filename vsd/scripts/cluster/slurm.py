#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import subprocess
import glob


import sys
import os
import time
import math

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../utilities" % current_path)

import file_handling as fh
import slurm_configuration as s_config

# Globals
sl = "\n"      # Single new line
dl = "\n\n"    # Double new line
tl = "\n\n\n"  # Triple new line


###############################################################################
def create_batch_config(slurm_config):
    """
    * creates the string that will be set into the batch job file.

    keyword arguments
    :param slurm_config : Slurm configuration parameters.
    :rtype              : String representing the batch job config. header.
    """

    # magic number
    b = "#!/bin/bash%s" % sl

    #########################
    # auto-generated header #
    #########################
    b += "######################################################%s" % sl
    b += "# WARNING - AUTO GENERATED FILE%s" % sl
    b += "# Please don't modify that file manually%s" % sl
    b += "######################################################%s" % sl

    ######################
    # node configuration #
    ######################
    # job name
    b += "#SBATCH --job-name=\"%s%d\"%s" % (slurm_config.job_name,
                                            slurm_config.job_number, sl)

    # number of nodes required to execute the job
    b += "#SBATCH --nodes=%s%s" % (slurm_config.num_nodes, sl)

    # number of cpus per tasks
    b += "#SBATCH --cpus-per-task=%s%s" % (slurm_config.num_cpus_per_task, sl)

    # number of tasks
    b += "#SBATCH --ntasks=%s%s" % (slurm_config.num_tasks_per_node, sl)

    # memory required per task in Mbytes
    b += "#SBATCH --mem=%s%s" % (slurm_config.memory_mb, sl)

    # slurm session time
    b += "#SBATCH --time=%s%s" % (slurm_config.session_time, sl)

    # job partition
    b += "#SBATCH --partition=%s%s" % (slurm_config.partition, sl)

    # job account
    b += "#SBATCH --account=%s%s" % (slurm_config.project_name, sl)

    # On which nodes, this job will be executed
    # This option is used if the required modules are installed on a specific
    # node
    # b += "#SBATCH --nodelist=%s%s" % (slurm_config.node_list, sl)

    #####################
    # user notification #
    #####################
    if slurm_config.enable_email_notification:
        b += "#SBATCH --mail-type=ALL%s" % sl
        b += "#SBATCH --mail-user=%s%s" % (slurm_config.user_email, sl)

    ##################
    # log generation #
    ##################
    if slurm_config.enable_logs:
        std_out = "%s/slurm-stdout_%d.log" % \
                  (slurm_config.log_files_path, slurm_config.job_number)
        std_err = "%s/slurm-stderr_%d.log" % \
                  (slurm_config.log_files_path, slurm_config.job_number)
        b += "#SBATCH --output=%s%s" % (std_out, sl)
        b += "#SBATCH --error=%s%s" % (std_err, dl)

    ####################
    # System variables #
    ####################
    # slurm profile
    b += "# Loading profiles%s" % sl
    b += "%s%s" % (slurm_config.profile, dl)

    # job home
    b += "#JOB_HOME=\"%s\"%s" % (slurm_config.execution_path, sl)

    # KERBEROS renewal
    b += "# Renewal of KERBEROS periodically for the length of the job%s" % sl
    b += "krenew -b -K 30%s" % dl

    # slurm modules
    b += "# Loading the modules.%s" % sl
    b += "%s%s" % (slurm_config.modules, dl)

    # environmental variables
    b += "# Setting the environmental variables.%s" % sl
    b += "export PATH=%s:$PATH%s" % (slurm_config.env_path, sl)
    b += "export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH%s" % \
         (slurm_config.env_ld_library_path, sl)
    b += "export PYTHONPATH=%s:$PYTHONPATH%s" % (slurm_config.env_python_path,
                                                 dl)
    # node list
    b += "echo \"On which node your job has been scheduled :\"%s" % sl
    b += "echo $SLURM_JOB_NODELIST%s" % dl

    # shell limits
    b += "echo \"Print current shell limits :\"%s" % sl
    b += "ulimit -a%s" % dl

    # running the serial tasks.
    b += "echo \"Now run your serial tasks ...\"%s" % sl
    b += "cd %s%s" % (slurm_config.execution_path, dl)
    ####################################################################

    return b


###############################################################################
def submit_batch_scripts(scripts_directory):
    """
    * lists all the generated slurm scripts and execute them one by one via
    the <sbatch> command.

    keyword arguments
    :param scripts_directory : where the batch scripts are generated.
    """

    # list all the scripts in the script directory.
    script_list = fh.list_files_in_directory(scripts_directory, "sh")

    if len(script_list) == 0:
        print("There are no scripts to execute")

    # submit each script to the batch partition in the cluster.
    for script in script_list:

        # change the permissions
        shell_command = "chmod +x %s/%s" % (scripts_directory, script)
        subprocess.call(shell_command, shell=True)

        # execute the script
        shell_command = "sbatch %s/%s" % (scripts_directory, script)
        print(shell_command)

        # execute the shell command.
        subprocess.call(shell_command, shell=True)
