#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015 - 2016
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
# Description : This script is used to run the VSDI simulation framework.
# For help, see the argument parser help or type
# ./run_in_silico_vsdi.py --help in your console.
################################################################################

# system imports
import os
import sys
import subprocess
import glob

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../pbrt" % current_path)
sys.path.append("%s/../utilities" % current_path)
sys.path.append("%s/../geometry" % current_path)
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../cluster" % current_path)

import file_handling as fh
import pbrt_configuration as pbrt
import sensor_data as s_data
import spatial_setup as s_setup
import sprite_bounds as s_bounds
import sprite_reader as s_reader
import arguments
import slurm
import slurm_configuration as s_config


################################################################################
# @create_executable_command
################################################################################
def create_executable_command(args, psh_file):
    """
    Create executables.

    :param args:
    :param psh_file:
    :return:
    """

    # get the time step for the psh file
    time_step = s_reader.get_time_step_from_psh(
        args.input_directory + "/" + psh_file)

    print("Time frame : %s" % time_step)
    results_directory = "%s/results/%s" % (args.output_directory, time_step)
    fh.clean_and_create_new_directory(results_directory)

    vsd_sensor_data = 0
    if (not (args.data_config_file == 'NO_FILE_PROVIDED')):

        # configure the spatial setup for the experiment
        vsd_sensor_data, mosaic_data = \
            s_setup.configure_spatial_setup(args.data_config_file)
        vsd_sensor_data.print_data()
    else:
        print('ERROR: a configuration file for the circuit MUST be given !')
        exit(0)

    # @note: if the selected method uses sprites, proceed with sprite reading &
    # data generation from the python workflow, else use the volumizesprite
    # tool to sample the space and create volume grid that reflects the
    # aggregation of the events in volumes
    sensor_configuration = ""
    voxelization_command = ""
    simulation_command = ""

    ### the vsd signal is a sprite, forward method
    if (args.simulation_method == "forward-direct-sprite") or \
            (args.simulation_method == "forward-linear-sprite") or \
            (args.simulation_method == "forward-scattering-sprite"):
        # create the configuration for sprite
        sensor_configuration = pbrt.generate_configuration(psh_file,
            args.input_directory, results_directory,
            args.pbrt_sprite_sensor_config, args.simulation_method,
            vsd_sensor_data, time_step)

    ### the vsd signal is a volume, backward method
    elif (args.simulation_method == "backward-direct-grid") or \
         (args.simulation_method == "backward-linear-grid") or \
         (args.simulation_method == "backward-scattering-grid"):

        sprite_x_shift = vsd_sensor_data.x_center
        sprite_y_shift = vsd_sensor_data.y_center
        sprite_z_shift = vsd_sensor_data.z_center

        # voxelize the sprite
        voxelization_command = pbrt.get_voxelization_command(psh_file,
            args.input_directory, results_directory, args.simulation_method,
            args.volumizer_executable, args.grid_resolution, vsd_sensor_data,
            sprite_x_shift, sprite_y_shift, sprite_z_shift)

        # create the configuration for the volume
        sensor_configuration = pbrt.generate_configuration(psh_file,
            args.input_directory, results_directory,
            args.pbrt_volume_sensor_config, args.simulation_method,
            vsd_sensor_data, time_step)

    else:
        print("error, the selected method [%s] is not know!" %
              args.simulation_method)
        exit(0)

    # if a template of the sensor configuration was given, generate the
    # corresponding sensor configuration file based on the current
    # distribution of the point sprites
    if sensor_configuration:

        # configuration
        sensor_configuration_path = "%s/%s.pbrt" % (
            results_directory, time_step)
        fh.save_file_to_disk(sensor_configuration_path,
            sensor_configuration)

        # get the simulation command
        simulation_command = "%s %s" % \
                        (args.pbrt_executable, sensor_configuration_path)
    else:
        print('ERROR: pbrt configuration file is not created !')
        exit(0)

    execution_command = voxelization_command + "\n" + simulation_command
    return execution_command


################################################################################
if __name__ == "__main__":

    print("*****************************************************")
    print("In Silico VSDI (Taylor H. Newton and Marwan Abdellah)")
    print("*****************************************************")

    # parse the arguments
    args = arguments.parse_arguments_for_single_machine()

    # a list that will keep all the paths of the time series
    psh_file_list = []

    # simulate the entire time series if you can and if stated 'series'
    if (args.psh_file == "series"):
        # load all the available time steps
        os.chdir(args.input_directory)
        for psh_file in glob.glob("*.psh"):
            psh_file_list.append(psh_file)
    else:
        psh_file_list.append(args.psh_file)

    # clean the output directory
    fh.clean_and_create_new_directory(args.output_directory)

    # create slurm directory
    slurm_directory = "%s/%s" % (args.output_directory, "slurm")
    fh.clean_and_create_new_directory(slurm_directory)

    logs_directory = "%s/%s" % (args.output_directory, "logs")
    fh.clean_and_create_new_directory(logs_directory)

    # execute the workflow on the cluster
    if args.node == 'cluster':

        # write all the batch jobs to the directory
        for i, psh_file in enumerate(psh_file_list):

            # create the slurm configuration.
            slurm_configuration = s_config.SlurmConfiguration()
            slurm_configuration.job_number = i
            slurm_configuration.log_files_path = logs_directory

            # create the slurm configuration string.
            slurm_batch_string = slurm.create_batch_config(slurm_configuration)

            # create the execution command and then append it to the batch string
            execution_command = create_executable_command(args, psh_file)
            slurm_batch_string += execution_command

            # save the script
            script_path = "%s/%s.sh" % (slurm_directory, str(i))
            fh.save_file_to_disk(script_path, slurm_batch_string)

        # run all the batch scripts
        slurm.submit_batch_scripts(slurm_directory)

    # execute the workflow locally
    elif args.node == 'local':

        # run the scripts directly
        for i, psh_file in enumerate(psh_file_list):

            # create the execution command
            execution_command = create_executable_command(args, psh_file)

            # run the execution command
            #print(execution_command)
            subprocess.call(execution_command, shell=True)
    else:
        print('ERROR:execution node [%s] is not specified or wrong!' %
              args.node)
        exit(0)
