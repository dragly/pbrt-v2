#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2016
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

# system imports
import os
import sys
import time
import sys
import subprocess
import argparse
import re

################################################################################
def parse_arguments_for_single_machine():
    """
    * argument parser for the script that runs on a single machine to simulate 
    the vsd activity in the cortex
    """

    # setup command line parser and parse arguments
    parser = argparse.ArgumentParser()

    arg_help = 'the path of the point sprite data (positions and intensities)' 
    parser.add_argument('--input-directory',
                        action='store',
                        dest='input_directory',
                        help=arg_help)
                        
    arg_help = 'the directory where the configuration files will be generated'
    parser.add_argument('--output-directory',
                        action='store',
                        dest='output_directory',
                        help=arg_help)

    arg_help = 'the point sprite header file (meta data of the sprite)'
    parser.add_argument('--psh-file',
                        action='store', default='NO_FILE_PROVIDED',
                        dest='psh_file',
                        help=arg_help)
                        
    arg_help = 'the data configuration file for the circuit and the sensor'
    parser.add_argument('--data-config-file',
                        action='store', default='NO_FILE_PROVIDED',
                        dest='data_config_file',
                        help=arg_help)

    arg_help = 'simulation method, direct-sprite, linear-sprite, ... '
    parser.add_argument('--simulation-method',
                        action='store', default='direct-sprite',
                        dest='simulation_method',
                        help=arg_help)

    arg_help = 'the template pbrt sensor configuration file for sprite'
    parser.add_argument('--pbrt-sprite-sensor-config',
                        action='store', default='NO_FILE_PROVIDED',
                        dest='pbrt_sprite_sensor_config',
                        help=arg_help)
                        
    arg_help = 'the template pbrt sensor configuration file for volume'
    parser.add_argument('--pbrt-volume-sensor-config',
                        action='store', default='NO_FILE_PROVIDED',
                        dest='pbrt_volume_sensor_config',
                        help=arg_help)

    arg_help = 'the path of the pbrt executable that will run the simulation'
    parser.add_argument('--pbrt-executable',
                        action='store', default='pbrt', # installed
                        dest='pbrt_executable',
                        help=arg_help)

    arg_help = 'the path of the sprite volumeizer executable that will ' \
               'convert the sprite to a volume '
    parser.add_argument('--volumizer-executable',
                        action='store', default='volumizesprite', # installed 
                        dest='volumizer_executable',
                        help=arg_help)
                        
    arg_help = 'the path of the spritebounds executable that will ' \
               'quicky extract the bounds of the sprite to get the sensor data'
    parser.add_argument('--sprite-bounds-executable',
                        action='store', default='spritebounds', # installed 
                        dest='spritebounds_executable',
                        help=arg_help)

    arg_help = 'resolution of the grid converted from the sprite'
    parser.add_argument('--grid-resolution',
                    action='store', default='512',
                    dest='grid_resolution',
                    help=arg_help)

    arg_help = 'running node, cluster or local, cluster by default'
    parser.add_argument('--node',
        action='store', default='cluster',
        dest='node',
        help=arg_help)
                    
    arg_help = 'the base (maximum) resolution of the sensor'
    parser.add_argument('--sensor-resolution',
                        action='store', default='512', 
                        dest='sensor_resolution',
                        help=arg_help)

    # parse the arguments
    args = parser.parse_args()
    
    return args
