#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015 - 2016
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import sys
import os
import time
import re
import math
import random
import subprocess

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../geometry" %current_path)
sys.path.append("%s/../pbrt" %current_path)
sys.path.append("%s/../sprite" %current_path)
sys.path.append("%s/../utilities" % current_path)
sys.path.append("%s/../volume" %current_path)

import sprite_reader as s_reader
import bounding_box as bb
import geometric_utilities as g_utils
import sensor_data as s_data 
import mosaic_data as m_data

def get_sprite_bounds(input_directory, output_directory, psh_file, 
                      spritebounds_executable):
    """
    * read the sprite in C++ and return the bounding box before and after the 
    translation to the origin. this way faster than implementing it in python !
    
    * keyword arguments
    :param psh_file: name of the sprite header file
    :param input_directory: the directory where the vsd point sprites exist
    :param output_directory: the directory where the pbrt configurations will be
           generated
    """                      
    
    # extract the prefix from the file name 
    psh_prefix = psh_file.replace(".psh", "")

    # execute the spritebounds executable to generate a text file with the 
    # bounding box data of the sprite
    # @note: this process is executed in C++ since it is a lot faster than doing 
    # it in python ! 
    # @note the 000 are the shift vector that is not needed here 
    shell_command = "%s %s %s %s %s %s %s " % \
        (spritebounds_executable, psh_prefix, input_directory, output_directory,
         0, 0, 0)
    subprocess.call(shell_command, shell=True)

    # read the generated file that contain the bounds and get the data  
    original_bounds_path = "%s/%s.original.bounds" % (output_directory, psh_prefix)
    original_bounds_file = open(original_bounds_path,"r")
    original_sprite_bounds = original_bounds_file.readlines()
    original_sprite_bounds = original_sprite_bounds[0] 
    original_sprite_bounds = original_sprite_bounds.split(',')
    
    original_x_min = float(original_sprite_bounds[0])
    original_y_min = float(original_sprite_bounds[1])
    original_z_min = float(original_sprite_bounds[2])
    original_x_max = float(original_sprite_bounds[3])
    original_y_max = float(original_sprite_bounds[4])
    original_z_max = float(original_sprite_bounds[5])
    original_bounds_file.close()
    
    original_bb = bb.BoundingBox(original_x_min, original_y_min, original_z_min, 
                                 original_x_max, original_y_max, original_z_max)
    
    centered_bounds_path = "%s/%s.centered.bounds" % (output_directory, psh_prefix)
    centered_bounds_file = open(centered_bounds_path,"r")
    centered_sprite_bounds = centered_bounds_file.readlines()
    centered_sprite_bounds = centered_sprite_bounds[0]
    centered_sprite_bounds = centered_sprite_bounds.split(',')
    centerd_x_min = float(centered_sprite_bounds[0])
    centerd_y_min = float(centered_sprite_bounds[1])
    centerd_z_min = float(centered_sprite_bounds[2])
    centerd_x_max = float(centered_sprite_bounds[3])
    centerd_y_max = float(centered_sprite_bounds[4])
    centerd_z_max = float(centered_sprite_bounds[5])
    centered_bounds_file.close()
    
    centered_bb = bb.BoundingBox(centerd_x_min, centerd_y_min, centerd_z_min, 
                                 centerd_x_max, centerd_y_max, centerd_z_max)
    
    # return both bounding boxes 
    return original_bb, centered_bb
    
     
