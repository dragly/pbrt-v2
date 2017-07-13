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
import sprite_bounds as s_bounds
 

BASE_RESOLUTION = 256 # base resolution of sensor (default)
INFINITY = 1e32


################################################################################
# @generate_configuration
################################################################################
def generate_configuration(psh_file, input_directory, output_directory,
        sensor_configuration_template, simulation_method, sensor_data, 
        time_step):
    """
    Creates a string that represents the sensor configuration for pbrt.
    * NOTE: the sensor configuration assumes that the sprite is centered at the
    origin in the integrator code itself, otherwise, it will generate errors.
    * the integrator        renderes a simple image that reflects the bounding
    box of the volume and also a surface the represent the sensor.
    * this configuration also creates a camera that is located at the bottom of
    the volume. this camera will be used to render an image that reflect the
    projection of the volume towards the sensor.

    :param psh_file: the point sprite header file that contains all the
           meta-data of the sprite
    :param input_directory: the directory where the point sprite positions
           are stored
    :param output_directory: where the output images will be generated
    :param sensor_configuration_template: input template pbrt file that
           defines a basic scene configuration for pbrt    
    :param simulation_method: direct, linear or montecarlo
    :param sensor_data: the data of the sensor 
    :param time_step: the time step of the given sprite
    """
    
    sensor_x_resolution = 0
    sensor_y_resolution = 0
    sensor_z_resolution = 0
    
    sensor_width = 0
    sensor_height = 0
    sensor_depth = 0
    
    output_volume = ""
    output_config = ""
    output_vsd_image = ""
    output_pbrt_image = ""
    
    # if the sensor is aligned to the bounding box center, compute the sprite 
    # dimensions and move the sprite to be centered at the origin, otherwise 
    # movie the entire mosaic to be centered at the origin. in this case, the 
    # sprite itself won't be centered at the origin, but the mosaic will be.  
    #if (mosaic_mode):
    #    print ("mosaic mode: set the sensor to cover the meso-circuit")
    # use the data obtained from the configuration file directly
    sensor_x_resolution = sensor_data.x_resolution
    sensor_z_resolution = sensor_data.z_resolution
    
    sensor_width = sensor_data.width
    sensor_height = sensor_data.y_center
    sensor_depth = sensor_data.depth
    
    base_resolution = BASE_RESOLUTION
    sensor_largest_dimension = g_utils.compute_largest(sensor_width, sensor_depth, 0)
    sensor_x_resolution = int(base_resolution * sensor_width / sensor_largest_dimension)
    sensor_z_resolution = int(base_resolution * sensor_depth / sensor_largest_dimension) 

    # read the file data into a single string and update this string step
    # by step with the specified parameters from the simulation.
    config = ""
    with open(sensor_configuration_template, 'r') as pbrt_input_string:
        config = pbrt_input_string.read()

    # update the config with the simulation parameters.
    # camera parameters
    # this is an orthographic camera, the frustum will cover the width and
    # the height of the bounding volume of the vsd data.
    config = config.replace("CAMERA_Y_POSITION", str(sensor_data.y_center))
    config = config.replace("CAMERA_MIN_X", str(-sensor_width / 2.0))
    config = config.replace("CAMERA_MAX_X", str(sensor_width / 2.0))
    config = config.replace("CAMERA_MIN_Z", str(-sensor_depth / 2.0))
    config = config.replace("CAMERA_MAX_Z", str(sensor_depth / 2.0))

    # simulation method
    if(simulation_method == "forward-direct-sprite"):
        config = config.replace("SIMULATION_METHOD_SPRITE", "vsdfds")
    elif(simulation_method == "forward-linear-sprite"):
        config = config.replace("SIMULATION_METHOD_SPRITE", "vsdfls")
    elif(simulation_method == "forward-scattering-sprite"):
        config = config.replace("SIMULATION_METHOD_SPRITE", "vsdfss") 
    elif(simulation_method == "backward-direct-grid"): 
        config = config.replace("SIMULATION_METHOD_VOLUME", "vsdbdg") 
    elif(simulation_method == "backward-linear-grid"): 
        config = config.replace("SIMULATION_METHOD_VOLUME", "vsdblg")
    elif(simulation_method == "backward-scattering-grid"): 
        config = config.replace("SIMULATION_METHOD_VOLUME", "vsdbsg")  
    else:
        print("ERROR: unknown simulation method !")
        exit(0) 
        
    # camera-rendered image name
    camera_image_name = "%s/%s" % (output_directory, time_step)
    config = config.replace("OUTPUT_FILE_NAME", camera_image_name)

    # the resolution of the camera film that will capture the projection of the
    # volume towards the sensor
    config = config.replace("X_RESOLUTION", str(sensor_x_resolution))
    config = config.replace("Y_RESOLUTION", str(sensor_y_resolution))
    config = config.replace("Z_RESOLUTION", str(sensor_z_resolution))

    # input psh file
    config = config.replace("DATA_DIRECTORY", input_directory)
    config = config.replace("PSH_FILE", psh_file)

    # sensor attributes, rectangle sensors
    sensor_image_name = "%s/%s" % (output_directory, time_step)
    config = config.replace("SENSOR_IMAGE_NAME", sensor_image_name)
    config = config.replace("SENSOR_HEIGHT", str(sensor_data.y_center))
    config = config.replace("SENSOR_X", str(sensor_width))
    config = config.replace("SENSOR_Y", str(sensor_depth))
    config = config.replace("X_SENSOR_RESOLUTION", str(sensor_x_resolution))
    config = config.replace("Y_SENSOR_RESOLUTION", str(sensor_z_resolution))

    config = config.replace("X_SHIFT", str(sensor_data.x_center))
    config = config.replace("Y_SHIFT", str(0))
    config = config.replace("Z_SHIFT", str(sensor_data.z_center))
    
    config = config.replace("VOLUME_MIN_X", str(INFINITY))
    config = config.replace("VOLUME_MIN_Y", str(0))
    config = config.replace("VOLUME_MIN_Z", str(INFINITY))
    config = config.replace("VOLUME_MAX_X", str(INFINITY))
    config = config.replace("VOLUME_MAX_Y", str(sensor_data.y_center))
    config = config.replace("VOLUME_MAX_Z", str(INFINITY))
    
    # volume-specific configuration 
    volume_prefix = "%s/%s_volume" % (output_directory, time_step)
    config = config.replace("STEP_SIZE", str(20))

    config = config.replace("VOLUME_PREFIX", volume_prefix)
    
    return config

################################################################################
# @voxelize_sprite
################################################################################
def voxelize_sprite(psh_file,
                    input_directory,
                    output_directory,
                    simulation_method,
                    volumizer_executable,
                    grid_base_resolution,
                    sensor_data,
                    sprite_x_shift,
                    sprite_y_shift,
                    sprite_z_shift):
    """
    Voxelizes the point sprite and returns valid pbrt configuration.

    :param psh_file:
    :param input_directory:
    :param output_directory:
    :param simulation_method:
    :param volumizer_executable:
    :param grid_base_resolution:
    :param sensor_data:
    :param sprite_x_shift:
    :param sprite_y_shift:
    :param sprite_z_shift:
    :return:
    """

    psh_file_prefix = re.sub('.psh', '', psh_file)
    shell_command = " %s %s %s %s %s %s %s %s %s %s %s %s %s %s " % \
        (volumizer_executable,          # 0
         psh_file_prefix,               # 1
         input_directory,               # 2
         output_directory,              # 3
         simulation_method,             # 4
         str(grid_base_resolution),     # 5
         str(sensor_data.x_center),     # 6
         str(sensor_data.y_center),     # 7 
         str(sensor_data.z_center),     # 8
         str(sensor_data.width),        # 9
         str(sensor_data.depth),        # 10
         str(sprite_x_shift),           # 11
         str(sprite_y_shift),           # 12
         str(sprite_z_shift))           # 13
         
    print("Voxelizing the sprite")
    subprocess.call(shell_command, shell=True)


################################################################################
# @get_voxelization_command
################################################################################
def get_voxelization_command(psh_file,
                             input_directory,
                             output_directory,
                             simulation_method,
                             volumizer_executable,
                             grid_base_resolution,
                             sensor_data,
                             sprite_x_shift,
                             sprite_y_shift,
                             sprite_z_shift):
    """
    Gets a voxelization command.

    :param psh_file:
    :param input_directory:
    :param output_directory:
    :param simulation_method:
    :param volumizer_executable:
    :param grid_base_resolution:
    :param sensor_data:
    :param sprite_x_shift:
    :param sprite_y_shift:
    :param sprite_z_shift:
    :return:
    """

    psh_file_prefix = re.sub('.psh', '', psh_file)
    shell_command = " %s %s %s %s %s %s %s %s %s %s %s %s %s %s " % \
                    (volumizer_executable,  # 0
                     psh_file_prefix,  # 1
                     input_directory,  # 2
                     output_directory,  # 3
                     simulation_method,  # 4
                     str(grid_base_resolution),  # 5
                     str(sensor_data.x_center),  # 6
                     str(sensor_data.y_center),  # 7
                     str(sensor_data.z_center),  # 8
                     str(sensor_data.width),  # 9
                     str(sensor_data.depth),  # 10
                     str(sprite_x_shift),  # 11
                     str(sprite_y_shift),  # 12
                     str(sprite_z_shift))  # 13

    print("Sprite Voxelizer")
    return shell_command
