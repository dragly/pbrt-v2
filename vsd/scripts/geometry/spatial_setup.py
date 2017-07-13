#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2016
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import sys
import os
import time
import re
import math
import random
import subprocess
import ConfigParser 

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

def configure_spatial_setup(data_config_file):
    """
    * create a spatial setup for the sprite with respect to the sensor that 
    should cover the entire mosaic (a 5x5 meso-scale circuit)
    
    * keyword arguments
    :param data_config_file: an input data config file that contains all the 
           parameters of the sensor and each column in the meso-scale circuit. 
    
    """
    ### mosaic data 
    # open the circuit configuration file and get how many columns 
    config = ConfigParser.ConfigParser()
    config.readfp(open(data_config_file))
    columns = config.options(config.sections()[0])
    data = config.options(config.sections()[1])
    
    # left lower corner of the mosaic
    lower_left_corner = config.get(str(config.sections()[1]),
                                   "lower left corner").split(" ")
    lower_left_corner_x = float(lower_left_corner[0])
    lower_left_corner_y = float(lower_left_corner[1])
    
    # upper right corner of the mosaic
    upper_right_corner = config.get(config.sections()[1], 
                                    "upper right corner").split(" ")
    upper_right_corner_x = float(upper_right_corner[0])
    upper_right_corner_y = float(upper_right_corner[1])
    
    # mosaic dimensions 
    # @note: width and height are changed here 
    mosaic_width = config.getfloat(config.sections()[1], "mosaic width")    # x
    mosaic_height = config.getfloat(config.sections()[1], "mosaic height")   # y
    mosaic_depth = config.getfloat(config.sections()[1], "mosaic depth")   # z
    
    # mosaic center 
    mosaic_center = config.get(config.sections()[1], "mosaic center").split(" ")
    mosaic_x_center = float(mosaic_center[0])
    mosaic_y_center = mosaic_height / 2.0   # starts from 0 at the bottom 
    mosaic_z_center = float(mosaic_center[1])
    
    ### sensor data
    # sensor dimensions
    # @note: height and depth are interchangable for consistency  
    # @todo: update the resolution from elsewhere ! 
    sensor_resolution = config.getint(config.sections()[2], "sensor resolution")
    sensor_width = config.getfloat(config.sections()[2], "sensor width")
    sensor_depth = config.getfloat(config.sections()[2], "sensor height")
    sensor_position = config.get(config.sections()[2], 
                                                "sensor position").split(" ")
    sensor_position_x = float(sensor_position[0])
    sensor_position_y = float(sensor_position[1])
    sensor_position_z = float(sensor_position[2])
    
    sensor_largest_dimension = 0
    if(sensor_width > sensor_depth): sensor_largest_dimension = sensor_width
    else: sensor_largest_dimension = sensor_depth
    sensor_x_scale = float(sensor_width) / float(sensor_largest_dimension)
    sensor_z_scale = float(sensor_depth) / float(sensor_largest_dimension)
    
    # @note: width and height are changed here
    sensor_data = s_data.SensorData()
    sensor_data.width = sensor_width
    sensor_data.depth = sensor_depth
    sensor_data.x_resolution = int(sensor_resolution * sensor_x_scale)
    sensor_data.z_resolution = int(sensor_resolution * sensor_z_scale)
    sensor_data.x_center = sensor_position_x
    sensor_data.y_center = sensor_position_y
    sensor_data.z_center = sensor_position_z
    
    # @note: the mosaic has been shifted to the origin 
    mosaic_data = m_data.MosaicData()
    mosaic_data.width = mosaic_width
    mosaic_data.height = mosaic_height
    mosaic_data.depth = mosaic_depth
    # the original loaded configuration 
    mosaic_data.original_x_center = mosaic_x_center
    mosaic_data.original_y_center = mosaic_y_center
    mosaic_data.original_z_center = mosaic_z_center
    # the mosaic will be shifted to the origin (0, 0, 0)
    mosaic_data.current_x_center = 0.0
    mosaic_data.current_y_center = 0.0
    mosaic_data.current_z_center = 0.0
    # the original loaded configuration 
    mosaic_data.original_x_min = lower_left_corner_x
    mosaic_data.original_y_min = 0
    mosaic_data.original_z_min = lower_left_corner_y
    mosaic_data.original_x_max = upper_right_corner_x
    mosaic_data.original_y_max = mosaic_height
    mosaic_data.original_z_max = upper_right_corner_y
    # the bounding box after shifting the mosaic to the origin 
    mosaic_data.current_x_min = -mosaic_width / 2.0
    mosaic_data.current_y_min = -mosaic_height / 2.0
    mosaic_data.current_z_min = -mosaic_depth / 2.0
    mosaic_data.current_x_max =  mosaic_width / 2.0
    mosaic_data.current_y_max =  mosaic_height / 2.0
    mosaic_data.current_z_max =  mosaic_depth / 2.0
    return sensor_data, mosaic_data 
