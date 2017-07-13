#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import math

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" %current_path)

import geometric_utilities as g_utils

class MosaicData:
    """
    * A class to encapsulate the data of the mosaic that contains all the columns
    to be simulated.
    """
    def __init__(self):
        # the loaded configuration of the mosaic 
        self.original_x_center = 0.0
        self.original_y_center = 0.0
        self.original_z_center = 0.0
        self.original_x_min = 0.0
        self.original_y_min = 0.0
        self.original_z_min = 0.0
        self.original_x_max = 0.0
        self.original_y_max = 0.0
        self.original_z_max = 0.0
        # the updated configuration after shifting the mosaic at the origin 
        self.current_x_center = 0.0
        self.current_y_center = 0.0
        self.current_z_center = 0.0
        self.current_x_min = 0.0
        self.current_y_min = 0.0
        self.current_z_min = 0.0
        self.current_x_max = 0.0
        self.current_y_max = 0.0
        self.current_z_max = 0.0
        # dimensions 
        self.width = 0.0
        self.height = 0.0
        self.depth = 0.0
        # sampling resolutions (in case of volumizing the sprite)
        self.x_resolution = 0
        self.y_resolution = 0
        self.z_resolution = 0
        return
