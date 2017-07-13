#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015 - 2016
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import math

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" %current_path)

class SensorData:
    """
    * A class to encapsulate the data of the sensor that will record the 
    activity of cortex. The sensor is always assumed to be centered at the 
    origin in the xz plane in pbrt.
    """
    def __init__(self):
        self.width = 0.0
        self.depth = 0.0
        self.x_center = 0.0
        self.y_center = 0.0
        self.z_center = 0.0
        self.x_resolution = 0
        self.z_resolution = 0
        return
        
    def print_data(self):
        print("Sensor: [w=%f, h=%f, x=%f, y=%f, z=%f, res=%dx%d]") % \
            (self.width, self.depth, 
             self.x_center, self.y_center, self.z_center, 
             self.x_resolution, int(self.z_resolution))
