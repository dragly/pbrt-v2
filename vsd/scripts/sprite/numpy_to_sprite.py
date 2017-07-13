#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import struct
import random
import os
import sys
import time
import ntpath
import numpy 
from numpy import *


################################################################################
def load_intensity_array(file_path):

    # retrieve the time step from the 'shitty' file name
    time_step = file_path.split("sprites")[1]
    time_step = time_step.split("sprite")[1]
    intensity_array = numpy.load(file_path)
    return intensity_array, time_step


################################################################################
def load_position_array(file_path):

    # get the position array from the file 
    position_array = numpy.load(file_path)
    return position_array
