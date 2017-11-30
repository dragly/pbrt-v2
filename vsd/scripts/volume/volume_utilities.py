#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################


import numpy
from numpy import *
import struct
import os
import sys


def map_volume_to_positive_values(volume):
    """
    * if the volume was negative, change the values to positive for
    visualizing the distribution.

    """
    volume_width = volume.shape[0]
    volume_height = volume.shape[1]
    volume_depth = volume.shape[2]

    for i in range(0, volume_width):
        for j in range(0, volume_height):
            for k in range(0, volume_depth):
                if (volume[i][j][k] < 0):
                    volume[i][j][k] = -1 * volume[i][j][k]

    return volume


def map_volume_to_ubyte(volume):
    """
    * maps the volume data to ubyte.

    """
    volume_width = volume.shape[0]
    volume_height = volume.shape[1]
    volume_depth = volume.shape[2]

    # find the maximum value in the volume
    max_value = sys.float_info.min
    for i in range(0, volume_width):
        for j in range(0, volume_height):
            for k in range(0, volume_depth):
                if (volume[i][j][k] > max_value):
                    max_value = volume[i][j][k]

    print ("max value %f" % max_value)

    # map the volume between 0 and 255
    for i in range(0, volume_width):
        for j in range(0, volume_height):
            for k in range(0, volume_depth):
                volume[i][j][k] = int(1.0 * volume[i][j][k] / max_value * 255)

    return volume


################################################################################
def write_volume_to_raw_file(volume, prefix):
    """
    """

    volume_width = volume.shape[0]
    volume_height = volume.shape[1]
    volume_depth = volume.shape[2]

    # write the header file the contains the meta-data
    # (dimensions for the moment)
    header_file = open("%s.hdr" % prefix, 'w')
    header_file.write("%d %d %d" % (volume_width, volume_height, volume_depth))
    header_file.close()

    # write the raw (img) file that contains the data
    raw_file = open("%s.raw" % prefix, 'wb')
    print("Writing volume to", prefix, ".raw")
    for i in range(0, volume_width):
        for j in range(0, volume_height):
            for k in range(0, volume_depth):
                raw_file.write(struct.pack('B', volume[i][j][k]))
    raw_file.close()
