#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import numpy
from numpy import *

import os
import sys

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../geometry" % current_path)
sys.path.append("%s/../sprite" % current_path)


def create_volume_from_sprite(sprite_positions, sprite_intensities,
                              sprite_bounding_box, largest_side_size):
    """
    * creates a sampled volume from the point sprite
    """
    # find the volume dimensions from the sprite bounding box information
    volume_width, volume_height, volume_depth = \
        sprite_bounding_box.get_volumetric_grid_dimensions(largest_side_size)

    print ("grid dimensions %d X %d X %d " % (volume_width,
                                              volume_height,
                                              volume_depth))

    # create the volume
    volume = numpy.zeros((volume_width, volume_height, volume_depth))

    # append the relative power (from 0 to 1) of every event at its
    # corresponding voxel
    for i_position, i_intensity in zip(sprite_positions, sprite_intensities):

        # map the coordinates of the event to the index of the voxel
        x_volume, y_volume, z_volume = \
            sprite_bounding_box.convert_coordinate_to_volume_grid_index \
            (i_position[0], i_position[1], i_position[2], largest_side_size)

        # add the power of the sprite to the voxel
        volume[x_volume][y_volume][z_volume] = \
            volume[x_volume][y_volume][z_volume] + i_intensity

    return volume


