#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../geometry" % current_path)
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" %current_path)

import sprite_generation as sg
import sprite_writer as sw
import bounding_box as bb
import geometric_utilities as g_utils


################################################################################
if __name__ == "__main__":

    # create a single time step sprite
    position_array, intensity_array = sg.generate_random_point_sprite \
    (number_points=1000000,
     x_min=-10, y_min=-3, z_min=-1,
     x_max=10, y_max=20, z_max=4)

    print position_array

    # write the sprite to .psh, .psp, and .psi files
    sprite_prefix = 'nh-random'
    sw.write_point_sprite_header \
        (".", sprite_prefix, sprite_prefix, sprite_prefix,
         position_array, 0)
    sw.write_point_sprite_positions \
        (".", sprite_prefix, position_array)
    sw.write_point_sprite_intensity \
        (".", sprite_prefix, intensity_array)

    # write the sprite to an ASCII file
    sw.write_point_sprite_to_file(file_path=".",
                                  file_name=sprite_prefix,
                                  file_format="ascii",
                                  position_array=position_array,
                                  intensity_array=intensity_array)




