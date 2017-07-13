#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import time

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../geometry" % current_path)
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" %current_path)

import sprite_generation as sg
import bounding_box as bb
import geometric_utilities as g_utils
import volume_generation as vg
import volume_utilities as v_utils


################################################################################
def create_volume_from_sprite(volume_largest_side_size):

    # generate a random sprite
    print ("generating a random sprite ")
    position_array, intensity_array = \
        sg.generate_random_point_sprite(number_points=10000)

    # find the bounds of the sprite
    print ("compute bounding box of the sprite")
    x_min, y_min, z_min, x_max, y_max, z_max = \
        g_utils.compute_bounds_of_sprite(position_array)
    sprite_bounding_box = bb.BoundingBox(x_min, y_min, z_min,
                                         x_max, y_max, z_max)
    sprite_bounding_box.print_info()

    # create the volume from sprite
    print ("create the volume from the sprite")
    volume = vg.create_volume_from_sprite(position_array,
                                          intensity_array,
                                          sprite_bounding_box,
                                          volume_largest_side_size)

    # convert the data to positive values, for visualization
    print("* BEGIN [%s] Conversion ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    volume = v_utils.map_volume_to_positive_values(volume)
    end = time.clock()
    print("** DONE MAPPING to 8 bits [%s] in %f" %
          (sys._getframe().f_code.co_name, end - start))

    # map the volume between 0 and 255
    print("* BEGIN [%s] MAPPING ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    volume = v_utils.map_volume_to_ubyte(volume)
    v_utils.write_volume_to_raw_file(volume, "volume")
    end = time.clock()
    print("** DONE MAPPING to 8 bits [%s] in %f" %
          (sys._getframe().f_code.co_name, end - start))

    # write the volume
    print("* BEGIN [%s] WRITING ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    v_utils.write_volume_to_raw_file(volume, "volume")
    end = time.clock()
    print("** DONE WRITING [%s] in %f" % (sys._getframe().f_code.co_name,
                                          end - start))


# create the volume
create_volume_from_sprite(512)




