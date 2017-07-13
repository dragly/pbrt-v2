#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import time
import argparse

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../geometry" % current_path)
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" % current_path)

import sprite_generation as sg
import bounding_box as bb
import geometric_utilities as g_utils
import volume_generation as vg
import volume_utilities as v_utils
import sprite_reader as s_reader


################################################################################
def create_volume_from_sprite(input_directory, psh_file,
                              volume_largest_side_size):

    # read the psh file and return its data
    psh_file_path = "%s/%s" % (input_directory, psh_file)
    event_count, \
    x_center, y_center, z_center, \
    x_coi, y_coi, z_coi, \
    volume_width, volume_height, volume_depth, \
    psp_file, psi_file, time_step = s_reader.parse_point_sprite_header(psh_file_path)

    # compute the details of the bounding box based on the parameters reterived
    # from the header file
    x_min, y_min, z_min, x_max, y_max, z_max = \
        g_utils.convert_center_and_dimensions_to_bounds \
        (x_center, y_center, z_center, volume_width, volume_height, volume_depth)
    sprite_bounding_box = bb.BoundingBox(x_min, y_min, z_min,
                                         x_max, y_max, z_max)

    # read the psb file and return a list of the vsd fluorescent events
    # positions
    position_array = s_reader.parse_point_sprite_positions \
        (psp_file, event_count)

    # read the psi file and return a list of the vsd fluorescent events
    # intensities
    intensity_array = s_reader.parse_point_sprite_intensities \
        (psi_file, event_count)

    sprite_bounding_box.print_info()

    # create the volume from sprite
    print("* BEGIN [%s] creating grid ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    volume = vg.create_volume_from_sprite(position_array,
                                          intensity_array,
                                          sprite_bounding_box,
                                          volume_largest_side_size)
    end = time.clock()
    print("** DONE creating grid [%s] in %f" %
          (sys._getframe().f_code.co_name, end - start))

    # convert the data to positive values, for visualization
    print("* BEGIN [%s] +ve ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    volume = v_utils.map_volume_to_positive_values(volume)
    end = time.clock()
    print("** DONE +ve [%s] in %f" %
          (sys._getframe().f_code.co_name, end - start))

    # map the volume between 0 and 255
    print("* BEGIN [%s] mapping ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    volume = v_utils.map_volume_to_ubyte(volume)
    v_utils.write_volume_to_raw_file(volume, "volume")
    end = time.clock()
    print("** DONE mapping [%s] in %f" %
          (sys._getframe().f_code.co_name, end - start))

    # write the volume
    print("* BEGIN [%s] writing ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    v_utils.write_volume_to_raw_file(volume, "volume")
    end = time.clock()
    print("** DONE writing [%s] in %f" % (sys._getframe().f_code.co_name,
                                          end - start))


################################################################################
if __name__ == "__main__":

    # setup command line parser and parse arguments
    parser = argparse.ArgumentParser()

    arg_help = 'the path of the point sprite data (positions and intensities)'
    parser.add_argument('--input-directory',
                        action='store',
                        dest='input_directory',
                        help=arg_help)

    arg_help = 'the point sprite header file (meta data of the sprite)'
    parser.add_argument('--psh-file',
                        action='store',
                        dest='psh_file',
                        help=arg_help)

    arg_help = 'the largest dimension of the volume'
    parser.add_argument('--n-side',
                        action='store',
                        dest='n_side',
                        help=arg_help)

    # parse the arguments
    args = parser.parse_args()

    # create the volume from the sprite
    create_volume_from_sprite(args.input_directory, args.psh_file,
                              int(args.n_side))
