#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import argparse

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../../geometry" % current_path)
sys.path.append("%s/../../sprite" % current_path)
sys.path.append("%s/../../utilities" %current_path)
sys.path.append("%s/../../volume" %current_path)

import numpy_to_sprite as npsprite
import sprite_generation as sg
import sprite_writer as sw
import bounding_box as bb
import geometric_utilities as g_utils
import file_handling as fh


################################################################################
if __name__ == "__main__":

    # setup command line parser and parse arguments
    parser = argparse.ArgumentParser()
    
    # read the positions from a sample dataset to get the positions of the
    # vsd fluorescene events.
    
    arg_help = 'the directory where the input files will be read from'
    parser.add_argument('--input-directory',
                        action='store',
                        dest='input_directory',
                        help=arg_help)

    arg_help = 'the directory where the out files will be generated'
    parser.add_argument('--output-directory',
                        action='store', default='.',
                        dest='output_directory',
                        help=arg_help)

    # parse the arguments
    args = parser.parse_args() 

    # find all the files in the directory
    print args.input_directory
    files = os.listdir(args.input_directory)
    
    sprite_prefix = "simulation"

    # position array to be given to all the witers
    position_array = []
    # check that position file exists
    try:
        'xyz' in files
    except Exception as e:
        print( "<p>Error: no position file.  %s</p>" % str(e) )

    # load the position file into a numpy array
    position_array = npsprite.load_position_array\
        ("%s/%s" % (args.input_directory, 'xyz'))

    # save the numpy array in a psp format
    sw.write_point_sprite_positions \
        (args.output_directory, sprite_prefix, position_array)
    
    for i_file in files:
        if (i_file == "xyz"):
            
            # skip this file
            print "doing nothing cause this is the position file"
            
        elif (i_file == 'gids'):

            # skip this file
            print "doing nothing cause this is the gids file"

        else:

            # load the intensity file into a numpy array
            intensity_array, time_step = npsprite.load_intensity_array\
                ("%s/%s" % (args.input_directory, i_file))

            print "time step : %s" % time_step

             # save the numpy array in a psp format 
            sprite_prefix_with_time_step = "%s%s" %(sprite_prefix, time_step)
            sw.write_point_sprite_intensity \
                (args.output_directory, sprite_prefix_with_time_step, intensity_array)

            sw.write_point_sprite_header(args.output_directory,
                              psh_file_name=sprite_prefix_with_time_step,
                              psp_file_name=sprite_prefix,
                              psi_file_name=sprite_prefix_with_time_step,
                              position_point_sprite_array=position_array,
                              time_step=time_step)
    
    
