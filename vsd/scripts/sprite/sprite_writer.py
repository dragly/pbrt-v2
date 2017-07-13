#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import struct
import random
import numpy
import os
import time
import sys
import ntpath
from numpy import *

def get_sprite_header_data(position_point_sprite_array=None):

    events_count = 0    # number of events in the point sprite
    x_center = 0.0      # sprite bounding box x-center
    y_center = 0.0      # sprite bounding box y-center
    z_center = 0.0      # sprite bounding box z-center
    width = 0.0         # sprite bounding box width
    height = 0.0        # sprite bounding box height
    depth = 0.0         # sprite bounding box depth
    psp_file_name = ""  # position file name
    psi_file_name = ""  # intensity file name

    # pre-process the data to extract the bounding box parameters, before
    # writing them to file.
    x_min = sys.float_info.max
    y_min = sys.float_info.max
    z_min = sys.float_info.max
    x_max = -sys.float_info.max
    y_max = -sys.float_info.max
    z_max = -sys.float_info.max

    for i_point in position_point_sprite_array:
        x_pos = i_point[0]
        y_pos = i_point[1]
        z_pos = i_point[2]

        if x_pos < x_min:
            x_min = x_pos
        if y_pos < y_min:
            y_min = y_pos
        if z_pos < z_min:
            z_min = z_pos

        if x_pos > x_max:
            x_max = x_pos
        if y_pos > y_max:
            y_max = y_pos
        if z_pos > z_max:
            z_max = z_pos

        events_count += 1

    # find the bounding box dimensions
    width = x_max - x_min
    height = y_max - y_min
    depth = z_max - z_min

    # the center of interest
    x_coi = numpy.mean(position_point_sprite_array[:, 0])
    y_coi = numpy.mean(position_point_sprite_array[:, 1])
    z_coi = numpy.mean(position_point_sprite_array[:, 2])

    x_center = x_min + width / 2.0
    y_center = y_min + height / 2.0
    z_center = z_min + depth / 2.0

    return events_count, x_center, y_center, z_center, x_coi, y_coi, z_coi, width, height, depth


################################################################################
def write_point_sprite_header(file_path=".",
                              psh_file_name="sprite_header",
                              psp_file_name="sprite_position",
                              psi_file_name="sprite_intensity",
                              position_point_sprite_array=None,
                              time_step=0):
    """
    * writes the header (or meta-data) of the point sprite to a .psh (point
    sprite header) file.

    * keyword arguments:
    :param file_path: the directory where the file will be written to.
    :param psp_file_name: the prefix (or name) of the position file.
    :param psi_file_name: the prefix (or name) of the intensity file.
    :param position_point_sprite_array: the array that carries the values of the
    positions of the point sprite at the respective compartments.
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    events_count, \
    x_center, y_center, z_center, x_coi, y_coi, z_coi, \
    width, height, depth = get_sprite_header_data(position_point_sprite_array)

    # write the data to the file
    psh_file_path = "%s/%s.psh" % (file_path, psh_file_name)
    psh_file = open(psh_file_path, 'w')

    # @note : the center is always set to the origin
    psh_file.write("EventsCount=%d\n" % events_count)
    psh_file.write("XCenter=%f\n" % x_center)
    psh_file.write("YCenter=%f\n" % y_center)
    psh_file.write("ZCenter=%f\n" % z_center)
    psh_file.write("XCOI=%f\n" % x_coi)
    psh_file.write("YCOI=%f\n" % y_coi)
    psh_file.write("ZCOI=%f\n" % z_coi)
    psh_file.write("AABBWidth=%f\n" % width)
    psh_file.write("AABBHeight=%f\n" % height)
    psh_file.write("AABBDepth=%f\n" % depth)

    psp_file = "%s.psp" % psp_file_name
    psi_file = "%s.psi" % psi_file_name
    psh_file.write("VSDPositionFile=%s\n" % psp_file)
    psh_file.write("VSDIntensityFile=%s\n" % psi_file)
    psh_file.write("TimeStep=%s\n" % str(time_step))

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))


################################################################################
def write_point_sprite_positions(file_path=".", file_name="position",
                                 position_point_sprite_array=None):
    """
    * writes the position array of the VSD data to a .psp (point sprite
    position) file.

    * keyword arguments:
    :param file_path: the directory where the file will be written to.
    :param file_name: the prefix (or name) of the file.
    :param position_point_sprite_array: the array that carries the values of the
    positions of the point sprite at the respective compartments.
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # open the position point sprite file in a binary mode
    position_ps_file_path = "%s/%s.psp" % (file_path, file_name)
    position_ps_file = open(position_ps_file_path, 'wb')

    events_count = 0
    for i_point in position_point_sprite_array:
        events_count += 1

    # iterate and fill the file with the position data from the numpy array
    number_events = 0
    for i_point in position_point_sprite_array:
        position_ps_file.write(struct.pack('f', (i_point[0])))
        position_ps_file.write(struct.pack('f', (i_point[1])))
        position_ps_file.write(struct.pack('f', (i_point[2])))
        number_events += 1

    # close the position point sprite file
    position_ps_file.close()

    print("[%d] events have been written to %s" % (number_events,
                                                   position_ps_file_path))
    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))


################################################################################
def write_point_sprite_intensity(file_path=".", file_name="intensity",
                                 intensity_point_sprite_array=None):
    """
    * writes the intensity array of the VSD data to a .psi (point sprite
    intensity) file.

    * keyword arguments:
    :param file_path: the directory where the file will be written to.
    :param file_name: the prefix (or name) of the file.
    :param intensity_point_sprite_array: the array that carries the values of the
    intensities of the point sprite at the respective compartments.
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # open the position point sprite file in a binary mode
    intensity_ps_file_path = "%s/%s.psi" % (file_path, file_name)
    intensity_ps_file = open(intensity_ps_file_path, 'wb')

    # iterate and fill the file with the position data from the numpy array
    number_events = 0
    for i_point in intensity_point_sprite_array:
        intensity_ps_file.write(struct.pack('f', i_point))
        number_events += 1

    # close the position point sprite file
    intensity_ps_file.close()

    print("[%d] events have been written to %s" % (number_events,
                                                   intensity_ps_file_path))

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))


################################################################################
def write_point_sprite_to_file(file_path=".", file_name="sprite",
                               file_format="ascii",
                               position_array=None, intensity_array=None,
                               time_step=0):
    """
    * writes the point sprite in the form of [x y z i] to a file including.
    * if the ascii format is selected, all the data will be written to the same
    file, and if the binary format is selected, the meta-data will be written
    to a psh file and the sprite itself will be created in a .ps file.

    * keyword arguments:
    :param file_path: the directory where the file will be written to.
    :param file_name: the prefix (or name) of the file.
    :param file_format: either ascii or binary.
    :param position_array: the array that carries the values of the
           positions of the point sprite at the respective compartments.
    :param intensity_array: the array that carries the values of the
           intensities of the point sprite at the respective compartments.
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    if(file_format == "ascii"):
        ps_file_path = "%s/%s.sprite" % (file_path, file_name)
        ps_file = open(ps_file_path, 'w')

        # retreive the sprite meta-data and write them to the file
        events_count, \
        x_center, y_center, z_center, \
        x_coi, y_coi, z_coi, \
        width, height, depth =  get_sprite_header_data(position_array)

        ps_file.write("EventsCount=%d\n" % events_count)
        ps_file.write("XCenter=%f\n" % x_center)
        ps_file.write("YCenter=%f\n" % y_center)
        ps_file.write("ZCenter=%f\n" % z_center)
        ps_file.write("XCOI=%f\n" % x_coi)
        ps_file.write("YCOI=%f\n" % y_coi)
        ps_file.write("ZCOI=%f\n" % z_coi)
        ps_file.write("AABBWidth=%f\n" % width)
        ps_file.write("AABBHeight=%f\n" % height)
        ps_file.write("AABBDepth=%f\n" % depth)
        ps_file.write("TimeStep=%f\n" % time_step)

        # iterate and fill the file with the position data from the numpy array
        number_events = 0
        if((position_array != None) and (intensity_array != None)):
            for i_position, i_intensity in zip(position_array, intensity_array):
                ps_file.write("[%f %f %f] [%f]\n" % (i_position[0],
                                                     i_position[1],
                                                     i_position[2],
                                                     i_intensity))
                number_events += 1

        # close the position point sprite file
        ps_file.close()
    else:
        psh_file_path = "%s/%s.psh" % (file_path, file_name)
        ps_file_path = "%s/%s.ps" % (file_path, file_name)

        psh_file = open(psh_file_path, 'w')

        # retreive the sprite meta-data and write them to the file
        events_count, \
        x_center, y_center, z_center, \
        x_coi, y_coi, z_coi, \
        width, height, depth =  get_sprite_header_data(position_array)

        psh_file.write("EventsCount=%d\n" % events_count)
        psh_file.write("XCenter=%f\n" % x_center)
        psh_file.write("YCenter=%f\n" % y_center)
        psh_file.write("ZCenter=%f\n" % z_center)
        psh_file.write("XCOI=%f\n" % x_coi)
        psh_file.write("YCOI=%f\n" % y_coi)
        psh_file.write("ZCOI=%f\n" % z_coi)
        psh_file.write("AABBWidth=%f\n" % width)
        psh_file.write("AABBHeight=%f\n" % height)
        psh_file.write("AABBDepth=%f\n" % depth)
        psh_file.write("VSDFile=%s.psi\n" % file_name)
        psh_file.write("TimeStep=%f\n" % time_step)
        psh_file.close()

        ps_file = open(ps_file_path, 'wb')
        # iterate and fill the file with the from the numpy arrays
        number_events = 0
        if((position_array != None) and (intensity_array != None)):
            for i_position, i_intensity in zip(position_array, intensity_array):
                ps_file.write(struct.pack('f', (i_position[0])))
                ps_file.write(struct.pack('f', (i_position[1])))
                ps_file.write(struct.pack('f', (i_position[2])))
                ps_file.write(struct.pack('f', (i_intensity)))
                number_events += 1

        # close the position point sprite file
        ps_file.close()

    print("[%d] events have been written to %s" % (number_events, ps_file_path))

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))
