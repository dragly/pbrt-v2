#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import struct
import time
import ntpath
from numpy import *


################################################################################
# @parse_point_sprite_header
################################################################################
def parse_point_sprite_header(psh_file_path, sprite_type="psp-psi"):
    """
    Parses a point-sprite header file and returns its parameters.

    :param psh_file_path: the path of the psh file where the sprite meta data
    will be retrieved
    :param sprite_type: either "ps" which means that all the data of the sprite
    including the xyz-position and the intensities are integrated in the same file
    or "psp-psi" where the position data is saved to a psp file and the
    intensities are saved in psi file.
    :return: header file parameters
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    events_count = 0    # number of events in the point sprite
    x_center = 0.0      # sprite bounding box x-center
    y_center = 0.0      # sprite bounding box y-center
    z_center = 0.0      # sprite bounding box z-center
    x_coi = 0.0         # x-center of interest of circuit
    y_coi = 0.0         # y-center of interest of circuit
    z_coi = 0.0         # z-center of interest of circuit
    width = 0.0         # sprite bounding box width
    height = 0.0        # sprite bounding box height
    depth = 0.0         # sprite bounding box depth
    ps_file_name = ""   # sprite file name
    psp_file_name = ""  # position file name
    psi_file_name = ""  # intensity file name
    time_step = 0       # time step

    with open(psh_file_path, 'r') as psh_file:
        for i_line in psh_file:
            value = i_line.split("=")[1]
            if "EventsCount" in i_line:
                events_count = int(value.strip(" "))
            if "XCenter" in i_line:
                x_center = float(value.strip(" "))
            if "YCenter" in i_line:
                y_center = float(value.strip(" "))
            if "ZCenter" in i_line:
                z_center = float(value.strip(" "))
            if "XCOI" in i_line:
                x_coi = float(value.strip(" "))
            if "YCOI" in i_line:
                y_coi = float(value.strip(" "))
            if "ZCOI" in i_line:
                z_coi = float(value.strip(" "))
            if "AABBWidth" in i_line:
                width = float(value.strip(" "))
            if "AABBHeight" in i_line:
                height = float(value.strip(" "))
            if "AABBDepth" in i_line:
                depth = float(value)
            if "VSDFile" in i_line:
                ps_file_name = value.strip("\n")
            if "VSDPositionFile" in i_line:
                psp_file_name = value.strip("\n")
            if "VSDIntensityFile" in i_line:
                psi_file_name = value.strip("\n")
            if "TimeStep" in i_line:
                time_step = value.strip("\n")

    # extract the full path of the psb from the file prefix
    path_prefix = psh_file_path.strip(ntpath.basename(psh_file_path))
    psp_file_path = path_prefix + psp_file_name
    psi_file_path = path_prefix + psi_file_name
    ps_file_path = path_prefix + ps_file_name

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    if(sprite_type == "psp-psi"):
        return (events_count, x_center, y_center, z_center, x_coi, y_coi, z_coi,
                width, height, depth, psp_file_path, psi_file_path, time_step)
    else:
        return (events_count, x_center, y_center, z_center, x_coi, y_coi, z_coi,
                width, height, depth, ps_file_path, time_step)


################################################################################
# @parse_point_sprite_positions
################################################################################
def parse_point_sprite_positions(psp_file_path, event_count):
    """
    Parses a point-sprite positions file and returns an array with these
    positions.

    :param psp_file_path: the path of the psb file where the sprite position
    data are stored.
    :param event_count: number of events in the sprite.
    :return: sprite position data
    """
    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # open the point sprite positions file and read the respective xyz
    # coordinates into an array.
    position_array = []
    print psp_file_path
    psp_file = open(psp_file_path, 'rb')
    for i_values in range(0,  event_count):
        x = float(struct.unpack('f', psp_file.read(4))[0])
        y = float(struct.unpack('f', psp_file.read(4))[0])
        z = float(struct.unpack('f', psp_file.read(4))[0])
        position_array.append([x, y, z])
    psp_file.close()

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    return position_array

################################################################################
# @parse_point_sprite_intensities
################################################################################
def parse_point_sprite_intensities(psi_file_path, event_count):
    """
    Parses a point-sprite positions file and returns an array with these
    positions.

    :param psi_file_path: the path of the psi file where the sprite intensities
    are stored.
    :param event_count: number of events in the sprite.
    :return: sprite intensities
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # open the point sprite positions file and read the respective xyz
    # coordinates into an array.
    intensity_array = []
    psi_file = open(psi_file_path)
    for i_values in range(0,  event_count):
        intensity = float(struct.unpack('f', psi_file.read(4))[0])
        intensity_array.append(intensity)
    psi_file.close()

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    return intensity_array


################################################################################
# @parse_point_sprite
################################################################################
def parse_point_sprite(ps_file_path, event_count):
    """
    Parses a point-sprite file and returns an array with the positions and
    another one with intensities.

    :param psp_file_path: the path of the psb file where the sprite position
    data are stored.
    :param event_count: number of events in the sprite.
    :return: sprite positions, sprite intensites
    """
    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # open the point sprite positions file and read the respective xyz
    # coordinates into an array.
    position_array = []
    intensity_array = []
    ps_file = open(ps_file_path, 'rb')
    for i_values in range(0,  event_count):
        x = float(struct.unpack('f', ps_file.read(4))[0])
        y = float(struct.unpack('f', ps_file.read(4))[0])
        z = float(struct.unpack('f', ps_file.read(4))[0])
        i = float(struct.unpack('f', ps_file.read(4))[0])
        position_array.append([x, y, z])
        intensity_array.append(i)
    ps_file.close()

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    return position_array, intensity_array

################################################################################
# @get_time_step_from_psh
################################################################################
def get_time_step_from_psh(psh_file_path):
    """
    Gets the time step from the point sprite header.
    :param psh_file_path:
    :return:
    """
    with open(psh_file_path, 'r') as psh_file:
        for i_line in psh_file:
            if "TimeStep" in i_line:    
                time_step = i_line.strip("\n")
                time_step = time_step.split("=")[1]
                return time_step 
