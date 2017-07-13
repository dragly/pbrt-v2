#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import time
import sys

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" %current_path)


################################################################################
def compute_largest(x, y, z):
    """
    * computes the largest value of a given set of three numbers.

    keyword arguments
    :param x: first value
    :param y: second value
    :param z : third value
    :return: the largest value of the given set of values.
    """

    largest = x
    if(y > largest):
        largest = y
    if(z > largest):
        largest = z

    return float(largest)


################################################################################
def compute_smallest(x, y, z):
    """
    * computes the smallest value of a given set of three numbers.

    keyword arguments
    :param x: first value
    :param y: second value
    :param z : third value
    :return: the smallest value of the given set of values.
    """

    smallest = x
    if(y < smallest):
        smallest = y
    if(z < smallest):
        smallest = z

    return float(smallest)


################################################################################
def compute_bounds_of_sprite(position_array):
    """
    * computes the spatial bounds of the point sprite positions

    * keyword arguments
    :param position_array:  a three-dimensional numpy array with the actual
           positions of the point sprite.
    :return: the spatial bounds of the sprite in x, y and z in the form of
            [x_min, y_min, z_min, x_max, y_max, z_max].
    """
    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    
    # initial values at infinity and -infinity
    x_min = sys.float_info.max
    y_min = sys.float_info.max
    z_min = sys.float_info.max
    x_max = -sys.float_info.max
    y_max = -sys.float_info.max
    z_max = -sys.float_info.max

    # iterate and find the actual minimum and maximum coordinates
    for i_point in position_array:
        x_pos = i_point[0]
        y_pos = i_point[1]
        z_pos = i_point[2]

        if (x_pos < x_min):
            x_min = x_pos
        if (y_pos < y_min):
            y_min = y_pos
        if (z_pos < z_min):
            z_min = z_pos

        if (x_pos > x_max):
            x_max = x_pos
        if (y_pos > y_max):
            y_max = y_pos
        if (z_pos > z_max):
            z_max = z_pos

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    return x_min, y_min, z_min, x_max, y_max, z_max


################################################################################
def convert_center_and_dimensions_to_bounds(x_center, y_center, z_center,
                                            width, height, depth):
    """
    * compute the minima and maxima of the bounding box given by its center
    and dimensions

    * keyword arguments
    :param x_center: the x-coordinate of the center of the bounding box
    :param y_center: the y-coordinate of the center of the bounding box
    :param z_center: the z-coordinate of the center of the bounding box
    :param width: the width of the bounding box
    :param height: the height of the bounding box
    :param depth: the depth of the bounding box 

    :return: the spatial bounds of the sprite in x, y and z in the form of
            [x_min, y_min, z_min, x_max, y_max, z_max].
    """

    x_min = x_center - (width / 2.0)
    y_min = y_center - (height / 2.0)
    z_min = z_center - (depth / 2.0) 
    x_max = x_center + (width / 2.0)
    y_max = y_center + (height / 2.0)
    z_max = z_center + (depth / 2.0)

    return x_min, y_min, z_min, x_max, y_max, z_max
