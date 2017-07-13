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

import struct
import random
import time
import numpy
import sys
from numpy import *

###############################################################################
def generate_random_point_sprite_positions(number_points=100,
                                           x_min=0.0, y_min=0.0, z_min=0.0,
                                           x_max=1.0, y_max=1.0, z_max=1.0):
    """
    * creates a random point sprite positions as a numpy three-dimensional array.

    keyword arguments
    :param number_points: the number of points in the sprite.
    :param x_min: minimum x coordinate of the sprite.
    :param y_min: minimum y coordinate of the sprite.
    :param z_min: minimum z coordinate of the sprite.
    :param x_max: maximum x coordinate of the sprite.
    :param y_max: maximum y coordinate of the sprite.
    :param z_max: maximum z coordinate of the sprite.

    :return: a three-dimensional numpy array with random positions of the point sprite.
    """
    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # create a random position array, and fill it with random data between the
    # given spatial range
    array = numpy.zeros((number_points, 3))
    for i_point in array:
        i_point[0] = random.uniform(x_min, x_max)
        i_point[1] = random.uniform(y_min, y_max)
        i_point[2] = random.uniform(z_min, z_max)

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    return array


###############################################################################
def generate_random_point_sprite_intensities(number_points=100,
                                             intensity_min=0.0,
                                             intensity_max=100.0):
    """
    * creates a random point sprite intensities as a numpy one-dimensional array.

    keyword arguments
    :param number_points: the number of points in the sprite.
    :param intensity_min: minimum power value of the sprite.
    :param intensity_max: maximum power value of the sprite.

    :return: a one-dimensional numpy array with random intensities of the point sprite.
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()
    
    # create a random position array, and fill it with random data between the
    # given spatial range
    array = numpy.zeros((number_points, 1))
    for i_point in array:
        i_point[0] = random.uniform(intensity_min, intensity_max)

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))
    
    return array


###############################################################################
def generate_random_point_sprite(number_points=100,
                                 x_min=0.0, y_min=0.0, z_min=0.0,
                                 x_max=1.0, y_max=1.0, z_max=1.0,
                                 intensity_min=0.0, intensity_max=100.0):
    """
    * creates a random point sprite with two numpy arrays, a
    three-dimensional one representing the positions and a one-dimensional
    array reflecting the intensities.

     keyword arguments
    :param number_points: the number of points in the sprite.
    :param x_min: minimum x coordinate of the sprite.
    :param y_min: minimum y coordinate of the sprite.
    :param z_min: minimum z coordinate of the sprite.
    :param x_max: maximum x coordinate of the sprite.
    :param y_max: maximum y coordinate of the sprite.
    :param z_max: maximum z coordinate of the sprite.
    :param intensity_min: minimum power value of the sprite.
    :param intensity_max: maximum power value of the sprite.

    :return two numpy arrays, one representing the positions of the sprite
    and the other one represents the intensities of the points in the sprite.
    """

    print("* BEGIN [%s] ..." % sys._getframe().f_code.co_name)
    start = time.clock()

    # create the random positions array
    position_array = generate_random_point_sprite_positions(number_points,
            x_min, y_min, z_min, x_max, y_max, z_max)

    # create the random intensities array
    intensity_array =  generate_random_point_sprite_intensities(number_points,
            intensity_min, intensity_max)

    end = time.clock()
    print("** DONE [%s] in %f" % (sys._getframe().f_code.co_name, end - start))

    return position_array, intensity_array

