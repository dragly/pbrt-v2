#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

from pylab import imshow, show, get_cmap
from numpy import random

################################################################################
def map_intensity_to_color(intensity):
    """
    * maps the relative intensity values of the sprite (0-1) to colors defined
    by the Spectral color map

    * keywords
    :param intensity: the intensity value of the point in the sprite
    :return color
    """

    # use the Spectral color map, in full range
    color_map = get_cmap("Spectral")

    # map the range 0-255 to 0-100
    mapped_value = int(intensity * 256.0)
    
    return color_map(mapped_value) 


