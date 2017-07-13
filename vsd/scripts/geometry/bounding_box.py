#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import os
import sys
import math

# append the internal scripts directories into the system paths
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("%s/../sprite" % current_path)
sys.path.append("%s/../volume" %current_path)

import geometric_utilities as g_utils

class BoundingBox:
    """
    * A class to encapsulate the bounding box data of the point sprite.
    """
    def __init__(self):
        self.x_coi = 0.0
        self.y_coi = 0.0
        self.z_coi = 0.0
        return

    def __init__(self, x_min, y_min, z_min, x_max, y_max, z_max, x_coi=0.0, y_coi=0.0, z_coi=0.0):
        self.compute_bounding_box(x_min, y_min, z_min, x_max, y_max, z_max, x_coi, y_coi, z_coi)

    def compute_bounding_box(self, x_min, y_min, z_min, x_max, y_max, z_max, x_coi, y_coi, z_coi):
        self.x_min = x_min
        self.y_min = y_min
        self.z_min = z_min
        self.x_max = x_max
        self.y_max = y_max
        self.z_max = z_max
        self.width = x_max - x_min
        self.height = y_max - y_min
        self.depth = z_max - z_min
        self.x_center = x_min + self.width / 2.0
        self.y_center = y_min + self.height / 2.0
        self.z_center = z_min + self.depth / 2.0
        self.x_coi = x_coi
        self.y_coi = y_coi
        self.z_coi = z_coi

        # for scaling with respect to a 'normalized' and 'unit' bounding box
        # later
        self.compute_largest_dimension()
        self.compute_smallest_dimension()

        # compute the 'normalized' bounding box that extends between 0 and 1
        self.compute_normalized_bounding_box()

        # compute the 'unit' bounding box that is centered at the origin and
        # extends between -0.5 and 0.5 at the largest dimension
        self.compute_unit_bounding_box()

    def compute_largest_dimension(self):
        self.largest_dimension = g_utils.compute_largest(self.width,
                                                         self.height,
                                                         self.depth)

    def compute_smallest_dimension(self):
        self.smallest_dimension = g_utils.compute_smallest(self.width,
                                                           self.height,
                                                           self.depth)

    def compute_normalized_bounding_box(self):
        self.normalized_width = self.width / self.largest_dimension
        self.normalized_height = self.height / self.largest_dimension
        self.normalized_depth = self.depth / self.largest_dimension

        # the normalized bounding box starts at the origin and extends to
        # the maximum coordinates, where the largest one is 1.0
        largest_coordinate = g_utils.compute_largest(self.x_max,
                                                     self.y_max,
                                                     self.z_max)
        self.normalized_x_min = 0.0
        self.normalized_y_min = 0.0
        self.normalized_z_min = 0.0
        self.normalized_x_max = self.x_max / largest_coordinate
        self.normalized_y_max = self.y_max / largest_coordinate
        self.normalized_z_max = self.z_max / largest_coordinate

        self.normalized_x_center = self.normalized_width / 2.0
        self.normalized_y_center = self.normalized_height / 2.0
        self.normalized_z_center = self.normalized_depth / 2.0

    def compute_unit_bounding_box(self):
        self.unit_x_max = self.normalized_width / 2.0
        self.unit_y_max = self.normalized_height / 2.0
        self.unit_z_max = self.normalized_depth / 2.0
        self.unit_x_min = -self.normalized_width / 2.0
        self.unit_y_min = -self.normalized_height / 2.0
        self.unit_z_min = -self.normalized_depth / 2.0
        self.unit_x_center = 0.0
        self.unit_y_center = 0.0
        self.unit_z_center = 0.0

    def get_x_min(self):
        return self.x_min

    def get_y_min(self):
        return self.y_min

    def get_z_min(self):
        return self.z_min

    def get_x_max(self):
        return self.x_max

    def get_y_max(self):
        return self.y_max

    def get_z_max(self):
        return self.z_max

    def get_x_center(self):
        return self.x_center

    def get_y_center(self):
        return self.y_center

    def get_z_center(self):
        return self.z_center

    def get_x_coi(self):
        return self.x_coi

    def get_y_coi(self):
        return self.y_coi

    def get_z_coi(self):
        return self.z_coi

    def get_width(self):
        return self.width

    def get_height(self):
        return self.height

    def get_depth(self):
        return self.depth

    def get_normalized_x_min(self):
        return self.normalized_x_min

    def get_normalized_y_min(self):
        return self.normalized_y_min

    def get_normalized_z_min(self):
        return self.normalized_z_min

    def get_normalized_x_max(self):
        return self.normalized_x_max

    def get_normalized_y_max(self):
        return self.normalized_y_max

    def get_normalized_z_max(self):
        return self.normalized_z_max

    def get_normalized_x_center(self):
        return self.normalized_x_center

    def get_normalized_y_center(self):
        return self.normalized_y_center

    def get_normalized_z_center(self):
        return self.normalized_z_center

    def get_normalized_width(self):
        return self.normalized_width

    def get_normalized_height(self):
        return self.normalized_height

    def get_normalized_depth(self):
        return self.normalized_depth

    def get_unit_x_min(self):
        return self.unit_x_min

    def get_unit_y_min(self):
        return self.unit_y_min

    def get_unit_z_min(self):
        return self.unit_z_min

    def get_unit_x_max(self):
        return self.unit_x_max

    def get_unit_y_max(self):
        return self.unit_y_max

    def get_unit_z_max(self):
        return self.unit_z_max

    def get_unit_x_center(self):
        return self.unit_x_center

    def get_unit_y_center(self):
        return self.unit_y_center

    def get_unit_z_center(self):
        return self.unit_z_center

    def get_unit_width(self):
        return self.normalized_width

    def get_unit_height(self):
        return self.normalized_height

    def get_unit_depth(self):
        return self.normalized_depth

    def get_volumetric_grid_dimensions(self, volume_largest_side_size):
        """
        * gets the dimensions of a volume that correspond to this bounding box
        given the resolution of the volume in terms of its largest side.
        """

        # the volume largest side should correspond to the largest sprite
        # dimension.
        volume_width = int(self.normalized_width * volume_largest_side_size)
        volume_height = int(self.normalized_height * volume_largest_side_size)
        volume_depth = int(self.normalized_depth * volume_largest_side_size)

        return volume_width, volume_height, volume_depth

    def convert_coordinate_to_volume_grid_index(self, x, y, z,
                                                volume_largest_side_size):
        """
        * converts the spatial coordinates of a sprite into a volume grid index.
        """

        volume_width, volume_height, volume_depth = \
            self.get_volumetric_grid_dimensions(volume_largest_side_size)

        # find the normalized values of the given coordinates (between 0 and 1).
        x_normalized = (x - self.x_min) / self.width
        y_normalized = (y - self.y_min) / self.height
        z_normalized = (z - self.z_min) / self.depth

        # find the corresponding volume grid coordinates
        x_grid = int(x_normalized * volume_width) - 1
        y_grid = int(y_normalized * volume_height) - 1
        z_grid = int(z_normalized * volume_depth) - 1

        return x_grid, y_grid, z_grid

    def print_original_bounding_box_info(self):
        print("Original bounding box info: "
            "\n\tCenter = %f %f %f "
            "\n\tDimensions = %f %f %f"
            "\n\tPMin = %f %f %f"
            "\n\tPMax = %f %f %f"
            "\n\tLarget dimension = %f "
            "\n\tSmallest dimension = %f" % (self.x_center,
                                             self.y_center,
                                             self.z_center,
                                             self.width,
                                             self.height,
                                             self.depth,
                                             self.x_min,
                                             self.y_min,
                                             self.z_min,
                                             self.x_max,
                                             self.y_max,
                                             self.z_max,
                                             self.largest_dimension,
                                             self.smallest_dimension))

    def print_normalized_bounding_box_info(self):
        print("Normalized bounding box info: "
            "\n\tCenter = %f %f %f "
            "\n\tDimensions = %f %f %f"
            "\n\tPMin = %f %f %f"
            "\n\tPMax = %f %f %f" % (self.normalized_x_center,
                                     self.normalized_y_center,
                                     self.normalized_z_center,
                                     self.normalized_width,
                                     self.normalized_height,
                                     self.normalized_depth,
                                     self.normalized_x_min,
                                     self.normalized_y_min,
                                     self.normalized_z_min,
                                     self.normalized_x_max,
                                     self.normalized_y_max,
                                     self.normalized_z_max))

    def print_unit_bounding_box_info(self):
        print("Unit bounding box info: "
            "\n\tCenter = %f %f %f "
            "\n\tDimensions = %f %f %f"
            "\n\tPMin = %f %f %f"
            "\n\tPMax = %f %f %f" % (self.unit_x_center,
                                     self.unit_y_center,
                                     self.unit_z_center,
                                     self.normalized_width,
                                     self.normalized_height,
                                     self.normalized_depth,
                                     self.unit_x_min,
                                     self.unit_y_min,
                                     self.unit_z_min,
                                     self.unit_x_max,
                                     self.unit_y_max,
                                     self.unit_z_max))

    def print_info(self):
        self.print_original_bounding_box_info()
        self.print_normalized_bounding_box_info()
        self.print_unit_bounding_box_info()
