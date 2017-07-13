#!/usr/bin/python
################################################################################
# Copyright BBP/EPFL (c) 2015
# Author : Marwan Abdellah (marwan.abdellah@epfl.ch)
################################################################################

import sys
import os
import shutil



################################################################################
def clean_and_create_new_directory(path):
    """
    * creates a new directory after removing an existing one with the same name.

    keyword arguments
    :param path: the path to the created directory
    """

    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


################################################################################
def save_file_to_disk(file_path, data_string):
    """
    * saves string to a file on disk

    keyword arguments
    :param file_path: an output path
    :param data_string: the string that will be written to the file
    """

    file_handler = open(file_path, 'w')
    file_handler.write(data_string)
    file_handler.close()


################################################################################
def list_files_in_directory(path, extension):
    """
    * lists all the files that ends with a specific extension in a single
    directory.

    keyword arguments
    :param path: the given directory
    :param extension: the extension of the file
    """

    directory_list = []
    for i_file in os.listdir(path):
        if i_file.endswith(".%s" % extension):
            directory_list.append(i_file)

    return directory_list
