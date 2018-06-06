#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2011-2017 Vincent Heuveline
#
# HiFlow3 is free software: you can redistribute it and/or modify it under the
# terms of the European Union Public Licence (EUPL) v1.2 as published by the
#/ European Union or (at your option) any later version.
#
# HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
# details.
#
# You should have received a copy of the European Union Public Licence (EUPL) v1.2
# along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

# \author Thomas Gengenbach

################################################################################
#
# SCRIPT FOR CONVERSION OF BINARY VTU FILES TO ASCII VTU-FORMAT
# =============================================================
#
# Based on Vtk
#
################################################################################

import os
import sys
import vtk

# ************************************************************************
# check for correct number of arguments

if len(sys.argv) != 3:
    print 'Usage: binaryvtu2asciivtu.py <src_file.vtu> <dest_file.vtu>'
    sys.exit(1)

# ************************************************************************
# arguments

src_file  = sys.argv[1]
dest_file = sys.argv[2]

# ************************************************************************
# split src_file and set src_postfix to the file-specific ending

tmp = src_file.rpartition(".")
src_postfix = tmp[2]

if (src_postfix == "" or (src_postfix != "vtu")):
    print 'Error: src_file didn:t have a valid file specifier!'
    sys.exit(1)

# split dest_file and set dest_postfix to the file-specific ending
tmp = dest_file.rpartition(".")
dest_postfix = tmp[2]

if (dest_postfix == "" or (dest_postfix != "vtu")):
    print 'Error: dest_file didn:t have a valid file specifier!'
    sys.exit(1)

# ************************************************************************
# unstructures grid structure

ugrid_filter = vtk.vtkAppendFilter()

# ************************************************************************
# read file
if (src_postfix == "vtu"):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(src_file)
    ugrid_filter.SetInput(reader.GetOutput())
    reader.Update()

# ************************************************************************
# write file

if (dest_postfix == "vtu"):
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(dest_file)
    writer.SetInput(ugrid_filter.GetOutput())
    writer.Update()
