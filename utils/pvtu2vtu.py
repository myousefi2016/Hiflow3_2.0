#!/usr/bin/env python

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

################################################################################
# SCRIPT FOR CONVERSION OF PVTU FILES to sequential VTU FILE
# =================================================
#
# Based on Vtk
#
# Author Michael Schick
################################################################################

import sys
import vtk
from vtk.util.misc import vtkGetDataRoot

filename_read = sys.argv[1]
filename_write = sys.argv[2]

reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName(filename_read)
reader.Update()

#seq_grid = vtk.vtkUnstructuredGrid()
#seq_grid.DeepCopy(reader.GetOutput())

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetDataModeToAscii()
writer.SetFileName(filename_write)
writer.SetInput(reader.GetOutput())
writer.Write()
