// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
// details.
//
// You should have received a copy of the European Union Public Licence (EUPL) v1.2
// along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#ifndef MET_FLOW_VARS_H
#    define MET_FLOW_VARS_H

#    include "linear_algebra/la_descriptor.h"

///
/// \file met_flow_vars.h
/// \brief global variables and typedefs in met_flow application
///
/// \author Philipp Gerstner
///

// Assembler options
#    define OPT_ASM0

// with radial tehd for eigenvalue analysis
#    define nWANT_RADIAL

// Dimension
#    define APPDIM 2

// use *_naturalconvection.xml
#    define nNAT_CONV

// 0: Compute Initial Cond, 1: IC+ primal run 2: IC + DWR run 3:
#    define RUN_MODE2

// select *_development.xml or *_production.xml
#    define nDEVEL 

// Use Q1+Q0 pressure ansatz
#    define nAUGMENT_PRESS

// Use nested schur solver for Boussinesq TEHD
#    define nSCHUR_SOLVER

// Use Hypre Linear Algebra
#    define nUSE_HYPRE

// Use PXest mesh implementation
#    define USE_PXEST

// U=0 on top wall
#    define NOSLIP_TOP

// Use parmetis for initial mesh partitioning
#    define nPARALLEL_PARTITIONING

// Use roating frame of reference
#    define nROTATING_FOR

// Use periodic BC on top and bottom wall
#    define nINFINITE_CYLINDER

// Use
#    define nCART_COORD

// Nondimensionalize system for radial_tehd
#    define nNONDIM

//
#    define ADD_DEP_TERM

// Width of ghost layer (in number of cells), applies only for pXest mesh
// NOTE: for local refinement, this has to be 2 !!
const static int GHOST_LAYER_WIDTH = 2;

// Number of variables and space dimension
const static int DIM = 2;
const static int VARDIM = 1;

// used namespaces
using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;
using namespace hiflow::la;

// Linear algebra typedefs
#    ifndef WITH_HYPRE
#        undef USE_HYPRE
#    endif

#    ifndef WITH_P4EST
#        undef USE_PXEST
#    endif

// mesh implementation
#    ifdef USE_PXEST
const static mesh::IMPL MESHIMPL = IMPL_P4EST;
#    else
const static mesh::IMPL MESHIMPL = IMPL_DBVIEW;
#    endif

#    ifdef USE_HYPRE
typedef LADescriptorHypreD LAD;
#    else
typedef LADescriptorCoupledD LAD;
#    endif

typedef double DATATYPE;
typedef la::SeqDenseMatrix<double> LocalMatrix;
typedef std::vector<double> LocalVector;
typedef std::vector<double> Coord;

// Global variables
const static int MASTER_RANK = 0;
const static double pi = 3.14159265;
const static double eps0 = 8.854187817e-12;

#endif
