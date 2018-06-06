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

/// \author Thomas Gengenbach

#ifndef _CONVDIFF_H_
#    define _CONVDIFF_H_

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>

#    include "hiflow.h"
#    include "config.h"
#    include "common/log.h"

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

static const int DIM = 2;

typedef LADescriptorCoupledD LAD;

typedef std::vector<double> Coord;

int main ( int argc, char** argv );

// ConvectionDiffusionAssembler /////////////

class ConvectionDiffusionAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    ConvectionDiffusionAssembler ( const std::vector<double>& beta, double nu )
    : beta_ ( beta ),
    nu_ ( nu )
    {
    }

    ConvectionDiffusionAssembler ( )
    : beta_ ( 0 ),
    nu_ ( 0.0 )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

    const std::vector<double> beta_;
    const double nu_;
};

// ConvectionDiffusionApp

class ConvectionDiffusionApp
{
  public:
    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;

    ConvectionDiffusionApp ( PropertyTree &param, const std::vector<double>& beta, double nu );
    ~ConvectionDiffusionApp ( );

    void run ( );

  private:
    void read_and_distribute_mesh ( const std::string& filename );

    void prepare_space ( );
    void prepare_linear_solver ( );
    void prepare_lin_alg_structures ( );
    void prepare_bc ( );
    void assemble_system ( );
    void solve_system ( );
    void visualize_solution ( );

    PropertyTree* param_;

    MPI_Comm comm_; // must use concrete object here -> hence C interface
    int rank_;
    int num_partitions_;
    const int master_rank_;

    MeshPtr mesh_;
    VectorSpace<double> space_;
    int refinement_level_;

    ConvectionDiffusionAssembler* local_asm_;

    std::string method_;

    SYSTEM la_sys_;
    Couplings<double> couplings_;
    double rho_, nu_;
    double H_, Um_, B_;
    const std::vector<double> beta_;

    LAD::MatrixType matrix_;
    LAD::VectorType sol_, rhs_;
    std::vector<int> dirichlet_dofs_;
    std::vector<LAD::DataType> dirichlet_values_;

    // would be better to use pointer to base class here, but not
    // possible since these objects do not have sensible constructors...
    GMRES<LAD> solver_;
    // ilupp preconditioner
#    ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
#    endif

    std::string mesh_filename_;
};

#endif /* _CONVDIFF_H_ */
