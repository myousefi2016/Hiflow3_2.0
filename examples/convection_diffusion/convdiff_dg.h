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

// DGJumpMatrixAssembler /////////////

class DGJumpMatrixAssembler : DGAssemblyAssistant<DIM, double>
{
  public:

    DGJumpMatrixAssembler ( const std::vector<double>& beta, double nu, double theta, double gamma )
    : beta_ ( beta ),
    nu_ ( nu ), theta_ ( theta ), gamma_ ( gamma )
    {
    }

    DGJumpMatrixAssembler ( )
    : beta_ ( 0 ),
    nu_ ( 0.0 ), theta_ ( 0.0 ), gamma_ ( 0.0 )
    {
    }

    void operator() ( const Element<double>& left_elem, const Element<double>& right_elem,
            const Quadrature<double>& left_quad, const Quadrature<double>& right_quad,
            int left_facet_number, int right_facet_number,
            InterfaceSide left_if_side, InterfaceSide right_if_side,
            LocalMatrix& lm )
    {
        const bool is_boundary = ( right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY );

        if ( is_boundary )
        {
            lm.Zeros ( );
            return; // do nothing.
        }

        initialize_for_interface ( left_elem, right_elem,
                                   left_quad, right_quad,
                                   left_facet_number, right_facet_number,
                                   left_if_side, right_if_side );

        const int num_q = num_quadrature_points ( );
        const int num_dofs_trial = trial ( ).num_dofs ( 0 );
        const int num_dofs_test = test ( ).num_dofs ( 0 );

        assert ( lm.nrows ( ) >= num_dofs_trial );
        assert ( lm.ncols ( ) >= num_dofs_test );

        const double C = is_boundary ? 1 : 0.5; // TODO: not needed?

        const double h = std::pow ( std::abs ( trial ( ).detJ ( 0 ) ), 1. / double(DIM ) );
        const double p = left_elem.get_fe_type ( 0 )->get_fe_deg ( );
        const double alpha = p * p / ( std::pow ( h, 2. ) );

        // For convective term
        Vec<DIM, double> beta;
        beta[0] = beta_[0];
        beta[1] = beta_[1];

        for ( int q = 0; q < num_q; ++q )
        {
            const double dS = std::abs ( ds ( q ) );
            const double W = w ( q ) * dS;

            const double dot_n = dot ( test ( ).n ( q ), trial ( ).n ( q ) );

            for ( int i = 0; i < num_dofs_test; ++i )
            {
                for ( int j = 0; j < num_dofs_trial; ++j )
                {
                    lm ( test ( ).dof_index ( i, 0 ), trial ( ).dof_index ( j, 0 ) ) +=
                            W * (
                            // Contribution from diffusive term
                            -C * nu_ * dot ( test ( ).grad_phi ( i, q ), trial ( ).n ( q ) ) * trial ( ).phi ( j, q )
                            + C * theta_ * nu_ * dot ( trial ( ).grad_phi ( j, q ), test ( ).n ( q ) ) * test ( ).phi ( i, q )
                            + gamma_ * alpha * test ( ).phi ( i, q ) * trial ( ).phi ( j, q ) * dot_n
                            // Contribution from convective term
                            + gamma_ * dot ( beta, test ( ).n ( q ) ) * trial ( ).phi ( j, q ) * test ( ).phi ( i, q ) * dot_n
                            );
                }
            }
        }
    }

  private:
    const std::vector<double> beta_;
    const double nu_;
    const double theta_;
    const double gamma_;
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
    //DG parameters
    double dg_theta_, dg_gamma_;

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
