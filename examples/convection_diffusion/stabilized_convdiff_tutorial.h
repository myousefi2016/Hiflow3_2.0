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

/// \author Simon Gawlok, Carmen Straub

#ifndef _STABILIZED_CONVDIFF_TUTORIAL_H_
#    define _STABILIZED_CONVDIFF_TUTORIAL_H_

#    define COMPUTE_HESSIAN

#    define SUPG
//#define GLS

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>
#    include <limits>
#    include <algorithm>
#    include <iostream>
#    include <fstream>

#    include "hiflow.h"

// #define SUPG
// #define GLS

// All names are imported for simplicity
using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

// Dimension of the problem
const int DIM = 2;

// Shorten some datatypes with typedefs
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

typedef std::vector<double> Coord;

int main ( int argc, char** argv );

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Convection-Diffusion Assembler

class ConvectionDiffusionAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    ConvectionDiffusionAssembler ( const std::vector<double>& a,
                                   double nu, double peclet )
    : a_ ( a ),
    nu_ ( nu ),
    peclet_ ( peclet )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        // Compute local matrix.
        const int num_q = num_quadrature_points ( );
        const int total_dofs = num_dofs_total ( );

        lm.Resize ( total_dofs, total_dofs );
        lm.Zeros ( );

        // Loop over all quadrature points q.
        for ( int q = 0; q < num_q; ++q )
        {
            // Compute local matrix
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            { // Loop over test functions
                for ( int j = 0; j < n_dofs; ++j )
                { // Loop over ansatz functions
                    const double prod = a_[0] * grad_phi ( j, q )[0] +
                            a_[1] * grad_phi ( j, q )[1];
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * ( prod * phi ( i, q ) +
                            nu_ * dot ( grad_phi ( j, q ), grad_phi ( i, q ) ) ) * std::abs ( detJ ( q ) );
                }
            }

            // Compute stabilization parameter
            Scalar tau = 0.0;
            tau = ( 1.0 / std::tanh ( peclet_ ) ) - ( 1.0 / peclet_ );

            // Implement the SUPG stabilization
#    ifdef SUPG
            for ( int i = 0; i < n_dofs; ++i )
            {
                const double prod_test = a_[0] * grad_phi ( i, q )[0] +
                        a_[1] * grad_phi ( i, q )[1];
                for ( int j = 0; j < n_dofs; ++j )
                {
                    const double prod_ansatz = a_[0] * grad_phi ( j, q )[0] +
                            a_[1] * grad_phi ( j, q )[1];
                    // Get Laplace of ansatz function as trace of the Hessian matrix
                    double hess = 0.0;
                    for ( int k = 0; k < DIM; ++k )
                    {
                        hess += H_phi ( j, q )( k, k );
                    }
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * ( tau * prod_ansatz *
                            prod_test - tau * nu_ * hess * prod_test ) * std::abs ( detJ ( q ) );
                }
            }
#    endif

            // Implement the GLS stabilization
#    ifdef GLS
            for ( int i = 0; i < n_dofs; ++i )
            {
                const double prod_test = a_[0] * grad_phi ( i, q )[0] +
                        a_[1] * grad_phi ( i, q )[1];

                // Get Laplace of test function as trace of the Hessian matrix
                double hess_test = 0.0;
                for ( int k = 0; k < DIM; ++k )
                {
                    hess_test += H_phi ( i, q )( k, k );
                }
                for ( int j = 0; j < n_dofs; ++j )
                {
                    const double prod = a_[0] * grad_phi ( j, q )[0] +
                            a_[1] * grad_phi ( j, q )[1];
                    // Get Laplace of ansatz function as trace of the Hessian matrix
                    double hess = 0.0;
                    for ( int k = 0; k < 1; ++k )
                    {
                        hess += H_phi ( j, q )( k, k );
                    }

                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * ( tau * prod * prod_test -
                            tau * nu_ * hess * prod_test - tau * nu_ * prod * hess_test + tau *
                            nu_ * nu_ * hess * hess_test ) * std::abs ( detJ ( q ) );
                }
            }

#    endif
        }
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature, LocalVector& lv )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        const int total_dofs = num_dofs_total ( );

        lv.clear ( );
        lv.resize ( total_dofs, 0. );
    }

  private:
    std::vector<double> a_;
    const double nu_;
    const double peclet_;
};

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Convection-Diffusion Application

class ConvectionDiffusionApp
{
  public:
    ConvectionDiffusionApp ( const std::string& param_filename,
                             const std::string& path_mesh );
    ~ConvectionDiffusionApp ( );
    void run ( ); // Run the application.

  private:
    // Member functions.

    // Read and distribute mesh.
    void read_and_distribute_mesh ( );

    // Setup space, linear algebra, and compute Dirichlet values.
    void prepare_system ( );

    // Compute the matrix and rhs.
    void assemble_system ( );

    // Compute solution.
    void solve_system ( );

    // Visualize the results.
    void visualize ( );

    // Member variables.
    PropertyTree params_; // XML tree with parameters.

    MPI_Comm comm_; // MPI communicator.
    int rank_; // Local process rank.
    int num_partitions_; // Number of processes.
    const int master_rank_; // Master rank.

    int refinement_level_; // Refinement level.
    int fe_degree_; // Degree of local shape functions.

    std::string path_mesh; // Path to geometry data.
    MeshPtr mesh_; // Local mesh.
    VectorSpace<double> space_; // Solution space.
    Couplings<double> couplings_; // Linear algebra couplings helper object.
    CoupledMatrix<Scalar>* matrix_; // System matrix.
    CoupledVector<Scalar>* sol_, *rhs_; // Vectors for solution and load vector.

    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> dirichlet_dofs_;
    // Dof values for Dirichlet boundary conditions.
    std::vector<Scalar> dirichlet_values_;

    StandardGlobalAssembler<double> global_asm_; // Global assembler.

    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact; // Solver for the linear system.

    double nu_; // Coefficient of diffusivity.
    std::vector<Scalar> a_; // Convection velocity.
};

#endif /* _STABILIZED_CONVDIFF_TUTORIAL_H_ */
