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

#ifndef HIFLOW_YALAPLACE_H
#    define HIFLOW_YALAPLACE_H

/// \brief This program demonstrates the use of the AssemblyAssistant
/// to solve a simple Laplace problem with Neumann boundary
/// conditions. It also contains code to compute the L2 and H1 errors,
/// given the analytical solution of the problem. It computes the
/// solution on a sequence of uniformly refined meshes and outputs the
/// convergence rate (error quotient). This error quotient can be
/// compared to the theoretical value (p for H1 and p+1 for L2, where
/// p is the degree of the local basis functions).

/// \author Staffan Ronnas

#    include "hiflow.h"

#    include <cmath>
#    include <string>
#    include <vector>
#    include <algorithm>
#    include <mpi.h>

#    define DIMENSION 2
//#define WEAK_BOUNDARY_CONDITIONS

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::doffem;
using namespace hiflow::mesh;

static const char* DATADIR = MESHES_DATADIR;
static const int DEBUG_LEVEL = 0;

typedef LADescriptorCoupledD LAD;

typedef std::vector<double> Coord;

class LaplaceApp
{
  public:
    static const int DIM = DIMENSION;
    static const int MAX_REFINEMENT_LEVEL = 4;
    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;

    enum ErrorNorm
    {
        L2error, H1error
    };

    LaplaceApp ( int degree );
    ~LaplaceApp ( );

    void read_mesh ( const std::string& filename );
    void distribute_mesh ( );
    void prepare_system ( );
    void prepare_bc ( );
    void assemble_system ( );
    void solve_system ( );
    double compute_error ( ErrorNorm norm );
    void refine_mesh ( );
    void write_visualization ( );
    bool is_done ( ) const;

    int refinement_level ( ) const
    {
        return refinement_level_;
    }
    void debug_output ( );

    int get_rank ( ) const
    {
        return rank_;
    }

#    if 0

    struct ExactSol
    {

        double operator() ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;
            const double x = pt[0];
            const double y = pt[1];
#        if DIMENSION == 2
            return std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y );
#        else
            const double z = pt[2];
            return std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::cos ( 2. * pi * z );
#        endif
        }

        Vec<DIM, double> eval_grad ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;

            Vec<DIM, double> grad;
            const double x = pt[0];
            const double y = pt[1];
            const double z = pt[2];

#        if DIMENSION == 2
            grad[0] = -2. * pi * std::sin ( 2. * pi * x ) * std::cos ( 2. * pi * y );
            grad[1] = -2. * pi * std::cos ( 2. * pi * x ) * std::sin ( 2. * pi * y );
#        else
            const double z = pt[2];
            grad[0] = -2. * pi * std::sin ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::cos ( 2. * pi * z );
            grad[1] = -2. * pi * std::cos ( 2. * pi * x ) * std::sin ( 2. * pi * y ) * std::cos ( 2. * pi * z );
            grad[2] = -2. * pi * std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::sin ( 2. * pi * z );
#        endif
            return grad;
        }
    };
#    endif

    struct ExactSol
    {

        double operator() ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;
            const double x = pt[0];
            const double y = pt[1];
#    if DIMENSION == 2
            return std::sin ( 2. * pi * x ) * std::sin ( 2. * pi * y );
#    else
            const double z = pt[2];
            return std::sin ( 2. * pi * x ) * std::sin ( 2. * pi * y ) * std::sin ( 2. * pi * z );
#    endif
        }

        Vec<DIM, double> eval_grad ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;

            Vec<DIM, double> grad;
            const double x = pt[0];
            const double y = pt[1];

#    if DIMENSION == 2
            grad[0] = 2. * pi * std::cos ( 2. * pi * x ) * std::sin ( 2. * pi * y );
            grad[1] = 2. * pi * std::sin ( 2. * pi * x ) * std::cos ( 2. * pi * y );
#    else
            const double z = pt[2];
            grad[0] = 2. * pi * std::cos ( 2. * pi * x ) * std::sin ( 2. * pi * y ) * std::sin ( 2. * pi * z );
            grad[1] = 2. * pi * std::sin ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::sin ( 2. * pi * z );
            grad[2] = 2. * pi * std::sin ( 2. * pi * x ) * std::sin ( 2. * pi * y ) * std::cos ( 2. * pi * z );
#    endif
            return grad;

        }
    };

    class LocalLaplaceAssembler : private AssemblyAssistant<LaplaceApp::DIM, double>
    {
      public:

        LocalLaplaceAssembler ( const std::vector<int>& dirichlet_dofs, const std::vector<double>& dirichlet_values )
        : dirichlet_dofs_ ( dirichlet_dofs ),
        dirichlet_values_ ( dirichlet_values )
        {
        }

        LocalLaplaceAssembler ( )
        : dirichlet_dofs_ ( 0 ),
        dirichlet_values_ ( 0 )
        {
        }

        void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
        {
            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
            const int num_q = num_quadrature_points ( );
            const int num_local_dofs = num_dofs ( 0 );

            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = std::abs ( detJ ( q ) );

                for ( int i = 0; i < num_local_dofs; ++i )
                {
                    for ( int j = 0; j < num_local_dofs; ++j )
                    {
                        lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                                wq * dot ( grad_phi ( i, q ), grad_phi ( j, q ) ) * dJ;
                    }
                }
            }
        }
        //MATRIX BOUNDARY CONDITIONS

        void operator() ( const Element<double>& element, int facet_number, const Quadrature<double>& quadrature, LocalMatrix& lm )
        {
#    ifdef WEAK_BOUNDARY_CONDITIONS
            AssemblyAssistant<DIM, double>::initialize_for_facet ( element, quadrature, facet_number );
            if ( !dirichlet_dofs_.empty ( ) )
            {
                const int num_q = num_quadrature_points ( );
                const int num_local_dofs = num_dofs ( 0 );

                for ( int q = 0; q < num_q; ++q )
                {
                    const double wq = w ( q );
                    const double dJ = std::abs ( detJ ( q ) );
                    const double alpha = nitsche_regularization ( q );

                    for ( int i = 0; i < num_local_dofs; ++i )
                    {
                        for ( int j = 0; j < num_local_dofs; ++j )
                        {
                            lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                                    wq * ( phi ( i, q ) * ( alpha * phi ( j, q ) - dot ( n ( q ), grad_phi ( j, q ) ) ) - dot ( n ( q ), grad_phi ( i, q ) ) * phi ( j, q ) ) * ds ( q );
                        }
                    }
                }
            }
#    endif // WEAK_BOUNDARY_CONDITIONS
        }

        void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
        {
            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
            const int num_q = num_quadrature_points ( );
            const int num_local_dofs = num_dofs ( 0 );

            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = std::abs ( detJ ( q ) );

                for ( int i = 0; i < num_local_dofs; ++i )
                {
                    lv[dof_index ( i, 0 )] +=
                            wq * ( my_f ( x ( q ) ) * phi ( i, q ) ) * dJ;
                }
            }
        }
        //VECTOR BOUNDARY CONDITIONS

        void operator() ( const Element<double>& element, int facet_number, const Quadrature<double>& quadrature, LocalVector& lv )
        {
#    ifdef WEAK_BOUNDARY_CONDITIONS
            AssemblyAssistant<DIM, double>::initialize_for_facet ( element, quadrature, facet_number );
            if ( !dirichlet_dofs_.empty ( ) )
            {
                const int num_q = num_quadrature_points ( );
                const int num_local_dofs = num_dofs ( 0 );

                for ( int q = 0; q < num_q; ++q )
                {
                    const double wq = w ( q );
                    const double dJ = std::abs ( detJ ( q ) );
                    const double alpha = nitsche_regularization ( q );

                    for ( int i = 0; i < num_local_dofs; ++i )
                    {
                        lv[dof_index ( i, 0 )] += wq * dirichlet_value ( element, 0, q ) * ( -dot ( n ( q ), grad_phi ( i, q ) ) + alpha * phi ( i, q ) ) * ds ( q );
                    }
                }
            }
#    endif // WEAK_BOUNDARY_CONDITIONS
        }

        // driving force

        double my_f ( Vec<DIM, double> pt ) const
        {
            ExactSol sol;
            const double pi = M_PI;
#    if DIMENSION == 2
            return 8. * pi * pi * sol ( pt );
#    else
            return 12. * pi * pi * sol ( pt );
#    endif
        }

        // Neumann boundary condition

        double my_g ( Vec<DIM, double> pt ) const
        {
            return 0.0;
        }

        // Dirichlet boundary condition
#    ifdef WEAK_BOUNDARY_CONDITIONS

        double dirichlet_value ( const Element<double>& element, int var, int q ) const
        {
            assert ( element.is_boundary ( ) );

            std::vector<int> indices;
            // gets the global dof ids
            element.get_dof_indices ( var, indices );

            LOG_DEBUG ( 3, "indices == " << string_from_range ( indices.begin ( ), indices.end ( ) ) );
            LOG_DEBUG ( 3, "dirichlet_dofs == " << string_from_range ( dirichlet_dofs_.begin ( ), dirichlet_dofs_.end ( ) ) );
            LOG_DEBUG ( 3, "dirichlet_values == " << string_from_range ( dirichlet_values_.begin ( ), dirichlet_values_.end ( ) ) );

            std::vector<int>::const_iterator it;
            for ( std::vector<int>::const_iterator ind_it = indices.begin ( ); ind_it != indices.end ( ); ++ind_it )
            {
                it = std::find ( dirichlet_dofs_.begin ( ), dirichlet_dofs_.end ( ), *ind_it );
                // find only the first appearance
                if ( it != dirichlet_dofs_.end ( ) )
                {
                    break;
                }
            }
            assert ( it != dirichlet_dofs_.end ( ) );

            int index = int(it - dirichlet_dofs_.begin ( ) );
            LOG_DEBUG ( 3, "Dirichlet value at dof " << *it << " is " << dirichlet_values_[index] );
            return dirichlet_values_[index];
        }
#    endif // WEAK_BOUNDARY_CONDITIONS

        const std::vector<int> dirichlet_dofs_;
        const std::vector<double> dirichlet_values_;
    };

    template<class ExactSol>
    class L2ErrorIntegrator : private AssemblyAssistant<LaplaceApp::DIM, double>
    {
      public:

        L2ErrorIntegrator ( const LAD::VectorType& dof )
        : dof_ ( dof )
        {
        }

        void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
        {
            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
            evaluate_fe_function ( dof_, 0, approx_sol_ );
        }

        void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, double& value )
        {
            initialize_for_element ( element, quadrature );
            const int num_q = num_quadrature_points ( );

            for ( int q = 0; q < num_q; ++q )
            {
                // error
                const double delta = sol_ ( x ( q ) ) - approx_sol_[q];

                // multiply with weight and transformation factor
                value += w ( q ) * ( delta * delta ) * std::abs ( detJ ( q ) );
            }
        }

      private:

        const LAD::VectorType& dof_;
        ExactSol sol_;
        FunctionValues< double > approx_sol_;

    };

    template<class ExactSol>
    class H1ErrorIntegrator : private AssemblyAssistant<LaplaceApp::DIM, double>
    {
      public:

        H1ErrorIntegrator ( const LAD::VectorType& dof )
        : dof_ ( dof )
        {
        }

        void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
        {
            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
            evaluate_fe_function_gradients ( dof_, 0, approx_grad_u_ );
        }

        void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, double& value )
        {
            initialize_for_element ( element, quadrature );
            const int num_q = num_quadrature_points ( );

            for ( int q = 0; q < num_q; ++q )
            {
                // gradient of exact solution
                Vec<DIM, double> grad_u = sol_.eval_grad ( x ( q ) );

                value += w ( q ) * (
                        dot ( grad_u, grad_u )
                        - 2. * dot ( grad_u, approx_grad_u_[q] )
                        + dot ( approx_grad_u_[q], approx_grad_u_[q] ) )
                        * std::abs ( detJ ( q ) );
            }
        }

      private:

        const LAD::VectorType& dof_;
        ExactSol sol_;
        FunctionValues< Vec<LaplaceApp::DIM, double> > approx_grad_u_;

    };

  private:
    void prepare_linear_algebra ( );

    MPI_Comm comm_;
    mesh::MeshPtr mesh_, master_mesh_;
    doffem::DofPartition<double> dof_partition_;
    VectorSpace<double> space_;
    LAD::MatrixType matrix_;
    LAD::VectorType x_, rhs_;
    std::vector<double> L2_err_, H1_err_;
    CG<LAD> solver_;
    //    GMRES<LAD> solver_;
    SYSTEM la_sys_;
    Couplings<double> couplings_;
    mutable int refinement_level_;
    std::vector<int> fe_degree_;
    int rank_, num_partitions_;

    std::vector<int> dirichlet_dofs_;
    std::vector<LAD::DataType> dirichlet_values_;

};

#endif /* _YALAPLACE_H_ */
