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

/// \author Staffan Ronnas, Simon Gawlok

#ifndef HP_ASSEMBLY_TEST_H
#    define HP_ASSEMBLY_TEST_H

#    include "hiflow.h"

#    include <cmath>
#    include <iostream>
#    include <string>
#    include <vector>
#    include <mpi.h>

using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType Vector;
typedef LAD::MatrixType Matrix;

static const char* DATADIR = MESH_DATADIR;
static const double DEG = 2.;

class LaplaceApp
{
  public:
    static const int DIM = 2;
    static const int MAX_REFINEMENT_LEVEL = 5;
    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;
    static const int M = 1;
    static const int N = 1;
    static const int P = 1;

    LaplaceApp ( );
    ~LaplaceApp ( );

    void read_mesh ( const std::string& filename );
    void prepare_system ( );
    void assemble_system ( );
    void solve_system ( );
    double compute_error ( );
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
#    if 1

    struct ExactSol
    {

        double operator() ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;

            if ( LaplaceApp::DIM == 2 )
            {
                const double x = pt[0];
                const double y = pt[1];
                return std::sin ( M * pi * x ) * std::sin ( N * pi * y );
            }
            else if ( LaplaceApp::DIM == 3 )
            {
                const double x = pt[0];
                const double y = pt[1];
                const double z = pt[2];
                return std::sin ( M * pi * x ) * std::sin ( N * pi * y ) * std::sin ( P * pi * z );
            }
            else
            {
                assert ( false );
            }

        }

        Vec<DIM, double> eval_grad ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;

            if ( DIM == 2 )
            {
                Vec<DIM, double> grad;
                const double x = pt[0];
                const double y = pt[1];
                grad[0] = M * pi * std::cos ( M * pi * x ) * std::sin ( N * pi * y );
                grad[1] = N * pi * std::sin ( M * pi * x ) * std::cos ( N * pi * y );
                return grad;
            }
            else
            {
                assert ( false );
            }
        }
    };
#    endif

#    if 0

    struct ExactSol
    {

        double operator() ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;

            if ( LaplaceApp::DIM == 2 )
            {
                const double x = pt[0];
                const double y = pt[1];
                return std::cos ( M * pi * x ) * std::cos ( N * pi * y );
            }
            else if ( LaplaceApp::DIM == 3 )
            {
                const double x = pt[0];
                const double y = pt[1];
                const double z = pt[2];
                return std::cos ( M * pi * x ) * std::cos ( N * pi * y ) * std::cos ( P * pi * z );
            }
            else
            {
                assert ( false );
            }

        }

        Vec<DIM, double> eval_grad ( const Vec<DIM, double>& pt ) const
        {
            const double pi = M_PI;

            if ( DIM == 2 )
            {
                Vec<DIM, double> grad;
                const double x = pt[0];
                const double y = pt[1];
                grad[0] = -M * pi * std::sin ( M * pi * x ) * std::cos ( N * pi * y );
                grad[1] = -N * pi * std::cos ( M * pi * x ) * std::sin ( N * pi * y );
                return grad;
            }
            else
            {
                assert ( false );
            }
        }
    };
#    endif

#    if 0

    struct ExactSol
    {

        double operator() ( const Vec<DIM, double>& pt ) const
        {
            const double x = pt[0];
            const double y = pt[1];
            return std::pow ( x, DEG ) * std::pow ( y, DEG );
        }

        Vec<DIM, double> eval_grad ( const Vec<DIM, double>& pt ) const
        {
            const double x = pt[0];
            const double y = pt[1];
            Vec<DIM> grad;
            grad[0] = ( DEG ) * std::pow ( x, DEG - 1 ) * std::pow ( y, DEG );
            grad[1] = ( DEG ) * std::pow ( x, DEG ) * std::pow ( y, DEG - 1 );
            return grad;
        }
    };
#    endif

    class LocalLaplaceAssembler : private AssemblyAssistant<LaplaceApp::DIM, double>
    {
      public:

        void operator() ( const Element<double>& element,
                const Quadrature<double>& quadrature, GlobalAssembler<double>::LocalMatrix& lm )
        {

            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

            const int num_q = num_quadrature_points ( );
            const int num_local_dofs = num_dofs ( 0 );

            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = detJ ( q );

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

        void operator() ( const Element<double>& element,
                const Quadrature<double>& quadrature, GlobalAssembler<double>::LocalVector& lv )
        {

            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

            const int num_q = num_quadrature_points ( );
            const int num_local_dofs = num_dofs ( 0 );

            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = detJ ( q );

                for ( int i = 0; i < num_local_dofs; ++i )
                {
                    lv[dof_index ( i, 0 )] +=
                            wq * ( my_f ( x ( q ) ) * phi ( i, q ) ) * dJ;
                }
            }
        }
#    if 1

        double my_f ( Vec<DIM, double> pt ) const
        {
            ExactSol sol;
            const double pi = M_PI;
            if ( DIM == 2 )
            {
                return (M * M + N * N ) * pi * pi * sol ( pt );
            }
            else if ( DIM == 3 )
            {
                return (M * M + N * N + P * P ) * pi * pi * sol ( pt );
            }
        }
#    endif
#    if 0

        double my_f ( Vec<DIM, double> pt ) const
        {
            const double x = pt[0];
            const double y = pt[1];
            return -2. * ( x * x + y * y );
            //            return -DEG * (DEG - 1.) *
            //                (std::pow(x, DEG - 2.) * std::pow(y, DEG) + std::pow(x, DEG) * std::pow(y, DEG - 2.));
        }
#    endif
    };

    template<class ExactSol>
    class L2ErrorIntegrator : private AssemblyAssistant<LaplaceApp::DIM, double>
    {
      public:

        L2ErrorIntegrator ( const LAD::VectorType& dof )
        : dof_ ( dof )
        {
        }

        void operator() ( const Element<double>& element,
                const Quadrature<double>& quadrature, double& value )
        {

            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
            evaluate_fe_function ( dof_, 0, approx_sol_ );

            const int num_q = num_quadrature_points ( );

            for ( int q = 0; q < num_q; ++q )
            {
                // error
                double delta = sol_ ( x ( q ) ) - approx_sol_[q];

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

        void operator() ( const Element<double>& element,
                const Quadrature<double>& quadrature, double& value )
        {
            AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
            evaluate_fe_function_gradients ( dof_, 0, approx_grad_u_ );

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
    void prepare_bc ( );
    void prepare_local_refinements ( );

    MPI_Comm comm_;

    mesh::MeshPtr mesh_, old_mesh_;
    doffem::DofPartition<double> dof_partition_;
    VectorSpace<double> space_;
    LAD::MatrixType matrix_;
    LAD::VectorType x_, rhs_;
    std::vector<double> err_;
    CG<LAD> solver_;
    //GMRES<LAD> solver_;
    SYSTEM la_sys_;
    Couplings<double> couplings_;
    mutable int refinement_level_;
    std::vector<int> fe_degree_;
    Visualization<double> visu_;
    int rank_, num_partitions_;
    std::vector<int> dirichlet_dofs_;
    std::vector<double> dirichlet_values_;
    HpFemAssembler<double> global_asm_;
};

#endif /* _YALAPLACE_H_ */
