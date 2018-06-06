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

/// \author Philipp Gerstner

// System includes.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "hiflow.h"

// Use matrix free implementation of eigensolver
#define nMATRIXFREE
// Use mixed matrix / matrix-free) implementation for eigensolver in case of generalized eigenvalue problem Ax= lambda Bx: explicit matrix for B, matrix-free for A
#define MATMIXED

// Make problem unsymmetric
#define ADD_CONVECTION

// Make problem complex
#define COMPLEX_OPERATOR

// Select Linear Algebra Platform
#define nUSE_PETSC
#define USE_HYPRE

#ifndef WITH_PETSC
#    undef USE_PETSC
#endif

#ifndef WITH_HYPRE
#    undef USE_HYPRE
#endif

// All names are imported for simplicity.
using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

// Shorten some datatypes with typedefs.
#ifdef USE_PETSC
typedef LADescriptorPETSc LAD;
#else
#    ifdef USE_HYPRE
typedef LADescriptorHypreD LAD;
#    else
typedef LADescriptorCoupledD LAD;
#    endif
#endif

typedef LAD::DataType Scalar;
typedef LAD::VectorType VectorType;
typedef LAD::MatrixType MatrixType;

typedef std::vector<double> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 2;

// Functor used to impose u = 0 on the boundary.

struct DirichletZero
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {

        // return array with Dirichlet values for dof:s on boundary face
        return std::vector<double>( coords_on_face.size ( ), 0. );
    }
};

// Functor used for the local assembly of the stiffness matrix and load vector.

class LocalPoissonAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    LocalPoissonAssembler ( )
    {
        this->mode_ = 1;
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        // compute local matrix
        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lm.Clear ( );
        lm.Resize ( total_dofs, total_dofs );
        const double u_x = 1.1;
        const double u_y = 0.2;

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                {
                    if ( this->mode_ == 1 )
                    {
                        lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * dot ( grad_phi ( j, q ), grad_phi ( i, q ) ) * std::abs ( detJ ( q ) );
#ifdef ADD_CONVECTION
                        lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * ( u_x * grad_phi ( j, q )[0] + u_y * grad_phi ( j, q )[1] ) * phi ( i, q ) * std::abs ( detJ ( q ) );
#endif
                    }
                    else
                    {
                        //				lm(dof_index(i, 0), dof_index(j, 0)) += wq * 0.5 * dot(grad_phi(j, q), grad_phi(i, q)) * std::abs(detJ(q));
                        lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * phi ( j, q ) * phi ( i, q ) * std::abs ( detJ ( q ) );
                    }
                }
            }
        }
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lv.clear ( );
        lv.resize ( total_dofs, 0. );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                lv[dof_index ( i, 0 )] += wq * f ( x ( q ) ) * phi ( i, q ) * std::abs ( detJ ( q ) );
            }
        }
    }

    void set_mode_to_real ( )
    {
        this->mode_ = 1;
    }

    void set_mode_to_imag ( )
    {
        this->mode_ = -1;
    }

    double f ( hiflow::Vec<DIMENSION, double> pt )
    {
        return 1.0;
    }

  private:
    int mode_;

};

class LocalMassAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        // compute local matrix
        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lm.Clear ( );
        lm.Resize ( total_dofs, total_dofs );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * phi ( j, q ) * phi ( i, q ) * std::abs ( detJ ( q ) );
                }
            }
        }
    }
};
