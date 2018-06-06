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
#include <mpi.h>
#include "hiflow.h"

// All names are imported for simplicity.
using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

// Shorten some datatypes with typedefs.
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

typedef std::vector<double> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 2;

// Choose mesh implementation
#define USE_MESH_P4EST

#ifndef WITH_P4EST
#    undef USE_MESH_P4EST
#endif

const int DEBUG_LEVEL = 2;

struct QuadratureSelection
{
    /// Constructor
    /// \param[in] order Desired order of quadrature rule

    QuadratureSelection ( const int order ) : order_ ( order )
    {
    }

    /// Operator to obtain quadrature rule on given element
    /// \param[in] elem Element on which quadrature should be done
    /// \param[out] quadrature Quadrature rule on given Element elem

    void operator() ( const Element<double>& elem,
            Quadrature<double>& quadrature )
    {
        // Get ID of FE type
        const FEType<double>::FiniteElement fe_id =
                elem.get_fe_type ( 0 )->get_my_id ( );

        // Switch by FE type for quadrature selection
        switch ( fe_id )
        {
            case FEType<double>::LAGRANGE_TRI:
            {
                quadrature.set_cell_type ( 2 );
                quadrature.set_quadrature_by_order
                        ( "GaussTriangle", order_ );
                break;
            }
            case FEType<double>::LAGRANGE_QUAD:
            {
                quadrature.set_cell_type ( 3 );
                quadrature.set_quadrature_by_order
                        ( "EconomicalGaussQuadrilateral", order_ );
                break;
            }
            case FEType<double>::LAGRANGE_TET:
            {
                quadrature.set_cell_type ( 4 );
                quadrature.set_quadrature_by_order
                        ( "EconomicalGaussTetrahedron", order_ );
                break;
            }
            case FEType<double>::LAGRANGE_HEX:
            {
                quadrature.set_cell_type ( 5 );
                quadrature.set_quadrature_by_order
                        ( "GaussHexahedron", order_ );
                break;
            }
            default:
            {
                assert ( false );
            }
        };
    }
    // Order of quadrature
    const int order_;
};

// Jump Term Assembler for the jumps over the inner edges

class JumpTermAssembler : private DGAssemblyAssistant <DIMENSION, double>
{
  public:

    JumpTermAssembler ( CoupledVector<Scalar>& u_h )
    : u_h_ ( u_h )
    {
    }

    void operator() ( const Element<double>& left_elem,
            const Element<double>& right_elem,
            const Quadrature<double>& left_quad,
            const Quadrature<double>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            double& local_val )
    {
        const bool is_boundary =
                (
                right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY
                );

        if ( is_boundary
             || ( left_if_side == right_if_side ) )
        {
            local_val = 0.;
            return;
        }

        this->initialize_for_interface ( left_elem, right_elem,
                                         left_quad, right_quad,
                                         left_facet_number, right_facet_number,
                                         left_if_side, right_if_side );

        left_grad_u_h.clear ( );
        right_grad_u_h.clear ( );

        this->trial ( ).evaluate_fe_function_gradients ( u_h_, 0, left_grad_u_h );
        this->test ( ).evaluate_fe_function_gradients ( u_h_, 0, right_grad_u_h );
        const int num_q = num_quadrature_points ( );

        double h_E = 0.;
        std::vector<double> y_coord;
        for ( int i = 0; i < num_q; ++i )
        {
            for ( int j = 0; j < num_q; ++j )
            {
                h_E = std::max ( h_E, norm ( this->x ( i ) - this->x ( j ) ) );
            }
            y_coord.push_back ( this->x ( i )[1] );
        }

        // Loop over quadrature points on each edge
        int rank;
        MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

        LOG_DEBUG ( 2, "[" << rank << "] trial cell = " << left_elem.get_cell_index ( ) << " test cell = " << right_elem.get_cell_index ( )
                    << " trial normal = " << this->trial ( ).n ( 0 )[0] << " , " << this->trial ( ).n ( 0 )[1]
                    << " test normal = " << this->test ( ).n ( 0 )[0] << " , " << this->test ( ).n ( 0 )[1]
                    << " trial grad = " << left_grad_u_h[0][0] << " , " << left_grad_u_h[0][1]
                    << " test grad = " << right_grad_u_h[0][0] << " , " << right_grad_u_h[0][1]
                    //<< " if lenght = " << h_E
                    << " " << string_from_range ( y_coord.begin ( ), y_coord.end ( ) ) );

        for ( int q = 0.; q < num_q; ++q )
        {
            const double dS = std::abs ( this->ds ( q ) );
            local_val += this->w ( q )

                    * std::abs (
                                 dot ( this->trial ( ).n ( q ), left_grad_u_h[q] )
                                 + dot ( this->test ( ).n ( q ), right_grad_u_h[q] )
                                 )
                    * std::abs (
                                 dot ( this->trial ( ).n ( q ), left_grad_u_h[q] )
                                 + dot ( this->test ( ).n ( q ), right_grad_u_h[q] )
                                 )
                    * dS;
        }
    }
    const CoupledVector<Scalar>& u_h_;
    FunctionValues< Vec<DIMENSION, double> > left_grad_u_h;
    FunctionValues< Vec<DIMENSION, double> > right_grad_u_h;
};

