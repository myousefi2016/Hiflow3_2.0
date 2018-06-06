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

#ifndef HIFLOW_ADAPTIVITY_FE_INTERPOLATION
#    define HIFLOW_ADAPTIVITY_FE_INTERPOLATION

/// \author Philipp Gerstner

#    include <map>
#    include <string>
#    include <vector>

#    include <boost/function.hpp>
#    include <mpi.h>
#    include "space/vector_space.h"
#    include "mesh/mesh.h"
#    include "mesh/types.h"

#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/couplings.h"

using namespace hiflow::mesh;
using namespace hiflow::la;

// TODO: check whether is_cg = false works
namespace hiflow
{
    ///
    /// \class FEInterpolation fe_interpolation.h
    /// \brief abstract base class for interpolating FE functions from one space to another one 

    template<class LAD, int DIM>
    class FEInterpolation
    {
        typedef typename LAD::DataType DataType;
        typedef typename LAD::VectorType VectorType;
        typedef std::vector<DataType> Coord;

      public:
        FEInterpolation ( );

        ~FEInterpolation ( )
        {
            this->clear ( );
        }

        /// \brief Initialize fe interpolation, this function involves the follwoing steps: <br> 
        /// check whether input space is suitable for patching <br>
        /// build dof interpoaltion map from interpolating space to input space
        /// @param[in] from_space space to be interpolated from
        /// @param[in] to_space space to be interpolated to
        virtual void init ( VectorSpace<DataType>* from_space, VectorSpace<DataType>* to_space );

        /// \brief apply fe interpolation 
        /// @param[in] from_vector vector to be interpolated
        /// @param[out] to_vector interpolated vector

        virtual void interpolate ( const VectorType& from_vector, VectorType& to_vector ) const
        {
        };

        /// \brief set tolerance for identifying dofs by physical coordinates
        /// @param[in] eps tolerance

        virtual void set_dof_tolerance ( DataType eps )
        {
            this->eps_ = eps;
        }

        /// \brief get space to be interpolated from
        /// @return reference to space

        virtual VectorSpace<DataType>* get_from_space ( )
        {
            return this->from_space_;
        }

        /// \brief get space to be interpolated to
        /// @return reference to space

        virtual VectorSpace<DataType>* get_to_space ( )
        {
            return this->to_space_;
        }

        /// \brief clear all data structs
        virtual void clear ( );

        virtual bool is_initialized ( )
        {
            return this->initialized_;
        }

      protected:
        /// \brief checks if input space is suitable for patch interpolation

        virtual bool check_space ( ) const
        {
        };

        /// pointers to input and output space
        VectorSpace<DataType>* from_space_;
        VectorSpace<DataType>* to_space_;

        /// pointers to input and output meshes
        MeshPtr from_mesh_;
        MeshPtr to_mesh_;

        /// FE degrees of interpolating space
        std::vector<int> degrees_;

        /// CG / DG flags for interpolating space
        std::vector<bool> is_cg_;

        /// number of variables in input space
        int num_var_;

        /// tolerance for dof identification
        DataType eps_;

        /// flag indicating whether patchinterpolation is initialized for given space
        bool initialized_;

        /// MPI rank
        int rank_;

        /// topological dimension of input space
        int tdim_;

        /// geometrical dimension of input space
        int gdim_;
    };

    ///
    /// \class FEInterpolationLagrange fe_interpolation.h
    /// \brief class for interpolating Lagrange FE functions from one space to another one \br
    /// interpolation works in the following way: given u_f (x) = \sum_i a_i \phi_i(x) \in V_f (from_space) \br
    /// compute u_t (x) = \sum_j b_j \psi_j(x) \in V_t (to_space) such that \bre
    /// u_t(x_j) = u_f(x_j) for all dof positions x_j (t-> to_space) with \psi_k(x_j) = delta_{kj}

    template<class LAD, int DIM>
    class FEInterpolationLagrange : virtual public FEInterpolation<LAD, DIM>
    {
        typedef typename LAD::DataType DataType;
        typedef typename LAD::VectorType VectorType;
        typedef std::vector<DataType> Coord;

      public:
        FEInterpolationLagrange ( );

        ~FEInterpolationLagrange ( )
        {
            this->clear ( );
        }

        /// \brief Initialize fe interpolation, this function involves the follwoing steps: <br> 
        /// check whether input space is suitable for patching <br>
        /// build dof interpoaltion map from interpolating space to input space
        /// @param[in] from_space space to be interpolated from
        /// @param[in] to_space space to be interpolated to
        virtual void init ( VectorSpace<DataType>* from_space, VectorSpace<DataType>* to_space );

        /// \brief apply fe interpolation 
        /// @param[in] from_vector vector to be interpolated
        /// @param[out] to_vector interpolated vector
        virtual void interpolate ( const VectorType& from_vector, VectorType& to_vector ) const;

        /// \brief clear all data structs
        virtual void clear ( );

      protected:

        /// \brief checks if input space is suitable for patch interpolation
        virtual bool check_space ( ) const;

        /// \brief build dof interpolation mapping
        virtual void create_dof_mapping ( );

        /// \brief dof_weights(j) = \{ (i_1, w_1), ... (i_n(j), w_n(j) \} such that \br
        /// b_j = \sum_{l} a_dof_weight(j).i_l * dof_weights(j).w_l  
        std::vector< std::map<int, std::map<int, double> > > dof_weights_;
    };

}
#endif
