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

#ifndef HIFLOW_ADAPTIVITY_PATCH_INTERPOLATION
#    define HIFLOW_ADAPTIVITY_PATCH_INTERPOLATION

/// \author Philipp Gerstner

#    include <map>
#    include <string>
#    include <vector>

#    include <boost/function.hpp>
#    include <mpi.h>
#    include "space/vector_space.h"
#    include "mesh/mesh.h"
#    include "mesh/types.h"
#    include "mesh/mesh_pXest.h"
#    include "mesh/entity.h"
#    include "assembly/assembly.h"

#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/couplings.h"

using namespace hiflow::mesh;
using namespace hiflow::la;

// TODO: check whether is_cg = false works
namespace hiflow
{
    ///
    /// \class SpacePatchInterpolation patch_interpolation.h
    /// \brief class for interpolating FE functions from one space to another one 
    /// with a coarser mesh and higher polynomial degree
    /// 

    template<class LAD, IMPL MESH_IMPL>
    class SpacePatchInterpolation
    {
        typedef typename LAD::DataType DataType;
        typedef typename LAD::VectorType VectorType;
        typedef std::vector<DataType> Coord;

      public:
        SpacePatchInterpolation ( );

        ~SpacePatchInterpolation ( )
        {
            this->clear ( );
        }

        /// \brief Initialize patch interpolation, this function involves the follwoing steps: <br> 
        /// check whether input space is suitable for patching <br>
        /// copy mesh of input space and uniformly coarsen it <br>
        /// setup interpolating space <br>
        /// build dof interpoaltion map from interpolating space to input space
        /// @param[in] input_space space to be patched
        virtual void init ( VectorSpace<DataType>const * input_space );

        /// \brief apply patch interpolation 
        /// @param[in] input_vector vector to be interpoalted
        /// @param[out] vector interpolated vector
        virtual void interpolate ( const VectorType& input_vector, VectorType& vector ) const;

        /// \brief set tolerance for identifying dofs by physical coordinates
        /// @param[in] eps tolerance

        virtual void set_dof_tolerance ( DataType eps )
        {
            this->eps_ = eps;
        }

        /// \brief get patch interpolating vector space
        /// @return reference to space

        virtual VectorSpace<DataType>* get_space ( )
        {
            return &this->space_;
        }

        /// \brief get dof interpoaltion mapping
        /// @return reference to mapping 

        virtual std::map<int, int>& get_mapping ( )
        {
            return this->high2low_;
        }

        /// \brief get linear algebra couplings, used for seeting up interpoalting vector by user
        /// @return reference to couplings

        virtual const Couplings<DataType>& get_couplings ( ) const
        {
            return this->couplings_;
        }

        /// \brief clear all data structs
        virtual void clear ( );

      protected:

        /// \brief checks if input space is suitable for patch interpolation
        virtual bool check_space ( ) const;

        /// \brief build dof interpolation mapping
        virtual void create_dof_mapping ( );

        /// \brief copy input mesh
        virtual void copy_mesh ( );

        /// \brief compute ghost cells of interpolating mesh
        virtual void compute_ghost ( );

        /// pointers to input space and input mesh
        VectorSpace<DataType> const * input_space_;
        MeshPtr input_mesh_;

        /// global assembler used for computing couplings
        StandardGlobalAssembler<DataType> global_asm_;

        /// interpolating space
        VectorSpace<DataType> space_;

        /// linear algebra couplings for interpolating space
        Couplings<DataType> couplings_;

        /// patched mesh
        MeshPtr mesh_;

        /// FE degrees of interpolating space
        std::vector<int> degrees_;

        /// CG / DG flags for interpolating space
        std::vector<bool> is_cg_;

        /// number of variables in input space
        int num_var_;

        /// dof mapping from interpolating space to input space
        std::map<int, int> high2low_;

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

    // TODO generalize for arbitrary degree
    /// \class TimePatchInterpolation patch_interpolation.h
    /// \brief class for evaluating Time-FE functions with degree <= 1 on time interval patch [t_{n-1}, t_n] \cup [t_n, t_{n+1}]
    ///

    template<typename T, class DataType>
    class TimePatchInterpolation
    {
      public:

        TimePatchInterpolation ( )
        {
        }

        ~TimePatchInterpolation ( )
        {
        }

        /// \brief Evaluate piecewise constant function
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1} (NOT USED)
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1}
        virtual T constant ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Evaluate piecewise linear function
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1}
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1}
        virtual T linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Evaluate temporal derivative of piecewise linear function
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1}
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1}
        virtual T dt_linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Evaluate linear function on complete patch [t_{n-1}, t_n] \cup [t_n, t_{n+1}]
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1} (NOT USED)
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1}    
        virtual T patch_linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Evaluate temporal derivative of  linear function on complete patch [t_{n-1}, t_n] \cup [t_n, t_{n+1}]
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1}
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1} 
        virtual T patch_dt_linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Evaluate quadrative function on complete patch [t_{n-1}, t_n] \cup [t_n, t_{n+1}]
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1}
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1}  
        virtual T patch_quadratic ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Evaluate temporal derivative of quadrative function on complete patch [t_{n-1}, t_n] \cup [t_n, t_{n+1}]
        /// @param[in] c relative time in active interval, c \in [0,1]
        /// @param[in] val_prev time dof value at time t_{n-1}
        /// @param[in] val time dof value at time t_{n}
        /// @param[in] val_next time dof value at time t_{n+1} 
        virtual T patch_dt_quadratic ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const;

        /// \brief Set time step size and active interval
        /// @param[in] dT_pc (t_{n} - t_{n-1})
        /// @param[in] dT_nc (t_{n+1} - t_{n})
        /// @param[in] rel_time 0: active interval is [t_{n-1}, t_n], 1: active interval is [t_n, t_{n+1}]
        virtual void set_time_steps ( DataType dT_pc, DataType dT_cn, int rel_time );

        /// \brief compute coefficients for evaluating functions
        virtual void compute_weight_tau_coeff ( );

        /// \brief clear allocated data
        virtual void clear ( );

      protected:
        /// \brief get absolute time in patch interval
        /// @param[in c relatvie time \in [0,1]
        /// @return absolute time
        DataType get_absolut_time ( DataType c ) const;

        DataType dT_pc_;
        DataType dT_cn_;
        int rel_time_;

        // coefficients for patchwise quadratic time interpolation for continuous, piecewise linear functions 
        DataType a_n2_;
        DataType a_n1_;
        DataType a_c2_;
        DataType a_c1_;
        DataType a_p2_;
        DataType a_p1_;
        DataType a_p0_;

        // coefficients for patchwise linear time interpolation for disconinuous, piecewsie constant functions
        DataType b_n1_;
        DataType b_n0_;
        DataType b_c1_;
        DataType b_c0_;
    };
}
#endif
