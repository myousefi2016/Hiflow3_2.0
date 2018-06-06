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

#ifndef HIFLOW_ADAPTIVITY_DYNAMIC_MESH
#    define HIFLOW_ADAPTIVITY_DYNAMIC_MESH

/// \author Philipp Gerstner

#    include <map>
#    include <string>
#    include <vector>

#    include <boost/function.hpp>
#    include <mpi.h>
#    include "space/vector_space.h"
#    include "mesh/mesh.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/couplings.h"
#    include "adaptivity/time_mesh.h"
#    include "adaptivity/fe_interpolation.h"

using namespace hiflow::mesh;
using namespace hiflow::la;

namespace hiflow
{
    template<class LAD, IMPL MESH_IMPL, int DIM>
    class DynamicMeshProblem;

    ///
    /// \class DynamicMeshHandler dynamic_mesh.h
    /// \brief Class that manages mesh transfer operations in a time stepping loop
    ///
    ///

    template<class LAD, IMPL MESH_IMPL, int DIM>
    class DynamicMeshHandler
    {
        typedef typename LAD::DataType DataType;
        typedef typename LAD::VectorType VectorType;

      public:
        DynamicMeshHandler ( DynamicMeshProblem<LAD, MESH_IMPL, DIM>* problem, MPI_Comm comm );

        ~DynamicMeshHandler ( )
        {
            this->clear ( );
        }

        virtual void clear ( );

        ////////////////////////////////////////////////    
        ///////////// Use in time loop /////////////////
        /// \brief update mesh, vector spaces and LA in case of mesh change
        /// @param[in] time_step  index of considered time step
        /// @param[in] mode 1: primal, -1: dual, 0: error estimation
        /// @param[in] initial_call indicates whether functions is called for the first time 
        /// @return true if mesh was updated
        virtual bool update ( int time_step, int mode, bool initial_call );

        /// \brief refine mesh specified by index
        /// @param[in] mesh_index index of considered mesh
        /// @param[in] markers flags for mesh refinement
        virtual void refine ( int mesh_index, std::vector<int>& markers );

        /// \brief initialize fe spaces for all meshes in mesh_list
        virtual void init_fe_spaces ( );

        /// \brief interpolate vector between different fe spaces
        /// @param[in] in_vec input vector for interpolation
        /// @param[in] in_index index of input space
        /// @param[in] out_vec output vector for interpolation
        /// @param[in] out_index index of output space
        virtual void interpolate_vector ( const VectorType& in_vec, int in_index, VectorType& out_vec, int out_index );

        /// \brief load specific mesh
        /// @param[in] adapt_counter iteration index of outer adaption loop
        /// @param[in] num_mesh number of meshes to load from file
        /// @param[in] prefix path to saved mesh files
        virtual void load_mesh_from_file ( int adapt_counter, int num_mesh, std::string& prefix );

        /// \brief save specific mesh
        /// @param[in] adapt_counter iteration index of outer adaption loop
        /// @param[in] prefix path to save mesh files
        virtual void save_mesh_to_file ( int adapt_counter, std::string& prefix );

        /// \brief visualize specific mesh
        /// @param[in] adapt_counter iteration index of outer adaption loop
        /// @param[in] prefix path to save visualization files
        virtual void visualize_mesh ( int adapt_counter, std::string& prefix );

        /// \brief add points when to change the mesh
        /// @param[in] additional changes time points to add for mesh change
        /// @param[in,out] indicator_mesh_indices map: indicator-time-step to mesh index 
        virtual void add_mesh_change_times ( const std::vector<DataType>& additional_changes, std::vector<int>& indicator_mesh_indices );

        ///////////////////////////////////////////////////
        //////////// setup ////////////////////////////////
        /// \brief set point to time mesh
        /// @param[in] mesh pointer to time mesh

        virtual void set_time_mesh ( TimeMesh<DataType>* mesh )
        {
            this->time_mesh_ = mesh;
        }

        /// \brief set initial time points for mesh change
        /// @param[in] times initial mesh change times

        virtual void set_initial_mesh_change_times ( const std::vector<DataType>& times )
        {
            assert ( this->mesh_change_times_.size ( ) == 0 );
            for ( int l = 0; l < times.size ( ); ++l )
            {
                this->mesh_change_times_.insert ( times[l] );
            }
        }

        /// \brief set base mesh
        /// @param[in] mesh pointer to space mesh

        virtual void set_initial_mesh ( MeshPtr mesh )
        {
            assert ( mesh != 0 );
            this->initial_mesh_ = mesh;
            this->tdim_ = mesh->tdim ( );
            this->gdim_ = mesh->gdim ( );

            this->fill_mesh_list ( );
        }

        /// \brief set type of mesh change adaption mode
        /// @param[in] type "None": single mesh, "Fixed": mesh changes according to initial mesh change times, \br
        /// "FixedFraction": change mesh at those (fixed number of ) times steps with largest difference in subsequent error indicators \br
        /// "FixedError": change mesh if difference of subsequent error indicators is above a certain threshold

        virtual void set_adapt_type ( std::string type )
        {
            this->adapt_type_ = type;
        }

        /// \brief
        /// @param[in] mode 1: primal loop, -1: dual loop, 0: error estimation 
        /// @param[in] type type of vectors, 1: primal vector, -1: dual vectors
        /// @param[in] vectors pointers to upate vectors
        virtual void set_update_vectors ( int mode, int type, std::vector<VectorType*> vectors );

        ////////////////////////////////////////////////////
        ////////////// get functions ///////////////////////
        /// \brief get mesh change points
        /// @return points

        virtual std::vector<DataType> get_mesh_change_times ( ) const
        {
            return this->mesh_change_times_.data ( );
        }

        /// \brief get mesh change points
        /// @return indices of points

        virtual std::vector<int> get_mesh_change_steps ( )
        {
            this->update_mesh_change_steps ( );
            return this->mesh_change_steps_.data ( );
        }

        /// \brief get index of active mesh
        /// @return index

        virtual int get_active_mesh_index ( )
        {
            return this->active_mesh_index_;
        }

        /// \brief get pointer to vector space corresponding to specified mesh
        /// @param[in] time_step  index of considered time step
        /// @param[in] mode 1: primal -1: dual
        /// @return pointer to vector space
        virtual VectorSpace<DataType>* get_space_by_step ( int time_step, int mode ) const;

        /// \brief get pointer to vector space corresponding to specified mesh
        /// @param[in] mesh_index index of considered mesh
        /// @param[in] mode 1: primal -1: dual
        /// @return pointer to vector space
        virtual VectorSpace<DataType>* get_space_by_index ( int mesh_index, int mode ) const;

        /// \brief get pointer to vector space corresponding to active mesh
        /// @param[in] mode 1: primal -1: dual
        /// @return pointer to vector space

        virtual VectorSpace<DataType>* get_active_space ( int mode ) const
        {
            return this->get_space_by_index ( this->active_mesh_index_, mode );
        }

        /// \brief get pointer to mesh for specific time step
        /// @param[in] time_step  index of considered time step
        /// @return mesh pointer
        virtual MeshPtr get_mesh_by_step ( int time_step ) const;

        /// \brief get pointer to mesh for specific index
        /// @param[in] mesh_index index of considered mesh
        /// @return mesh pointer
        virtual MeshPtr get_mesh_by_index ( int mesh_index ) const;

        /// \brief get pointer to active mesh
        /// @return mesh pointer

        virtual MeshPtr get_active_mesh ( ) const
        {
            return this->mesh_;
        }

        /// \brief get underlying MPI communicator
        /// @return communicator

        virtual MPI_Comm& get_mpi_comm ( )
        {
            return this->comm_;
        }

        /// \brief get LA couplings for active mesh
        /// @param[in] mode 1: primal -1: dual
        /// @return couplings

        virtual Couplings<DataType>* get_active_couplings ( int mode ) const
        {
            return this->get_couplings_by_index ( this->active_mesh_index_, mode );
        }

        /// \brief get LA couplings for specified mesh index
        /// @param[in] mesh_index index of considered mesh
        /// @param[in] mode 1: primal -1: dual
        /// @return couplings
        virtual Couplings<DataType>* get_couplings_by_index ( int mesh_index, int mode ) const;

        /// \brief get number of different meshes
        /// @return number

        virtual int num_mesh ( ) const
        {
            return this->mesh_list_.size ( );
        }

        /// \brief get number of mesh change points
        /// @return number

        virtual int num_mesh_change ( ) const
        {
            return this->mesh_change_times_.size ( );
        }

        /// \brief get index of mesh corresponding to specified time point
        /// @param[in] time point in time
        /// @return mesh index
        virtual int mesh_index ( DataType time ) const;

        /// \brief get index of mesh corresponding to specified time point index
        /// @param[in] time_step  index of considered time step
        /// @return mesh index
        virtual int mesh_index ( int time_step ) const;

        /// \brief get index of first time point of specified mesh
        /// @param[in] mesh_index index of considered mesh
        /// @return time step
        virtual int first_step_for_mesh ( int mesh_index ) const;

        /// \brief get index of last time point of specified mesh
        /// @param[in] mesh_index index of considered mesh
        /// @return time step
        virtual int last_step_for_mesh ( int mesh_index ) const;

        /// \brief get number of cells in mesh specified by time step
        /// @param[in] time_step  index of considered time step
        /// @return number of cells
        virtual int num_cells_by_step ( int time_step ) const;

        /// \brief get number of cells in mesh specified by index
        /// @param[in] mesh_index index of considered mesh
        /// @return number of cells
        virtual int num_cells_by_index ( int mesh_index ) const;

      protected:
        /// \brief fill mesh list by copying initial mesh
        virtual void fill_mesh_list ( );

        /// \brief create fe interpolator objects
        /// @param[in] num_mesh number of meshes
        virtual void create_fe_interpolator ( int num_mesh );

        /// \brief update indices of mesh change time points
        virtual void update_mesh_change_steps ( );

        /// \brief checks whether mesh changes 
        /// @param[in] last_time last time point 
        /// @param[in] cur_time current time point
        /// @return true if mesh changes between time points
        virtual bool need_to_change_mesh ( DataType last_time, DataType cur_time ) const;

        /// \brief checks whether mesh changes
        /// @param[in] last_step index of last time point
        /// @param[in] cur_step index of current time point
        /// @return true if mesh changes between time steps
        virtual bool need_to_change_mesh ( int last_step, int cur_step ) const;

        /// \brief checks whether mesh changes
        /// @param[in] time_step  index of considered time step
        /// @return true if mesh changes at time step
        virtual bool need_to_change_mesh ( int time_step ) const;

        TimeMesh<DataType>* time_mesh_;

        int tdim_;

        int gdim_;

        int active_mesh_index_;

        MeshPtr initial_mesh_;

        MeshPtr mesh_;

        int rank_;

        MPI_Comm comm_;

        std::string adapt_type_;

        VectorSpace<DataType>* space_;
        VectorSpace<DataType>* space_dual_;

        std::vector<MeshPtr> mesh_list_;

        SortedArray<DataType> mesh_change_times_;

        SortedArray<int> mesh_change_steps_;

        std::vector<VectorSpace<DataType>* > space_list_;

        std::vector<VectorSpace<DataType>* > space_list_dual_;

        std::vector<VectorSpace<DataType>* > space_list_tmp_;

        std::vector< std::vector< FEInterpolation<LAD, DIM>* > > fe_interpolator_;

        std::vector<Couplings<DataType>* > couplings_list_;
        std::vector<Couplings<DataType>* > couplings_list_dual_;

        std::vector< std::vector< std::vector< bool > > > coupling_vars_list_;
        std::vector< std::vector< std::vector< bool > > > coupling_vars_list_dual_;

        std::vector<VectorType*> vec_primal_primal_;
        std::vector<VectorType*> vec_dual_primal_;
        std::vector<VectorType*> vec_dual_dual_;
        std::vector<VectorType*> vec_est_primal_;
        std::vector<VectorType*> vec_est_dual_;

        DynamicMeshProblem<LAD, MESH_IMPL, DIM>* problem_;
    };

    ///
    /// \class  DynamicMeshProblem dynamic_mesh.h
    /// \brief base class for problems involvoing a mesh that changes during a time loop
    ///
    /// 

    template<class LAD, IMPL MESH_IMPL, int DIM>
    class DynamicMeshProblem
    {
      public:

        DynamicMeshProblem ( MPI_Comm comm )
        : dmh_ ( new DynamicMeshHandler<LAD, MESH_IMPL, DIM>( this, comm ) ),
        active_mesh_index_ ( -1 )
        {
        }

        ~DynamicMeshProblem ( )
        {
        }

        /// \brief routine which is called in init_fe_spaces() of DynamicMeshHandler to setup primal and dual fe space for a given mesh
        /// @param[out] space reference to space to be set up
        /// @param[out] coupling_vars vector containing information of coupled variables
        /// @param[in] mesh pointer to specific mesh
        /// @param[in] mode 1: primal space, -1: dual space 
        virtual void setup_space ( VectorSpace<typename LAD::DataType>& space, std::vector< std::vector<bool> >& coupling_vars, MeshPtr mesh, int mode ) = 0;

        /// \brief routine which is called in udpate() of DynamicMeshHandler to setup required LA objects for specific fe space in case of primal problem time loop
        /// @param[in] space reference to primal space
        /// @param[in] space_dual referencde to dual space
        /// @param[in] coupling_vars information about coupled variables in primal space
        /// @param[in] coupling_vars_dual information about coupled variables in dual space
        /// @param[out] couplings LA couplings for primal LA objects
        /// @param[out] couplings_dual LA couplings for dual LA objects
        virtual void setup_LA_primal ( VectorSpace<typename LAD::DataType>& space, VectorSpace<typename LAD::DataType>& space_dual,
                                       std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                       Couplings<typename LAD::DataType>& couplings, Couplings<typename LAD::DataType>& couplings_dual ) = 0;

        /// \brief routine which is called in udpate() of DynamicMeshHandler to setup required LA objects for specific fe space in case of dual problem time loop
        /// @param[in] space reference to primal space
        /// @param[in] space_dual referencde to dual space
        /// @param[in] coupling_vars information about coupled variables in primal space
        /// @param[in] coupling_vars_dual information about coupled variables in dual space
        /// @param[out] couplings LA couplings for primal LA objects
        /// @param[out] couplings_dual LA couplings for dual LA objects
        virtual void setup_LA_dual ( VectorSpace<typename LAD::DataType>& space, VectorSpace<typename LAD::DataType>& space_dual,
                                     std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                     Couplings<typename LAD::DataType>& couplings, Couplings<typename LAD::DataType>& couplings_dual ) = 0;

        /// \brief routine which is called in udpate() of DynamicMeshHandler to setup required LA objects for specific fe space in case of error estimator problem time loop
        /// @param[in] space reference to primal space
        /// @param[in] space_dual referencde to dual space
        /// @param[in] coupling_vars information about coupled variables in primal space
        /// @param[in] coupling_vars_dual information about coupled variables in dual space
        /// @param[out] couplings LA couplings for primal LA objects
        /// @param[out] couplings_dual LA couplings for dual LA objects
        virtual void setup_LA_est ( VectorSpace<typename LAD::DataType>& space, VectorSpace<typename LAD::DataType>& space_dual,
                                    std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                    Couplings<typename LAD::DataType>& couplings, Couplings<typename LAD::DataType>& couplings_dual ) = 0;

        /// \brief routine which is called in update() of DynamicMeshHandler to set vector space pointer in application to active space 
        /// @paam[in] space pointer to active space object
        /// @paam[in] mode 1: primal space, -1: dual space
        virtual void set_active_space ( VectorSpace<typename LAD::DataType>* space, int mode ) = 0;

        /// \brief routine which is called in update() of DynamicMeshHandler to set mesh pointer in application to active mesh  
        /// @param[in] mesh pointer to active mesh
        virtual void set_active_mesh ( MeshPtr mesh ) = 0;

        /// \brief routine which is called in init_dmh() of DynamicMeshProblem to get initial points for mesh changes and pass them to dmh object 
        virtual void init_mesh_change_list ( ) = 0;

        /// \brief routine which is called in init_dmh() of DynamicMeshProblem to pass those vectors to dmh object, which have to be interpolated from \br
        /// one space to another one inside of update() of DynamicMeshHandler
        virtual void set_update_vectors ( ) = 0;

        /// \brief set index of active mesh
        /// @param[in] index active index

        virtual void set_active_mesh_index ( int index )
        {
            this->active_mesh_index_ = index;
        }

        /// \brief main intialization routine for dmh object. This has to be called in the very beginning of an application
        /// @paam[in] adapt_type type of mesh change points adaption. see set_adapt_type() of DynamiMeshHandler
        virtual void init_dmh ( std::string& adapt_type );

        /// \brief interpolate vector between different fe spaces
        /// @param[in] in_vec input vector for interpolation
        /// @param[in] in_index index of input space
        /// @param[out] out_vec output vector for interpolation
        /// @param[in] out_index index of output space
        /// @param[in] mode 1: interpolate primal vector, -1: interpolate dual vector 
        void interpolate_vector ( const typename LAD::VectorType& in_vec, int in_index, typename LAD::VectorType& out_vec, int out_index, int mode );

      protected:
        /// \brief routine which selects point for changing the mesh according to spatial error indicators and passes to the dmh object
        /// @param[in] adapt_type type of mesh change adaption, see set_adapt_type() of DynamiMeshHandler
        /// @param[in] adapt_counter iteration index of outer adaption loop
        /// @param[in] min_steps_for_mesh minimum number of time steps between two successive mesh change points
        /// @param[in] max_mesh_number maximum number of mesh change points
        /// @param[in] fixed_fraction_rate number of newly created mesh changes in case of adapt_type = FixedFraction
        /// @param[in] fixed_error_tol threshold when to crate mesh change point in case of adapt_type = Fixed Error
        /// @param[in] space_indicator spatial error indicators 
        /// @param[in,out] indicator_mesh_indices map indicaotr[t] to mesh index
        virtual void adapt_mesh_change_list ( std::string adapt_type, int adapt_counter,
                                              int min_steps_for_mesh, int max_mesh_number,
                                              int fixed_fraction_rate, typename LAD::DataType fixed_error_tol,
                                              const std::vector< std::vector< typename LAD::DataType> >& space_indicator, std::vector<int>& indicator_mesh_indices );

        int active_mesh_index_;

        DynamicMeshHandler<LAD, MESH_IMPL, DIM>* dmh_;

        TimeMesh<typename LAD::DataType> t_mesh_;
    };
}

#endif

