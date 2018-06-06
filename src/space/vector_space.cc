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

#include "vector_space.h"

#include "fem/fecellwisemanager.h"
#include "fem/fetype.h"
#include "common/log.h"

/// @author Martin Baumann, Chandramowli Subramanian, Michael Schick

namespace hiflow
{

    const int DEBUG_LEVEL = 0;

    /// standard constructor

    template<class DataType>
    VectorSpace<DataType>::VectorSpace ( )
    {
        dof_ = NULL;
        fe_manager_ = NULL;
        mesh_ = NULL;
        comm_ = MPI_COMM_WORLD;
    }

    /// MPI constructor

    template<class DataType>
    VectorSpace<DataType>::VectorSpace ( const MPI_Comm& comm )
    {
        dof_ = NULL;
        fe_manager_ = NULL;
        mesh_ = NULL;
        comm_ = comm;
    }

    /// destructor

    template<class DataType>
    VectorSpace<DataType>::~VectorSpace ( )
    {
        this->Clear ( );
    }

    /// Initialize vector space for scalar problems
    /// \param[in] degree defines the polynomial degree of the Lagrangian elements

    template<class DataType>
    void VectorSpace<DataType>::Init ( int degree, Mesh& mesh, hiflow::doffem::DOF_ORDERING order )
    {
        std::vector<int> degree_vec ( 1, degree );
        std::vector<bool> cg_vec ( 1, true );

        LOG_DEBUG ( 3, "This is mesh " << &mesh << " in VectorSpace::Init(int degree, Mesh& mesh)\n" );

        Init ( degree_vec, mesh, cg_vec, order );
    }

    /// Initialize vector space for scalar problems
    /// \param[in] degree defines the polynomial degree of the Lagrangian elements
    ///                   for each scalar, i.e. degree.size() should be the number
    ///                   of scalars used -> nval
    /// \param[in] is_cg sets a discontinuous ansatz space if false and continuous else.
    ///                  By default a continuous Galerkin method is used.

    template<class DataType>
    void VectorSpace<DataType>::Init ( std::vector<int> degree, Mesh& mesh, std::vector<bool> is_cg, hiflow::doffem::DOF_ORDERING order )
    {
        // clear
        this->Clear ( );

        LOG_INFO ( "Finite Elements [Q,P]", string_from_range ( degree.begin ( ), degree.end ( ) ) );

        // set parameters
        int nvar = degree.size ( );
        interminable_assert ( nvar > 0 );

        int dim = mesh.tdim ( );
        mesh_ = &mesh;

        //fill missing variables up with continuous ansatz space
        if ( is_cg.size ( ) < nvar )
        {
            is_cg.resize ( nvar, true );
        }

        assert ( is_cg.size ( ) == nvar );
        LOG_DEBUG ( 3, "This is mesh " << &mesh << " in VectorSpace::Init(std::vector<int> degree, Mesh& mesh)\n" );

        // create finite element manager
        fe_manager_ = new FEManager ( dim, nvar );
        fe_manager_->set_mesh ( mesh );

        std::vector<int> data ( 1, 0 );
        for ( int var = 0; var < fe_manager_->get_nb_var ( ); ++var )
        {
            fe_manager_->set_ca ( var, is_cg[var] );
            data[0] = degree[var];
            fe_manager_->init_fe_tank ( var, FEType::LAGRANGE, data );
        }
        LOG_INFO ( "Finite Element Type", "Lagrange" );

        // create degrees of freedom
        dof_ = new DofPartition; // (sequential and parallel)

        dof_->set_mpi_comm ( comm_ );
        dof_->set_mesh ( mesh_ );
        dof_->set_fe_manager ( fe_manager_ );

        dof_->number ( order );
    }

    /// Initialize p-refinement ready vector space for scalar problems
    /// \param[in] degree defines the polynomial degree of the Lagrangian elements
    ///                   for scalar 'i' and cell 'j' in case of degree[i][j]

    template<class DataType>
    void VectorSpace<DataType>::Init_p ( std::vector<std::vector<int> > degree, mesh::Mesh& mesh, hiflow::doffem::DOF_ORDERING order )
    {
        this->Clear ( );
        // set parameters
        int nvar = degree.size ( );
        interminable_assert ( nvar > 0 );
        for ( int var = 0; var < degree.size ( ); ++var )
            interminable_assert ( mesh.num_entities ( mesh.tdim ( ) ) == degree[var].size ( ) )

            int dim = mesh.tdim ( );
        mesh_ = &mesh;

        LOG_DEBUG ( 3, "This is mesh " << &mesh << " in VectorSpace::Init(std::vector<std::vector<int> > degree, Mesh& mesh)\n" );
        // create finite element manager
        fe_manager_ = new typename doffem::FECellwiseManager<DataType>( dim, nvar );
        fe_manager_->set_mesh ( mesh );
        for ( int var = 0; var < fe_manager_->get_nb_var ( ); ++var )
        {
            fe_manager_->init_fe_tank ( var, FEType::LAGRANGE, degree.at ( var ) );
        }
        LOG_INFO ( "Finite Element Type", "Lagrange" );

        // create degrees of freedom
        dof_ = new DofPartition;
        dof_->set_mpi_comm ( comm_ );
        dof_->set_mesh ( mesh_ );
        dof_->set_fe_manager ( fe_manager_ );
        dof_->number ( order );
    }

    /// \details Compute Solution via transformation to reference cell and
    ///          summation over the local indices (fem ordering strategy)
    ///          of the shapefuntions

    template<class DataType>
    template<class VectorType>
    DataType VectorSpace<DataType>::get_solution_value ( int var, const MeshEntity& cell,
                                                         const std::vector<DataType>& coord,
                                                         const VectorType& sol ) const
    {
        DataType sum = 0.0;
        // Global DoF Ids on the given mesh cell
        std::vector<int> global_dof_ids;
        GetDofIndices ( var, cell, &global_dof_ids );

        // Determine corresponding coordinates on reference cell via transformation
        std::vector<DataType> coord_ref ( coord.size ( ) );
        if ( get_dim ( ) == 2 )
            GetCellTransformation ( cell ).inverse ( coord[0], coord[1],
                                                     coord_ref[0], coord_ref[1] );
        else
            GetCellTransformation ( cell ).inverse ( coord[0], coord[1], coord[2],
                                                     coord_ref[0], coord_ref[1], coord_ref[2] );
        std::vector<DataType> weights ( global_dof_ids.size ( ) );
        fe_manager ( ).get_fe_on_cell ( cell.index ( ), var )->N ( coord_ref, weights );
        // Summation over weights of shapefunctions on reference cell
        for ( int i_loc = 0; i_loc < global_dof_ids.size ( ); ++i_loc )
        {
            sum += sol[global_dof_ids[i_loc]] * weights[i_loc];
        }
        return sum;
    }

    /// Clears allocated dof and femanager

    template<class DataType>
    void VectorSpace<DataType>::Clear ( )
    {
        if ( dof_ != NULL )
            delete dof_;
        dof_ = NULL;

        if ( fe_manager_ != NULL )
            delete fe_manager_;
        fe_manager_ = NULL;

        // don't delete mesh object as someone else owns it
        mesh_ = NULL;
    }

    template<class DataType>
    std::vector<int> VectorSpace<DataType>::get_dof_func ( std::vector<int> vars ) const
    {
        // check vars vector for correct input
        assert ( vars.size ( ) >= 0 );
        assert ( vars.size ( ) < this->get_nb_var ( ) );

        // if standard configuration is desired, the mapping is computed for all
        // variables
        if ( vars.size ( ) == 0 )
        {
            vars.resize ( this->get_nb_var ( ) );
            for ( size_t i = 0; i != vars.size ( ); ++i )
            {
                vars[i] = i;
            }
        }

        // consecutive numbering of vars
        std::sort ( vars.begin ( ), vars.end ( ) );
        std::map<int, int> vars2consec;
        for ( int i = 0; i < vars.size ( ); ++i )
        {
            vars2consec[vars[i]] = i;
        }

        std::vector<int> dof_func;
        std::map<int, int> dof_map;
        //dof_func_.resize ( space_.dof ( ).ndofs_on_sd ( rank_ ) );
        //const int my_dof_offset = space_.dof ( ).get_my_dof_offset ( );
        for ( hiflow::mesh::EntityIterator it = this->mesh_->begin ( this->mesh_->tdim ( ) ),
              e_it = this->mesh_->end ( this->mesh_->tdim ( ) );
              it != e_it;
              ++it )
        {
            // get global dof indices on current cell
            for ( int v_i = 0, v_i_e = vars.size ( ); v_i != v_i_e; ++v_i )
            {
                const int var = vars[v_i];
                std::vector<int> var_indices;
                this->dof ( ).get_dofs_on_cell ( var, it->index ( ), var_indices );
                // reduce dof indices to only subdomain indices
                for ( int i = 0, i_e = var_indices.size ( );
                      i != i_e;
                      ++i )
                {
                    if (
                         this->dof ( ).is_dof_on_sd
                         (
                           var_indices[i]
                           )
                         )
                    {
                        dof_map[var_indices[i]] = vars2consec[var];
                    }
                }
            }
        }

        // create return vector
        dof_func.resize ( dof_map.size ( ) );
        size_t i = 0;
        for ( std::map<int, int>::const_iterator it = dof_map.begin ( ),
              it_e = dof_map.end ( ); it != it_e; ++it )
        {
            dof_func[i] = it->second;
            ++i;
        }

        dof_map.clear ( );
        return dof_func;
    }

    /// Prints dof information.

    template<class DataType>
    void VectorSpace<DataType>::Print ( ) const
    {
        dof_->print_numer ( );
    }

    // template instantiation
    template class VectorSpace<double>;
    template class VectorSpace<float>;

    template double VectorSpace<double>::get_solution_value<std::vector<double> > (
            int var, const MeshEntity& cell,
            const std::vector<Coordinate>& coord,
            const std::vector<double>& sol ) const;
    template float VectorSpace<float>::get_solution_value<std::vector<float> > (
            int var, const MeshEntity& cell,
            const std::vector<float>& coord,
            const std::vector<float>& sol ) const;
} // namespace hiflow
