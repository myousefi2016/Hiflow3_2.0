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

/// @author Staffan Ronnas, Thomas Gengenbach, Philipp Gerstner

#include "mesh_tools.h"

#include <mpi.h>

// public headers
#include "common/macros.h"
#include "common/log.h"
#include "attributes.h"
#include "entity.h"
#include "iterator.h"
#include "partitioning.h"
#include "reader.h"
#include "mesh/geometric_tools.h"
#include "mesh_db_view.h"
#include "refined_mesh_db_view.h"
#include "mesh/mesh_database.h"

#include "mesh_pXest.h"
#include "pXest_utils.h"
#ifdef WITH_P4EST
#    include "p4est_extended.h"
#    include "p8est_extended.h"
#endif

// private headers
#include "communication.h"
#include "mpi_communication.h"

const int DEBUG_LEVEL = 1;

namespace hiflow
{

    using namespace mesh;
    // TODO implement p4est reader for non .inp files

    //////////////// General interface /////////////////////
    /// \details The file type is determined from the filename
    /// suffix. MPI communicator is necessary for reading parallel Vtk
    /// files. In this case, the returned mesh is only the local part
    /// of the mesh.
    ///
    /// \param filename  Name of the file to read.
    /// \param tdim      Topological dimension of mesh.
    /// \param gdim      Geometrical dimension of mesh.
    /// \param comm      MPI communicatior for parallel input. Can be 0 for sequential input.
    /// \return  Shared pointer to created mesh object.

    MeshPtr read_mesh_from_file ( const std::string& filename,
                                  TDim tdim,
                                  GDim gdim,
                                  const MPI_Comm* comm,
                                  std::vector<MasterSlave> period,
                                  IMPL impl )
    {
        // find position of last '.'
        int dot_pos = filename.find_last_of ( "." );
        assert ( dot_pos != static_cast < int > ( std::string::npos ) );
        std::string suffix = filename.substr ( dot_pos );

        // choose reader according to filename suffix
        Reader* reader;
        MeshBuilder* builder;
        assert ( impl == IMPL_P4EST || impl == IMPL_DBVIEW );

        if ( impl == IMPL_P4EST )
        {
#ifdef WITH_P4EST
            builder = new MeshPXestBuilder ( tdim, gdim, period );
#else
            LOG_ERROR ( "Not compiled with p4ests support!" );
            exit ( -1 );
#endif
        }
        else if ( impl == IMPL_DBVIEW )
            builder = new MeshDbViewBuilder ( tdim, gdim, period );

        if ( suffix == std::string ( ".inp" ) )
        {
            assert ( impl == IMPL_P4EST || impl == IMPL_DBVIEW );

            if ( impl == IMPL_P4EST )
            {
#ifdef WITH_P4EST
                reader = new UcdPXestReader ( builder );
#else
                LOG_ERROR ( "Not compiled with p4ests support!" );
                exit ( -1 );
#endif
            }
            else if ( impl == IMPL_DBVIEW )
                reader = new UcdReader ( builder );
        }
        else if ( suffix == std::string ( ".vtu" ) )
        {
            assert ( impl == IMPL_DBVIEW );
            reader = new VtkReader ( builder );
        }
        else if ( suffix == std::string ( ".pvtu" ) )
        {
            assert ( comm != 0 );
            assert ( impl == IMPL_DBVIEW );
            reader = new PVtkReader ( builder, *comm );
        }

        // try to read the mesh
        MeshPtr mesh;
        try
        {
            reader->read ( filename.c_str ( ), mesh );
            delete reader;
        }
        catch ( const ReadGridException& exc )
        {
            delete reader;

            std::cerr << "Failed to read " << filename << "\n";
            std::cerr << exc.what ( ) << "\n";

            throw exc;
        }
        MeshPtr shared_mesh ( mesh );

        // TODO fix segmentation fault that occurs when uncommenting the following line
        // delete builder;
        return shared_mesh;
    }

    void save_mesh ( std::string filename, const MeshPtr mesh, const MPI_Comm& comm )
    {
#ifdef WITH_HDF5
        int tdim = mesh->tdim();
        int gdim = mesh->gdim();
        std::size_t pos_suffix = filename.find_last_of ( "." );

        std::string filename_without_suffix = filename.substr ( 0, pos_suffix );

        assert ( !filename_without_suffix.empty ( ) );

        std::string suffix = ".h5";

        filename = filename_without_suffix + suffix;

        const Mesh* mesh_ptr = mesh.get ( );
        const MeshDbView* mesh_db_view_ptr;
        const RefinedMeshDbView* refined_mesh_ptr;
        const MeshPXest* pXest_ptr;

        {
            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "Mesh";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );
            write_array_parallel ( group_ptr, "tdim", &tdim, 1, comm );
            write_array_parallel ( group_ptr, "gdim", &gdim, 1, comm );
        }

        if ( ( pXest_ptr = dynamic_cast < const MeshPXest* > ( mesh_ptr ) ) != 0 )
        {
            int history_index = pXest_ptr->get_history_index ( );
            pXest_ptr->get_db ( )->save ( filename, comm );
            std::map< Id, MeshPtr> mesh_history = pXest_ptr->get_db ( )->get_mesh_history ( );
            int counter = 0;
            for ( std::map< Id, MeshPtr>::const_iterator it = mesh_history.begin ( ); it != mesh_history.end ( ); ++it )
            {
                std::stringstream filenamei;
                filenamei << filename_without_suffix << "." << counter << suffix;
                it->second->save ( filenamei.str ( ), comm );
                ++counter;
            }
            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "HistoryMeshs";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );
            write_array_parallel ( group_ptr, "num_history_meshes", &counter, 1, comm );
            write_array_parallel ( group_ptr, "history_index", &history_index, 1, comm );
        }
        else if ( ( refined_mesh_ptr = dynamic_cast < const RefinedMeshDbView* > ( mesh_ptr ) ) != 0 )
        {
            refined_mesh_ptr->MeshDbView::save ( filename, comm );
            refined_mesh_ptr->get_db ( )->save ( filename, comm );
            std::vector<ConstMeshPtr> history = refined_mesh_ptr->get_all_ancestors ( );
            int num_meshes = history.size ( );
            for ( int i = 0; i < num_meshes; ++i )
            {
                std::stringstream filenamei;
                filenamei << filename_without_suffix << "." << i << suffix;
                const Mesh* conv_ptr = history[i].get ( );
                assert ( dynamic_cast < const MeshDbView* > ( conv_ptr ) != NULL );
                dynamic_cast < const MeshDbView* > ( conv_ptr )->MeshDbView::save ( filenamei.str ( ), comm );
            }

            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "RefinedMeshDbView";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );
            write_array_parallel ( group_ptr, "num_history_meshes", &num_meshes, 1, comm );

        }
        else if ( ( mesh_db_view_ptr = dynamic_cast < const MeshDbView* > ( mesh_ptr ) ) != 0 )
        {
            mesh_db_view_ptr->MeshDbView::save ( filename, comm );
            mesh_db_view_ptr->get_db ( )->save ( filename, comm );
        }
        else
        {
            mesh->Mesh::save ( filename, comm );
        }
#else
        throw "HiFlow was not compiled with HDF5 support!\n";
#endif
    }

    mesh::MeshPtr load_mesh ( std::string filename, const MPI_Comm& comm, mesh::IMPL impl )
    {
#ifdef WITH_HDF5
        std::size_t pos_suffix = filename.find_last_of ( "." );

        std::string filename_without_suffix = filename.substr ( 0, pos_suffix );

        assert ( !filename_without_suffix.empty ( ) );

        std::string suffix = ".h5";
        filename = filename_without_suffix + suffix;

        int* tdimptr;
        int* gdimptr;
        int tdim, gdim;
        {
            H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "Mesh";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );
            int dummy;
            read_array_parallel ( group_ptr, "tdim", tdimptr, dummy, comm );
            assert(dummy == 1);
            read_array_parallel ( group_ptr, "gdim", gdimptr, dummy, comm );
            assert(dummy == 1);
        }
        tdim = *tdimptr;
        gdim = *gdimptr;
        delete tdimptr;
        delete gdimptr;

        MeshPtr mesh;

        if ( impl == mesh::IMPL_P4EST )
        {

            MeshPXestDatabasePtr pXest_database ( new MeshPXestDatabase ( tdim, gdim ) );
            pXest_database->load ( filename, comm );
            int num_meshes = 0;
            int history_index = -1;

            {
                H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );
                //SETTING UP HDF5 GROUP
                std::stringstream groupname;
                groupname << "HistoryMeshs";
                H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );
                int dummy;
                int* buffer;
                read_array_parallel ( group_ptr, "num_history_meshes", buffer, dummy, comm );
                assert(dummy == 1);
                num_meshes = *buffer;
                delete buffer;
                assert ( num_meshes > 0 );

                read_array_parallel ( group_ptr, "history_index", buffer, dummy, comm );
                assert(dummy == 1);
                history_index = *buffer; 
                delete buffer;
                assert ( history_index > -1 );
            }

            for ( int i = 0; i < num_meshes; ++i )
            {
                MeshPXest* pXest_ptr = new MeshPXest ( tdim, gdim );
                std::stringstream filenamei;
                filenamei << filename_without_suffix << "." << i << suffix;
                pXest_ptr->load ( filenamei.str ( ), comm );
                pXest_ptr->set_db ( pXest_database );

                int temp_history_index = pXest_ptr->get_history_index ( );
                if ( temp_history_index == history_index )
                    mesh = pXest_ptr;

                pXest_database->set_mesh ( MeshPtr ( pXest_ptr ), temp_history_index );

            }
            if(pXest_database->get_layer_width ( ) > 0)
                compute_ghost_cells_pXest ( *mesh, comm, pXest_database->get_layer_width ( ) );
        }
        else if ( impl == mesh::IMPL_REFINED )
        {
            MeshDatabasePtr database ( new MeshDatabase ( tdim, gdim ) );
            database->load ( filename, comm );
            RefinedMeshDbView* refined_mesh_ptr = new RefinedMeshDbView ( tdim, gdim );
            refined_mesh_ptr->MeshDbView::load ( filename, comm );
            refined_mesh_ptr->set_db ( database );

            int num_meshes = 0;
            {
                H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );
                //SETTING UP HDF5 GROUP
                std::stringstream groupname;
                groupname << "RefinedMeshDbView";
                H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );
                int dummy;
                int* buffer;
                read_array_parallel ( group_ptr, "num_history_meshes", buffer, dummy, comm );
                assert(dummy == 1);
                num_meshes = *buffer;
                delete buffer;
            }
            std::vector<ConstMeshPtr> history_meshes ( num_meshes );
            for ( int i = 0; i < num_meshes; ++i )
            {
                std::stringstream filenamei;
                filenamei << filename_without_suffix << "." << i << suffix;

                MeshDbView* mesh_db_view_ptr = new MeshDbView ( tdim, gdim );
                mesh_db_view_ptr->load ( filenamei.str ( ), comm );
                mesh_db_view_ptr->set_db ( database );
                history_meshes[i] = mesh_db_view_ptr;
            }
            refined_mesh_ptr->set_history ( history_meshes );

            mesh = refined_mesh_ptr;
        }
        else if ( impl == mesh::IMPL_DBVIEW )
        {
            MeshDatabasePtr database ( new MeshDatabase ( tdim, gdim ) );
            database->load ( filename, comm );
            MeshDbView* mesh_db_view_ptr = new MeshDbView ( tdim, gdim );
            mesh_db_view_ptr->load ( filename, comm );
            mesh_db_view_ptr->set_db ( database );
            mesh = mesh_db_view_ptr;
        }

        return mesh;
#else
        throw "HiFlow was not compiled with HDF5 support!\n";
#endif
    }
    /// \details Subfunction that returns a distributed mesh, given a master_mesh with corresponding partitioning

    mesh::MeshPtr create_distributed_mesh ( const mesh::MeshPtr master_mesh,
                                            const int master_rank,
                                            const MPI_Comm& comm,
                                            std::vector<int>& partitioning,
                                            mesh::IMPL impl )
    {
        int my_rank = -1, num_ranks = -2;
        MPI_Comm_rank ( comm, &my_rank );
        MPI_Comm_size ( comm, &num_ranks );
        assert ( 0 <= my_rank );
        assert ( 0 <= master_rank );
        assert ( master_rank < num_ranks );
        assert ( my_rank < num_ranks );

        // check that master_mesh is non-null iff we are on master process
        assert ( my_rank != master_rank || master_mesh != 0 );
        assert ( my_rank == master_rank || master_mesh == 0 );
        assert ( impl == mesh::IMPL_DBVIEW || impl == mesh::IMPL_P4EST );

        // Communicate period
        std::vector<MasterSlave> period ( 0 );
        if ( my_rank == master_rank )
        {
            period = master_mesh->get_period ( );
        }
        int num_period = period.size ( );
        MPI_Bcast ( &num_period, 1, MPI_INT, master_rank, comm );
        for ( int k = 0; k < num_period; ++k )
        {
            std::vector<double> period_dbl_values ( 3, 0. );
            int period_index;
            if ( my_rank == master_rank )
            {
                period_dbl_values[0] = period[k].master ( );
                period_dbl_values[1] = period[k].slave ( );
                period_dbl_values[2] = period[k].h ( );
                period_index = period[k].index ( );
            }
            MPI_Bcast ( &period_dbl_values[0], 3, MPI_DOUBLE, master_rank, comm );

            MPI_Bcast ( &period_index, 1, MPI_INT, master_rank, comm );
            if ( my_rank != master_rank )
                period.push_back ( MasterSlave ( period_dbl_values[0], period_dbl_values[1], period_dbl_values[2], period_index ) );
        }

        std::vector<EntityPackage> sent_entities;
        EntityPackage recv_entities;
        std::vector<EntityPackage> sent_facet_entities;
        EntityPackage recv_facet_entities;
        std::vector<int> num_entities_on_proc;
        std::vector<int> num_facet_entities_on_proc;

        if ( my_rank == master_rank )
        {
            const TDim tdim = master_mesh->tdim ( );

            // compute distribution of cells from the partitioning
            CellDistribution distribution;
            compute_cell_distribution ( num_ranks, partitioning, &distribution );

            assert ( static_cast < int > ( distribution.num_cells.size ( ) ) == num_ranks );

            // compute distribution of facets from cell distribution
            CellDistribution facet_distribution;
            compute_facet_from_cell_distribution ( *master_mesh, distribution, &facet_distribution );

            // prepare sent entities
            pack_distributed_entities ( *master_mesh, tdim, distribution, sent_entities );
            pack_distributed_entities ( *master_mesh, tdim - 1, facet_distribution, sent_facet_entities );

            // use num_cells vector from distribution for num_entities_on_proc
            num_entities_on_proc.swap ( distribution.num_cells );
            num_facet_entities_on_proc.swap ( facet_distribution.num_cells );
            LOG_DEBUG ( 2, "num_entities_on_proc = "
                        << string_from_range ( num_entities_on_proc.begin ( ), num_entities_on_proc.end ( ) ) );
            LOG_DEBUG ( 2, "num_facet_entities_on_proc = "
                        << string_from_range ( num_facet_entities_on_proc.begin ( ), num_facet_entities_on_proc.end ( ) ) );
        }

        // Communicate cells
        MpiScatter scatter ( comm, master_rank, num_entities_on_proc );
        scatter.communicate ( sent_entities, recv_entities );

        LOG_DEBUG ( 2, "Received " << recv_entities.offsets.size ( ) << " cells on proc " << my_rank );

        // unpack the received part of the mesh into a builder
        MeshBuilder* builder;
        if ( impl == mesh::IMPL_DBVIEW )
        {
            builder = new MeshDbViewBuilder ( recv_entities.tdim, recv_entities.gdim, period );
        }
        if ( impl == mesh::IMPL_P4EST )
        {
#ifdef WITH_P4EST
            builder = new MeshPXestBuilder ( recv_entities.tdim, recv_entities.gdim, period );
#else
            LOG_ERROR ( "Not compiled with p4ests support!" );
            exit ( -1 );
#endif
        }

        unpack_entities ( recv_entities, *builder );

        // Communicate facets
        MpiScatter facet_scatter ( comm, master_rank, num_facet_entities_on_proc );
        facet_scatter.communicate ( sent_facet_entities, recv_facet_entities );

        LOG_DEBUG ( 2, "Received " << recv_facet_entities.offsets.size ( ) << " facets on proc " << my_rank );

        LOG_DEBUG ( 3, "Received material numbers "
                    << string_from_range ( recv_facet_entities.material_numbers.begin ( ),
                                           recv_facet_entities.material_numbers.end ( ) )
                    << " on proc " << my_rank << "\n" );

        // unpack the received facets into a builder
        unpack_entities ( recv_facet_entities, *builder );

        // build the local mesh
        MeshPtr recv_mesh ( builder->build ( ) );
        delete builder;
        return recv_mesh;
    }

    mesh::MeshPtr extract_local_mesh ( const mesh::MeshPtr master_mesh,
                                       const int master_rank,
                                       const MPI_Comm& comm,
                                       std::vector<int>& partitioning,
                                       mesh::IMPL impl )
    {
        int my_rank = -1, num_ranks = -2;
        MPI_Comm_rank ( comm, &my_rank );
        MPI_Comm_size ( comm, &num_ranks );
        assert ( 0 <= my_rank );
        assert ( 0 <= master_rank );
        assert ( master_rank < num_ranks );
        assert ( my_rank < num_ranks );

        // check that master_mesh is non-null iff we are on master process
        assert ( master_mesh != 0 );
        assert ( impl == mesh::IMPL_DBVIEW || impl == mesh::IMPL_P4EST );

        // get period
        std::vector<MasterSlave> period ( 0 );
        period = master_mesh->get_period ( );

        std::vector<EntityPackage> sent_entities;
        EntityPackage my_entities;
        std::vector<EntityPackage> sent_facet_entities;
        EntityPackage my_facet_entities;
        std::vector<int> num_entities_on_proc;
        std::vector<int> num_facet_entities_on_proc;

        const TDim tdim = master_mesh->tdim ( );

        // compute distribution of cells from the partitioning
        CellDistribution distribution;
        compute_cell_distribution ( num_ranks, partitioning, &distribution );

        assert ( static_cast < int > ( distribution.num_cells.size ( ) ) == num_ranks );

        // compute distribution of facets from cell distribution
        CellDistribution facet_distribution;
        compute_facet_from_cell_distribution ( *master_mesh, distribution, &facet_distribution );

        // prepare sent entities
        pack_distributed_entities ( *master_mesh, tdim, distribution, sent_entities );
        pack_distributed_entities ( *master_mesh, tdim - 1, facet_distribution, sent_facet_entities );

        // use num_cells vector from distribution for num_entities_on_proc
        num_entities_on_proc.swap ( distribution.num_cells );
        num_facet_entities_on_proc.swap ( facet_distribution.num_cells );

        LOG_DEBUG ( 2, "num_entities_on_proc = "
                    << string_from_range ( num_entities_on_proc.begin ( ), num_entities_on_proc.end ( ) ) );
        LOG_DEBUG ( 2, "num_facet_entities_on_proc = "
                    << string_from_range ( num_facet_entities_on_proc.begin ( ), num_facet_entities_on_proc.end ( ) ) );

        int my_num_entities = num_entities_on_proc[my_rank];
        int my_num_facet_entities = num_facet_entities_on_proc[my_rank];

        // Extract cells
        my_entities = sent_entities[my_rank];

        LOG_DEBUG ( 2, "Received " << my_entities.offsets.size ( ) << " cells on proc " << my_rank );

        // unpack the received part of the mesh into a builder
        MeshBuilder* builder;
        if ( impl == mesh::IMPL_DBVIEW )
        {
            builder = new MeshDbViewBuilder ( my_entities.tdim, my_entities.gdim, period );
        }
        if ( impl == mesh::IMPL_P4EST )
        {
#ifdef WITH_P4EST
            builder = new MeshPXestBuilder ( my_entities.tdim, my_entities.gdim, period );
#else
            LOG_ERROR ( "Not compiled with p4ests support!" );
            exit ( -1 );
#endif
        }

        unpack_entities ( my_entities, *builder );

        // Extract facets
        my_facet_entities = sent_facet_entities[my_rank];

        LOG_DEBUG ( 2, "Received " << my_facet_entities.offsets.size ( ) << " facets on proc " << my_rank );

        LOG_DEBUG ( 3, "Received material numbers "
                    << string_from_range ( my_facet_entities.material_numbers.begin ( ),
                                           my_facet_entities.material_numbers.end ( ) )
                    << " on proc " << my_rank << "\n" );

        // unpack the received facets into a builder
        unpack_entities ( my_facet_entities, *builder );

        // build the local mesh
        MeshPtr recv_mesh ( builder->build ( ) );
        delete builder;
        return recv_mesh;
    }

    /// \details Partition and distribute mesh that has been read on the master
    /// process to all processes in the given MPI communicator. Uses a MetisGraphPartitioner if
    /// available, otherwise a NaiveGraphPartitioner. The function returns the part of the mesh
    /// belonging to the local process.
    ///
    /// \param master_mesh  MeshPtr to the mesh to be distributed on master process, 0 on other processes.
    /// \param master_rank  Rank of master process.
    /// \param comm         MPI communicator used for communication.
    /// \return  Shared pointer to local part of the distributed mesh.

    MeshPtr partition_and_distribute ( const MeshPtr master_mesh,
                                       const int master_rank,
                                       const MPI_Comm& comm,
                                       int* uniform_ref_steps,
                                       IMPL impl )
    {
        assert ( impl == IMPL_DBVIEW || impl == IMPL_P4EST );

        if ( impl == IMPL_P4EST )
        {
#ifdef WITH_P4EST
            assert ( uniform_ref_steps != NULL );
            return partition_and_distribute_pXest ( master_mesh, master_rank, comm, *uniform_ref_steps );
#else
            LOG_ERROR ( "Not compiled with p4ests support!" );
            exit ( -1 );
#endif
        }

#ifdef WITH_PARMETIS
        ParMetisGraphPartitioner gpartitioner;
#else
#    ifdef WITH_METIS
        MetisGraphPartitioner gpartitioner;
#    else
        NaiveGraphPartitioner gpartitioner;
#    endif
#endif

        if ( uniform_ref_steps != NULL )
        {
            *uniform_ref_steps = 0;
        }
        return partition_and_distribute_dbview ( master_mesh,

                                                 master_rank,
                                                 comm,
                                                 &gpartitioner );

    }

    MeshPtr partition_and_distribute ( const MeshPtr master_mesh,
                                       const int master_rank,
                                       const MPI_Comm& comm,
                                       const GraphPartitioner* gpartitioner,
                                       int* uniform_ref_steps,
                                       IMPL impl )
    {
        assert ( impl == IMPL_DBVIEW );
        if ( impl == IMPL_DBVIEW )
        {
            if ( uniform_ref_steps != NULL )
            {
                *uniform_ref_steps = 0;
            }
            return partition_and_distribute_dbview ( master_mesh,
                                                     master_rank,
                                                     comm,
                                                     gpartitioner );
        }
    }

    MeshPtr partition_and_distribute ( const MeshPtr master_mesh,
                                       const int master_rank,
                                       const MPI_Comm& comm,
                                       const mesh::MeshPartitioner* mpartitioner,
                                       int* uniform_ref_steps,
                                       IMPL impl )
    {
        assert ( impl == IMPL_DBVIEW );
        if ( impl == IMPL_DBVIEW )
        {
            if ( uniform_ref_steps != NULL )
            {
                *uniform_ref_steps = 0;
            }
            return partition_and_distribute_dbview ( master_mesh,
                                                     master_rank,
                                                     comm,
                                                     mpartitioner );
        }
    }

    /// \details The ghost cells are computed in two steps. First, a
    /// search to identify the vertices shared with other processes is
    /// performed, and the results are stored in the shared_verts
    /// object. Second, the cells containing shared vertices are
    /// exchanged between all neighboring processes (i.e. processes
    /// that share at least one vertex).
    ///
    /// The returned mesh object contains an integer attribute named
    /// "__sub_domain__" which contains the owner process number for
    /// all cells. Another integer attribute "__remote_index__"
    /// contains the index of the cell for all ghost cells, and -1 for
    /// all local cells.
    ///
    /// Currently, this function is limited to work only with meshes
    /// whose concrete type is MeshDbView. The function will throw a
    /// "bad_cast" exception if this is not the case.
    ///
    /// \param local_mesh            The mesh on the current process.
    /// \param comm                  MPI communicator over which communication is performed.
    /// \param[in,out] shared_verts  Table with information about shared vertices.
    /// \return  Shared pointer to mesh object containing ghost cells.

    MeshPtr compute_ghost_cells ( const Mesh& local_mesh, const MPI_Comm& comm, SharedVertexTable& shared_verts, IMPL impl, int layer_width )
    {
        assert ( impl == IMPL_DBVIEW || impl == IMPL_P4EST );

        if ( impl == IMPL_P4EST )
        {
#ifdef WITH_P4EST
            return compute_ghost_cells_pXest ( local_mesh, comm, layer_width );
#else
            LOG_ERROR ( "Not compiled with p4ests support!" );
            exit ( -1 );
#endif
        }
        if ( impl == IMPL_DBVIEW )
        {
            return compute_ghost_cells_dbview ( local_mesh, comm, shared_verts );
        }
        //assert (impl == IMPL_DBVIEW);
        //return compute_ghost_cells_dbview(local_mesh, comm, shared_verts );
    }

    bool interpolate_attribute ( const MeshPtr parent_mesh, const std::string parent_attribute_name, MeshPtr child_mesh )
    {
        // same as copy, but interpolate values from the neighbors
        assert ( child_mesh != 0 );
        assert ( parent_mesh != 0 );
        assert ( child_mesh->tdim ( ) == parent_mesh->tdim ( ) );

        for ( int i = 0; i <= parent_mesh->tdim ( ); ++i )
        {
            // TODO(Thomas): double for now, if not, return false ->
            // to check (dynamic_cast?!)
            AttributePtr double_parent_attribute;
            try
            {
                double_parent_attribute = parent_mesh->get_attribute ( parent_attribute_name, i );
            }
            catch ( MissingAttributeException )
            {
                // Attribute <name> is missing!
                std::cerr << "...for topological dimension " << i << ".\n";
                continue;
            }
            std::vector<double> double_attribute_vector ( child_mesh->num_entities ( i ), 0.0 );

            // set this, to use later
            double mean = 0.0; //double_parent_attribute->get_double_value(0);

            std::vector<EntityNumber> index_to_interpolate;

            for ( int j = 0; j < child_mesh->num_entities ( i ); ++j )
            {
                // get index of entity in parent mesh
                EntityNumber index_parent;
                Id child_id = child_mesh->get_id ( i, j );
                // assume same database -> ids fit!
                // if not, return false, if they don:t fit, but exist
                //-> undefined behavior...
                bool found = parent_mesh->find_entity ( i, child_id, &index_parent );

                if ( found )
                {
                    double_attribute_vector[j] = double_parent_attribute->get_double_value ( index_parent );
                }
                else
                {
                    // store index and interpolate later...
                    index_to_interpolate.push_back ( j );
                }
            }

            for ( std::vector<EntityNumber>::const_iterator iti = index_to_interpolate.begin ( );
                  iti != index_to_interpolate.end ( ); ++iti )
            {
                Entity entity_to_interpolate = child_mesh->get_entity ( i, *iti );

                // iterate over incidents to this entity, gather
                // attributes, and take mean value (perhaps consider
                // weighted mean as well)

                // treat this in a special way, since 0->0
                // connectivity doesn:t exist in mesh
                int counter;
                double gather;

                if ( i == 0 )
                {
                    IncidentEntityIterator iei = entity_to_interpolate.begin_incident ( child_mesh->tdim ( ) );
                    counter = 0;
                    gather = 0.0;
                    for (; iei != entity_to_interpolate.end_incident ( child_mesh->tdim ( ) ); ++iei )
                    {
                        VertexIdIterator vertex_id_it = iei->begin_vertex_ids ( );
                        for (; vertex_id_it != iei->end_vertex_ids ( ); ++vertex_id_it )
                        {
                            if ( entity_to_interpolate.id ( ) == *vertex_id_it )
                            {
                                continue;
                            }
                            EntityNumber vertex_index;
                            parent_mesh->find_entity ( 0, *vertex_id_it, &vertex_index );
                            Entity other_vertex = parent_mesh->get_entity ( 0, vertex_index );
                            double tmp;

                            try
                            {
                                ++counter;
                                other_vertex.get ( parent_attribute_name, &tmp );
                            }
                            catch ( MissingAttributeException )
                            {
                                --counter;
                                tmp = 0.0;
                            }
                            gather += tmp;
                        }
                    }
                }
                else
                {
                    IncidentEntityIterator iei = entity_to_interpolate.begin_incident ( i );
                    counter = 0;
                    gather = 0.0;
                    for (; iei != entity_to_interpolate.end_incident ( i ); ++iei )
                    {
                        double tmp;
                        try
                        {
                            iei->get ( parent_attribute_name, &tmp );
                            ++counter;
                        }
                        catch ( MissingAttributeException )
                        {
                            tmp = 0.0;
                        }
                        gather += tmp;
                    }
                }
                // reset mean if possible, take old value else
                if ( counter != 0 ) mean = gather / counter;
                double_attribute_vector[*iti] = mean;
            }
            assert ( static_cast < int > ( double_attribute_vector.size ( ) ) == child_mesh->num_entities ( i ) );
            AttributePtr double_attribute = AttributePtr ( new DoubleAttribute ( double_attribute_vector ) );
            child_mesh->add_attribute ( parent_attribute_name, i, double_attribute );
            LOG_DEBUG ( 1, "Added attribute " << parent_attribute_name << " of dimension " << i << " to mesh." );
        }
        return true;
    }

    mesh::MeshPtr read_partitioned_mesh ( const std::string& filename,
                                          mesh::TDim tdim,
                                          mesh::GDim gdim,
                                          const MPI_Comm& comm )
    {
        MeshDbViewBuilder mb ( tdim, gdim );

        // Read parallel mesh file on all processes.
        MeshPtr local_mesh;
        PVtkReader reader ( &mb, comm );
        reader.read ( filename.c_str ( ), local_mesh );

        // Read boundary file on all processes, to get boundary material numbers.
        int point_pos = filename.find_last_of ( '.' );
        const std::string boundary_filename = filename.substr ( 0, point_pos ) + std::string ( "_bdy.vtu" );

        MeshDbViewBuilder bdy_mb ( mb.get_db ( ) );
        VtkReader bdy_reader ( &bdy_mb );
        MeshPtr bmesh;
        bdy_reader.read ( boundary_filename.c_str ( ), bmesh );

        return local_mesh;
    }

    void set_default_material_number_on_bdy ( mesh::MeshPtr mesh, MaterialNumber default_value )
    {
        MeshPtr bdy_mesh = mesh->extract_boundary_mesh ( );

        set_default_material_number_on_bdy ( mesh, bdy_mesh, default_value );
    }

    void set_default_material_number_on_bdy ( mesh::MeshPtr mesh, mesh::MeshPtr bdy_mesh, MaterialNumber default_value )
    {
        mesh::TDim mesh_tdim = mesh->tdim ( );

        const EntityIterator facet_begin = bdy_mesh->begin ( mesh_tdim - 1 );
        const EntityIterator facet_end = bdy_mesh->end ( mesh_tdim - 1 );

        for ( EntityIterator it = facet_begin; it != facet_end; it++ )
        {
            int mesh_facet_index;
            bdy_mesh->get_attribute_value ( "_mesh_facet_index_", mesh_tdim - 1, it->index ( ), &mesh_facet_index );
            MaterialNumber former_number = mesh->get_material_number ( mesh_tdim - 1, mesh_facet_index );
            if ( former_number == -1 )
            {
                mesh->set_material_number ( mesh_tdim - 1, mesh_facet_index, default_value );
            }
        }

    }

    void adapt_boundary_to_function ( mesh::MeshPtr mesh, const BoundaryDomainDescriptor& bdd )
    {
        MeshPtr bdry_mesh = mesh->extract_boundary_mesh ( );
        // go through the boundary mesh and project every point on the
        // boundary given through the BoundaryDomainDescriptor
        // Then replace the initial vertex with the projected vertex
        for ( EntityIterator it = bdry_mesh->begin ( bdry_mesh->tdim ( ) );
              it != bdry_mesh->end ( bdry_mesh->tdim ( ) );
              ++it )
        {
            mesh::MaterialNumber mat = it->get_material_number ( );

            for ( IncidentEntityIterator it_vert = it->begin_incident ( 0 );
                  it_vert != it->end_incident ( 0 ); ++it_vert )
            {
                std::vector<Coordinate> vert ( mesh->gdim ( ) );
                it_vert->get_coordinates ( vert );
                std::vector<Coordinate> new_vert = project_point ( bdd, vert, mat );
                int vert_id = it_vert->index ( );
                bdry_mesh->replace_vertex ( new_vert, vert_id );
            }
        }
    }

    //////////////// MeshDbView Implementation /////////////////////
    /// \details Partition and distribute mesh that has been read on the master
    /// process to all processes in the given MPI communicator using a given graph
    /// partitioner. The function returns the part of the mesh belonging to
    /// the local process.
    ///
    /// \param master_mesh  Mesh to be distributed on master process, 0 on other processes.
    /// \param master_rank  Rank of master process.
    /// \param comm         MPI communicator used for communication.
    /// \param gpartitioner The graph partitioner to be employed,
    /// e.g. NaiveGraphPartitioner, MetisGraphPartitioner or ParMetisGraphPartitioner if available.
    /// \return  Shared pointer to local part of the distributed mesh.

    MeshPtr partition_and_distribute_dbview ( const MeshPtr master_mesh,
                                              const int master_rank,
                                              const MPI_Comm& comm,
                                              const GraphPartitioner* gpartitioner )
    {
        assert ( gpartitioner != 0 );

        int my_rank = -1, num_ranks = -2;
        MPI_Comm_rank ( comm, &my_rank );
        MPI_Comm_size ( comm, &num_ranks );
        assert ( 0 <= my_rank );
        assert ( 0 <= master_rank );
        assert ( master_rank < num_ranks );
        assert ( my_rank < num_ranks );

        // check that master_mesh is non-null iff we are on master process
        assert ( my_rank != master_rank || master_mesh != 0 );
        assert ( my_rank == master_rank || master_mesh == 0 );

        Graph mesh_graph;
        mesh_graph.clear ( );

        // only master proc computes the graph
        if ( my_rank == master_rank )
        {
            compute_dual_graph ( *master_mesh, &mesh_graph );
        }

        std::vector<int> partitioning;
        if ( ( gpartitioner->is_collective ( ) ) || ( my_rank == master_rank ) )
        {
            gpartitioner->partition_graph ( mesh_graph, master_rank, 0, 0, comm, &partitioning );
        }

        return create_distributed_mesh ( master_mesh, master_rank, comm, partitioning, IMPL_DBVIEW );
    }

    /// \details Partition and distribute mesh that has been read on the master
    /// process to all processes in the given MPI communicator, using a given mesh
    /// partitioner. The function returns the part of the mesh belonging to
    /// the local process.
    ///
    /// \param master_mesh  Mesh to be distributed on master process, 0 on other processes.
    /// \param master_rank  Rank of master process.
    /// \param comm         MPI communicator used for communication.
    /// \param mpartitioner The mesh partitioner to be employed.
    /// \return  Shared pointer to local part of mesh

    MeshPtr partition_and_distribute_dbview ( const MeshPtr master_mesh,
                                              const int master_rank,
                                              const MPI_Comm& comm,
                                              const MeshPartitioner* mpartitioner )
    {
        LOG_DEBUG ( 1, "Not yet implemented." );
        quit_program ( );

        assert ( mpartitioner != 0 );

        int my_rank = -1, num_ranks = -2;
        MPI_Comm_rank ( comm, &my_rank );
        MPI_Comm_size ( comm, &num_ranks );
        assert ( 0 <= my_rank );
        assert ( 0 <= master_rank );
        assert ( master_rank < num_ranks );
        assert ( my_rank < num_ranks );

        // check that master_mesh is non-null iff we are on master process
        assert ( my_rank != master_rank || master_mesh != 0 );
        assert ( my_rank == master_rank || master_mesh == 0 );

        std::vector<int> partitioning;

        if ( mpartitioner == 0 )
        {
            LOG_DEBUG ( 1, "Error: User must provide a mesh partitioner!" );
            quit_program ( );
        }

        mpartitioner->partition_mesh ( master_mesh, master_rank, /*num_ranks,*/ 0, 0, comm, &partitioning );
        return create_distributed_mesh ( master_mesh, master_rank, comm, partitioning, mesh::IMPL_DBVIEW );
    }

    MeshPtr compute_ghost_cells_dbview ( const Mesh& local_mesh, const MPI_Comm& comm, SharedVertexTable& shared_verts )
    {
        int rank = -1;
        MPI_Comm_rank ( comm, &rank );

        // TODO this function only works if the concrete type of the mesh is MeshDbView.
        // otherwise exception is thrown here
        const MeshDbView& mesh_db_view = dynamic_cast < const MeshDbView& > ( local_mesh );

        // exchange shared vertices with other processes
        update_shared_vertex_table ( mesh_db_view, comm, shared_verts );

        // objects needed for ghost cell communication
        const TDim tdim = local_mesh.tdim ( );
        std::vector<SubDomainId> sub_domains ( local_mesh.num_entities ( tdim ), rank );
        std::vector<EntityNumber> remote_indices ( local_mesh.num_entities ( tdim ), -1 );
        MeshDbViewBuilder builder ( mesh_db_view );

        // communicate ghost cells and their facets
        communicate_ghost_cells ( mesh_db_view, comm, shared_verts, builder, sub_domains, remote_indices );
        communicate_ghost_cell_facets ( mesh_db_view, comm, shared_verts, builder );
        MeshPtr mesh_with_ghost_cells ( builder.build ( ) );

        // add attributes to new mesh
        AttributePtr sub_domain_attr = AttributePtr ( new IntAttribute ( sub_domains ) );
        mesh_with_ghost_cells->add_attribute ( "_sub_domain_", tdim, sub_domain_attr );
        AttributePtr remote_index_attr = AttributePtr ( new IntAttribute ( remote_indices ) );
        mesh_with_ghost_cells->add_attribute ( "_remote_index_", tdim, remote_index_attr );

        return mesh_with_ghost_cells;
    }

    //////////////// MeshP4est Implementation /////////////////////

    MeshPtr partition_and_distribute_pXest ( const MeshPtr master_mesh,
                                             const int master_rank,
                                             const MPI_Comm& comm,
                                             int& uniform_ref_steps )
    {
#ifdef WITH_P4EST
        // *********************************************************************
        // 0. check some stuff
        int my_rank = -1, num_ranks = -2;
        MPI_Comm_rank ( comm, &my_rank );
        MPI_Comm_size ( comm, &num_ranks );
        assert ( 0 <= my_rank );
        assert ( 0 <= master_rank );
        assert ( master_rank < num_ranks );
        assert ( my_rank < num_ranks );

        // check that master_mesh is non-null iff we are on master process
        //assert ( my_rank != master_rank || master_mesh != 0 );
        //assert ( my_rank == master_rank || master_mesh == 0 );

        // cast pointer
        boost::intrusive_ptr<MeshPXest> master_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( master_mesh );

        // *********************************************************************
        // 0. check if all processes have read in the mesh
        int has_mesh = 0;
        int seq_rank = master_rank;
        bool parallel = false;
        if ( master_mesh != 0 )
        {
            has_mesh = 1;
        }
        int num_mesh = 0;
        MPI_Allreduce ( &has_mesh, &num_mesh, 1, MPI_INT, MPI_SUM, comm );
        if ( num_mesh == num_ranks )
        {
            // everyone has the coarse mesh available
            seq_rank = my_rank;
            parallel = true;
        }
        else
        {
            seq_rank = master_rank;
            parallel = false;
        }

        // *********************************************************************
        // 1. Create p4est forest structure -> yields initial partitioning

        // broadcast p4est connectivity data
        int num_conn_vertices;
        int num_conn_cells;
        int array1_size;
        int array2_size;
        std::vector<int> tree_to_vertices;
        std::vector<double> conn_vertices;
        int tdim;

        std::vector<int> broadcast_info ( 7, 0 );

        uniform_ref_steps = 0;
        MeshPtr ref_master_mesh;
        std::vector<int> graph_part;
        int part_size = 0;

        // If more processes than coarse cells are involved -> refine mesh uniformly until every process can have at least one cell
        if ( my_rank == seq_rank )
        {
            tdim = master_mesh->tdim ( );
            num_conn_cells = master_mesh_pXest->get_db ( )->get_num_conn_cells ( );
            LOG_DEBUG ( 1, "Number of coarse cells " << num_conn_cells );

            ref_master_mesh = master_mesh;

            while ( num_conn_cells < 8 * num_ranks )
            {
                /*int num_ref_steps = (int) std::ceil( 8 * std::log((double) num_ranks / (double) num_conn_cells)
                                    / std::log( static_cast<int>(std::pow(static_cast<double>(2), tdim)) ) ); */

                ref_master_mesh = ref_master_mesh->refine_uniform_seq ( 1 );
                master_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( ref_master_mesh );
                num_conn_cells = master_mesh_pXest->get_db ( )->get_num_conn_cells ( );
                LOG_DEBUG ( 1, "Number of coarse cells " << num_conn_cells );
                uniform_ref_steps++;
            }

            LOG_DEBUG ( 1, "Initial mesh was uniformly refined " << uniform_ref_steps << " times " );

            // use graph partitioner to optimize communication pattern of mesh -> obtain graph_part
            // this partitioning is used to partition the newly created forest
            Graph mesh_graph;
            mesh_graph.clear ( );

#    ifdef WITH_METIS
            MetisGraphPartitioner gpartitioner;
#    else
            NaiveGraphPartitioner gpartitioner;
#    endif
            compute_dual_graph ( *ref_master_mesh, &mesh_graph );
            gpartitioner.partition_graph ( mesh_graph, master_rank, 0, 0, comm, &graph_part );
            part_size = graph_part.size ( );

            // collect data for broadcasting
            num_conn_vertices = master_mesh_pXest->get_db ( )->get_num_conn_vertices ( );
            num_conn_cells = master_mesh_pXest->get_db ( )->get_num_conn_cells ( );

            // PERF_TODO: evtl copy vermeiden
            tree_to_vertices = *master_mesh_pXest->get_db ( )->get_tree_to_vertices ( );
            conn_vertices = *master_mesh_pXest->get_db ( )->get_conn_vertices ( );
            array1_size = tree_to_vertices.size ( );
            array2_size = conn_vertices.size ( );

            broadcast_info[0] = tdim;
            broadcast_info[1] = num_conn_cells;
            broadcast_info[2] = num_conn_vertices;
            broadcast_info[3] = array1_size;
            broadcast_info[4] = array2_size;
            broadcast_info[5] = uniform_ref_steps;
            broadcast_info[6] = part_size;
        }

        if ( !parallel )
        {
            // broadcast connectivity from master
            MPI_Bcast ( vec2ptr ( broadcast_info ), 7, MPI_INT, master_rank, comm );

            tdim = broadcast_info[0];
            num_conn_cells = broadcast_info[1];
            num_conn_vertices = broadcast_info[2];
            array1_size = broadcast_info[3];
            array2_size = broadcast_info[4];
            uniform_ref_steps = broadcast_info[5];
            part_size = broadcast_info[6];

            if ( my_rank != seq_rank )
            {
                tree_to_vertices.resize ( array1_size );
                conn_vertices.resize ( array2_size );
                graph_part.resize ( part_size );
            }
            MPI_Bcast ( vec2ptr ( tree_to_vertices ), array1_size, MPI_INT, master_rank, comm );
            MPI_Bcast ( vec2ptr ( conn_vertices ), array2_size, MPI_DOUBLE, master_rank, comm );
            MPI_Bcast ( &graph_part[0], part_size, MPI_INT, master_rank, comm );
        }
        // compute number of cells per processor, used by forest partitioner
        std::vector<int> cells_per_proc ( num_ranks, 0 );
        for ( int l = 0; l < graph_part.size ( ); ++l )
        {
            cells_per_proc[ graph_part[l] ]++;
        }

        for ( int l = 0; l < cells_per_proc.size ( ); ++l )
        {
            LOG_DEBUG ( 2, "part: " << string_from_range ( graph_part.begin ( ), graph_part.end ( ) ) );
            assert ( cells_per_proc[l] > 0 );
        }

        // build connectivity and forest on each process, get initial partition from forest
        p4est_connectivity_t* conn4;
        p4est_t* forest4;

        p8est_connectivity_t* conn8;
        p8est_t* forest8;

        std::vector<int> partitioning;
        treeId first_local_tree = -1;
        treeId last_local_tree = -2;

        ForestInitData init_data ( my_rank );

        std::vector<size_t> old_to_new_tree;
        if ( tdim == 2 )
        {
            // build connectivity
            LOG_DEBUG ( 1, "Build p4est conenctivity" );
            pXest_build_conn ( num_conn_vertices, num_conn_cells, conn_vertices, tree_to_vertices, conn4 );

            // reorder connectivity according to graph_part
            LOG_DEBUG ( 1, "Reorder p4est conenctivity" );
            pXest_reorder_conn ( graph_part, conn4, old_to_new_tree );

            // reorder entity to quad mapping
            if ( my_rank == seq_rank )
            {
                master_mesh_pXest->get_db ( )->permute_entity_to_quad_map ( old_to_new_tree, 2, 0 );
            }

            // build forest
            LOG_DEBUG ( 1, "Build p4est forest" );
            forest4 = p4est_new ( comm, conn4, QuadData_size ( tdim ), pXest_init_fn, &init_data );

            // partition forest according to graph_part
            LOG_DEBUG ( 1, "Partition p4est forest" );
            pXest_partition_forest ( forest4, 0, cells_per_proc );
            //p4est_partition_ext ( forest4, 1, NULL );

            // sort quadrant arrays in forest
            LOG_DEBUG ( 1, "Sort quad arrays in forest" );
            pXest_sort_quad_arrays_in_forest ( forest4 );

            // get partitioning of trees. CAUTION: tree_id != cell_index
            partitioning = pXest_get_partition_from_initial_forest ( forest4 );
            first_local_tree = pXest_get_first_local_treeId ( forest4 );
            last_local_tree = pXest_get_last_local_treeId ( forest4 );
        }
        if ( tdim == 3 )
        {
            // build connectivity
            LOG_DEBUG ( 1, "Build p4est conenctivity" );
            pXest_build_conn ( num_conn_vertices, num_conn_cells, conn_vertices, tree_to_vertices, conn8 );

            // reorder connectivity according to graph_part
            LOG_DEBUG ( 1, "Reorder p4est conenctivity" );
            pXest_reorder_conn ( graph_part, conn8, old_to_new_tree );

            // reorder entity to quad mapping
            if ( my_rank == seq_rank )
            {
                master_mesh_pXest->get_db ( )->permute_entity_to_quad_map ( old_to_new_tree, 3, 0 );
            }

            // build forest
            LOG_DEBUG ( 1, "Build p4est forest" );
            forest8 = p8est_new ( comm, conn8, QuadData_size ( tdim ), pXest_init_fn, &init_data );

            // partition forest according to graph_part
            LOG_DEBUG ( 1, "Partition p4est forest" );
            pXest_partition_forest ( forest8, 0, cells_per_proc );
            //p8est_partition_ext ( forest8, 1, NULL );

            // sort quadrant arrays in forest
            LOG_DEBUG ( 1, "Sort quad arrays in forest" );
            pXest_sort_quad_arrays_in_forest ( forest8 );

            // get partitioning of trees. CAUTION: tree_id != cell_index
            partitioning = pXest_get_partition_from_initial_forest ( forest8 );
            first_local_tree = pXest_get_first_local_treeId ( forest8 );
            last_local_tree = pXest_get_last_local_treeId ( forest8 );
        }

        // Some debug ouput
        if ( my_rank == master_rank )
        {
            LOG_DEBUG ( 2, "Graph Partitioning: " << string_from_range ( graph_part.begin ( ), graph_part.end ( ) ) );
            LOG_DEBUG ( 2, "Old_to_new_tree permutation: " << string_from_range ( old_to_new_tree.begin ( ), old_to_new_tree.end ( ) ) );
            LOG_DEBUG ( 2, "Forest Partitioning: " << string_from_range ( partitioning.begin ( ), partitioning.end ( ) ) );

            if ( DEBUG_LEVEL >= 2 )
            {
                EntityToQuadMap* ent2quad = master_mesh_pXest->get_db ( )->get_entity_to_quad_map ( tdim, 0 );
                for ( EntityToQuadMap::iterator it = ent2quad->begin ( ); it != ent2quad->end ( ); ++it )
                {
                    LOG_DEBUG ( 2, " Coarse Entity id " << it->first << " maps to tree id " << it->second.tree );
                }
            }
        }

        // set quad_to_entity map in forest
        treeId num_local_trees = last_local_tree - first_local_tree + 1;
        for ( treeId l = 0; l < num_local_trees; ++l )
        {
            if ( tdim == 2 )
            {
                p4est_tree_t* tree = pXest_get_tree_in_forest ( forest4, l + first_local_tree );
                p4est_quadrant_t* quad = pXest_get_local_quad_in_tree ( tree, 0 );
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );

                LOG_DEBUG ( 2, my_rank << ": cell id " << l << ", tree of quad " << q_ptr->get_tree_id ( ) << " owner " << q_ptr->get_owner_rank ( )
                            << " quad_pointer " << quad << " quad_dataptr " << q_ptr );
                q_ptr->set_cell_id ( l, 0 );
                q_ptr->set_remote_cell_id ( l, 0 );
            }
            if ( tdim == 3 )
            {
                p8est_tree_t* tree = pXest_get_tree_in_forest ( forest8, l + first_local_tree );
                p8est_quadrant_t* quad = pXest_get_local_quad_in_tree ( tree, 0 );
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );

                LOG_DEBUG ( 2, my_rank << ": cell id " << l << ", tree of quad " << q_ptr->get_tree_id ( ) << " owner " << q_ptr->get_owner_rank ( )
                            << " quad_pointer " << quad << " quad_dataptr " << q_ptr );

                q_ptr->set_cell_id ( l, 0 );
                q_ptr->set_remote_cell_id ( l, 0 );
            }
        }

        // *********************************************************************
        // 2. partition mesh database, as before
        MeshPtr recv_mesh;
        if ( !parallel )
        {
            // only master has coarse mesh
            recv_mesh = create_distributed_mesh ( master_mesh_pXest, master_rank, comm, graph_part /*partitioning*/, IMPL_P4EST );
        }
        else
        {
            // every process has coarse mesh
            recv_mesh = extract_local_mesh ( master_mesh_pXest, master_rank, comm, graph_part /*partitioning*/, IMPL_P4EST );
        }

        // *********************************************************************
        // 3. additional p4est stuff
        assert ( recv_mesh != NULL );
        boost::intrusive_ptr<MeshPXest> recv_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( recv_mesh );

        // pass connectivity and forest to mesh database
        if ( tdim == 2 )
        {
            recv_mesh_pXest->get_db ( )->set_p4est_conn ( conn4 );
            recv_mesh_pXest->get_db ( )->set_p4est_forest ( forest4 );
            recv_mesh_pXest->update_cell_index_in_quads ( true );

        }
        if ( tdim == 3 )
        {
            recv_mesh_pXest->get_db ( )->set_p8est_conn ( conn8 );
            recv_mesh_pXest->get_db ( )->set_p8est_forest ( forest8 );
            recv_mesh_pXest->update_cell_index_in_quads ( true );
        }

        // store mesh in database
        recv_mesh_pXest->get_db ( )->set_mesh ( recv_mesh, 0 );

        // distribute coarse_entity_to_quad maps for cells and pass to database
        EntityToQuadMap map;
        EntityToQuadMap* coarse_map_ptr;
        std::vector<Id> coarse_cell_ids ( partitioning.size ( ), -1 );
        std::vector<Id> local_cell_ids ( partitioning.size ( ), -1 );
        if ( my_rank == seq_rank )
        {
            //boost::intrusive_ptr<MeshPXest> master_mesh_pXest = boost::static_pointer_cast<MeshPXest> (master_mesh);
            coarse_map_ptr = master_mesh_pXest->get_db ( )->get_entity_to_quad_map ( tdim, 0 );

            for ( EntityIterator it = master_mesh_pXest->begin ( tdim ); it != master_mesh_pXest->end ( tdim ); ++it )
            {
                int jc = std::distance ( master_mesh_pXest->begin ( tdim ), it );
                assert ( jc < partitioning.size ( ) );
                coarse_cell_ids[jc] = it->id ( );
            }
            LOG_DEBUG ( 2, "Rank " << my_rank << ": coarse_cell_ids " << string_from_range ( coarse_cell_ids.begin ( ), coarse_cell_ids.end ( ) ) );
        }
        else
        {
            coarse_map_ptr = NULL;
        }

        // get local cell ids
        int local_cell_counter = 0;
        for ( int jc = 0; jc < partitioning.size ( ); ++jc )
        {
            if ( graph_part[jc] /*partitioning[jc]*/ == my_rank )
            {
                local_cell_ids[jc] = local_cell_counter;
                ++local_cell_counter;
            }
        }

        LOG_DEBUG ( 2, "Rank " << my_rank << ": local_cell_ids " << string_from_range ( local_cell_ids.begin ( ), local_cell_ids.end ( ) ) );

        // distribute entity_to_quad maps
        if ( !parallel )
        {
            pXest_distribute_coarse_maps ( comm, master_rank, graph_part /*partitioning*/, coarse_cell_ids, local_cell_ids, coarse_map_ptr, map );
        }
        else
        {
            pXest_extract_local_maps ( comm, master_rank, graph_part /*partitioning*/, coarse_cell_ids, local_cell_ids, coarse_map_ptr, map );
        }

        recv_mesh_pXest->get_db ( )->set_entity_to_quad_map ( tdim, 1, map );

        return recv_mesh;
#else
        return master_mesh;
#endif
    }

    MeshPtr compute_ghost_cells_pXest ( const Mesh& local_mesh, const MPI_Comm& comm, int layer_width )
    {
#ifdef WITH_P4EST
        if ( layer_width < 2 )
        {
            int rank;
            MPI_Comm_rank ( comm, &rank );
            if ( rank == 0 )
            {
                std::cout << "CAUTION: You have chosen a ghost layer width of 1. In order to make use of local mesh refinement, layer width has to be 2 !!!" << std::endl;
            }
        }

        ConstMeshPtr tmp_mesh = &local_mesh;
        boost::intrusive_ptr<const MeshPXest> local_mesh_pXest = boost::static_pointer_cast<const MeshPXest> ( tmp_mesh );
        return local_mesh_pXest->add_ghost_cells ( comm, layer_width );
#else
        LOG_ERROR ( "Not compiled with p4ests support!" );
        exit ( -1 );
#endif
    }

}
