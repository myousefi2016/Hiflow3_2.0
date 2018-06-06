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

/// \author Staffan Ronnas

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <mpi.h>

#include "hiflow.h"
#include "../src/mesh/communication.h"
#include "../src/mesh/mpi_communication.h"
#include "test.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESH_DATADIR;
const int DEBUG_LEVEL = 1;

int main ( int argc, char** argv )
{
    // Seems to work sequentially
    // TODO: try it in parallel!
    MPI_Init ( &argc, &argv );

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    int num_partitions = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_partitions );

    const int master_rank = 0;

    ostringstream rank_str;
    rank_str << rank;

    std::string debug_filename = string ( "dbg_output." ) + rank_str.str ( ) + string ( ".log" );
    std::ofstream debug_file ( debug_filename.c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    EntityPackage sent_entities, received_entities;
    std::vector<int> num_entities_on_proc;

    if ( rank == master_rank )
    {
        std::string filename = std::string ( datadir ) + std::string ( "dfg_bench3d_rect.inp" );
        LOG_DEBUG ( 1, "Reading mesh from " << filename << " on process " << rank );

        MeshDbViewBuilder read_builder ( tdim, gdim );
        //        VtkReader reader(&read_builder);
        UcdReader reader ( &read_builder );

        MeshPtr master_mesh;
        reader.read ( filename.c_str ( ), master_mesh );
        assert ( master_mesh != 0 );

        LOG_DEBUG ( 1, "Mesh has " << master_mesh->num_entities ( tdim ) << " cells" );

        int entities_per_proc = master_mesh->num_entities ( tdim ) / num_partitions;
        int total_entities = 0;
        for ( int p = 0; p < num_partitions - 1; ++p )
        {
            num_entities_on_proc.push_back ( entities_per_proc );
            total_entities += entities_per_proc;
        }
        num_entities_on_proc.push_back ( master_mesh->num_entities ( tdim ) - total_entities );

        pack_entities ( *master_mesh, tdim, &sent_entities );
    }

    LOG_DEBUG ( 1, "Entering communication on process " << rank );
    MpiScatter scatter ( MPI_COMM_WORLD, master_rank, num_entities_on_proc );
    scatter.communicate ( &sent_entities, &received_entities );

    LOG_DEBUG ( 1, "Received entity_package on process " << rank << ", num_coordinates = " << received_entities.coords.size ( ) );

    LOG_DEBUG ( 1, "Building received mesh part on process " << rank );
    MeshDbViewBuilder builder ( tdim, gdim );
    unpack_entities ( received_entities, builder );

    MeshPtr mesh = builder.build ( );

    assert ( mesh != 0 );

    LOG_DEBUG ( 1, "Iterating over vertices" );
    LOG_DEBUG ( 1, "Iterating over cells" );
    // Check that the same vertex id:s are obtained directly via
    // begin_vertex_ids() / end_vertex_ids() and by incidence iteration.
    for ( EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {
        SortedArray<Id> vertex_ids ( it->begin_vertex_ids ( ), it->end_vertex_ids ( ) );
        SortedArray<Id> incident_vertex_ids;
        for ( IncidentEntityIterator it2 = it->begin_incident ( 0 ); it2 != it->end_incident ( 0 ); ++it2 )
        {
            incident_vertex_ids.insert ( it2->id ( ) );
        }
        if ( vertex_ids.size ( ) != incident_vertex_ids.size ( ) )
        {
            LOG_DEBUG ( 1, "Not the same number of incident vertices: vertex_ids.size() == " << vertex_ids.size ( )
                        << "incident_vertex_ids.size() == " << incident_vertex_ids.size ( ) );
        }
        TEST_EQUAL ( vertex_ids.size ( ), incident_vertex_ids.size ( ) );
        for ( int i = 0; i < static_cast < int > ( vertex_ids.size ( ) ); ++i )
        {
            if ( vertex_ids[i] != incident_vertex_ids[i] )
            {
                LOG_DEBUG ( 1, "Different vertex ids: vertex_ids[i] = " << vertex_ids[i] << ", "
                            << "incident_vertex_ids[i] = " << incident_vertex_ids[i] );
            }
        }
        TEST_EQUAL ( vertex_ids.size ( ), incident_vertex_ids.size ( ) );
    }

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );

    MPI_Finalize ( );
    return 0;
}
