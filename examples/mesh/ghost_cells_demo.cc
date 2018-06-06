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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "hiflow.h"

/// \author Staffan Ronnas

/// \brief This program demonstrates how to communicate the ghost
/// cells of a local mesh, and create an extended local mesh
/// containing the ghost cells received from the neighboring
/// partitions.

using namespace std;
using hiflow::LogKeeper;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

const int DEBUG_LEVEL = 1;

static const char* datadir = MESHES_DATADIR;

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    ostringstream rank_str;
    rank_str << rank;

    std::string debug_filename = string ( "dbg_output." ) + rank_str.str ( ) + string ( ".log" );
    std::ofstream debug_file ( debug_filename.c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );
    std::string filename = std::string ( datadir ) + std::string ( "unitsquare_8_cells.pvtu" );
    //    std::string filename = std::string(datadir) + std::string("lunge_fein_3.pvtu");
    //    std::string filename = std::string(datadir) + std::string("partitioned_lung.pvtu");
    //    std::string filename = std::string("./metis_partitioned_lung.pvtu");

    // read in mesh in parallel
    LOG_DEBUG ( 1, "Reading mesh " << filename );
    MeshDbViewBuilder mb ( tdim, gdim );
    ScopedPtr<Reader>::Type reader ( new PVtkReader ( &mb, MPI_COMM_WORLD ) );
    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    // NB: this should be fixed once reader only fills builder instead of returning MeshPtr
    MeshDbView* local_mesh = dynamic_cast < MeshDbView* > ( mesh.get ( ) );
    assert ( local_mesh != 0 );

    LOG_DEBUG ( 1, "Read mesh with " << local_mesh->num_entities ( tdim ) << " cells" );

    // compute shared vertex table
    SharedVertexTable shared_verts;
    update_shared_vertex_table ( *local_mesh, MPI_COMM_WORLD, shared_verts );

    LOG_DEBUG ( 1, "Computed shared vertex table. Process " << rank
                << " has " << shared_verts.num_neighbor_domains ( ) << " neighbor domains" );

    std::vector<SubDomainId> sub_domains ( local_mesh->num_entities ( tdim ), rank );
    std::vector<EntityNumber> remote_indices ( local_mesh->num_entities ( tdim ), -1 );

    // build local mesh with ghost cells, with same database as local mesh
    LOG_DEBUG ( 1, "Building ghost mesh" );
    ;
    MeshDbViewBuilder ghost_builder ( *local_mesh );

    // call helper function for communicating the ghost cells
    communicate_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts, ghost_builder, sub_domains, remote_indices );
    MeshPtr mesh_with_ghost_cells = ghost_builder.build ( );

    // add attributes
    AttributePtr sub_domain_attr = AttributePtr ( new IntAttribute ( sub_domains ) );
    mesh_with_ghost_cells->add_attribute ( "_sub_domain_", tdim, sub_domain_attr );
    AttributePtr remote_index_attr = AttributePtr ( new IntAttribute ( remote_indices ) );
    mesh_with_ghost_cells->add_attribute ( "_remote_index_", tdim, remote_index_attr );

    // write mesh
    PVtkWriter writer ( MPI_COMM_WORLD );
    writer.add_attribute ( "_sub_domain_", tdim );
    writer.add_attribute ( "_remote_index_", tdim );
    writer.write ( "ghost_mesh.pvtu", *mesh_with_ghost_cells );

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );

    MPI_Finalize ( );
    return 0;
}
