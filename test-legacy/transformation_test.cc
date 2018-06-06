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

/// \author Michael Schick

#include <string>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESH_DATADIR;

int main ( int argc, char *argv[] )
{

    // mesh
    const string filename = string ( datadir ) + string ( "one_tet.inp" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );
    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    FEManager<double> fe_mgr ( gdim, 0 );
    fe_mgr.set_mesh ( *mesh );

    for ( EntityIterator it = mesh->begin ( gdim ); it != mesh->end ( gdim ); ++it )
    {
        std::vector<double> cv;
        it->get_coordinates ( cv );

        std::cout << "==================================" << std::endl;
        std::cout << "Coordinates before Transformation:" << std::endl;
        std::cout << "==================================" << std::endl;
        if ( gdim == 2 )
        {
            for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
                std::cout << "ID: " << vtx << "  Coord x: " << cv[vtx * gdim] << "  Coord y: " << cv[vtx * gdim + 1] << std::endl;
        }
        else
        {
            for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
                std::cout << "ID: " << vtx << "  Coord x: " << cv[vtx * gdim] << "  Coord y: " << cv[vtx * gdim + 1] << "  Coord z: " << cv[vtx * gdim + 2] << std::endl;
        }
        std::cout << "==================================" << std::endl;

        std::vector<double> cv_at ( cv.size ( ) );
        double x, y, z;
        for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
        {
            std::vector<double> tmp ( gdim );
            if ( gdim == 2 )
            {
                fe_mgr.get_cell_transformation ( it->index ( ) )->inverse ( cv[vtx * gdim], cv[vtx * gdim + 1], x, y );
                tmp[0] = x;
                tmp[1] = y;
                cv_at[vtx * gdim] = x;
                cv_at[vtx * gdim + 1] = y;
            }
            else
            {
                fe_mgr.get_cell_transformation ( it->index ( ) )->inverse ( cv[vtx * gdim], cv[vtx * gdim + 1], cv[vtx * gdim + 2], x, y, z );
                tmp[0] = x;
                tmp[1] = y;
                tmp[2] = z;
                cv_at[vtx * gdim] = x;
                cv_at[vtx * gdim + 1] = y;
                cv_at[vtx * gdim + 2] = z;
            }
        }

        std::cout << "==================================" << std::endl;
        std::cout << "Coordinates after Transformation:" << std::endl;
        std::cout << "==================================" << std::endl;
        if ( gdim == 2 )
        {
            for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
                std::cout << "ID: " << vtx << "  Coord x: " << cv_at[vtx * gdim] << "  Coord y: " << cv_at[vtx * gdim + 1] << std::endl;
        }
        else
        {
            for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
                std::cout << "ID: " << vtx << "  Coord x: " << cv_at[vtx * gdim] << "  Coord y: " << cv_at[vtx * gdim + 1] << "  Coord z: " << cv_at[vtx * gdim + 2] << std::endl;
        }
        std::cout << "==================================" << std::endl;

        std::vector<double> cv_bk ( cv.size ( ) );

        for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
        {
            std::vector<double> tmp0 ( gdim );
            std::vector<double> tmp ( gdim );

            if ( gdim == 2 )
            {
                tmp0[0] = cv_at[vtx * gdim];
                tmp0[1] = cv_at[vtx * gdim + 1];
                tmp[0] = fe_mgr.get_cell_transformation ( it->index ( ) )->x ( tmp0 );
                tmp[1] = fe_mgr.get_cell_transformation ( it->index ( ) )->y ( tmp0 );
                cv_bk[vtx * gdim] = tmp[0];
                cv_bk[vtx * gdim + 1] = tmp[1];
            }
            else
            {
                tmp0[0] = cv_at[vtx * gdim];
                tmp0[1] = cv_at[vtx * gdim + 1];
                tmp0[2] = cv_at[vtx * gdim + 2];
                tmp[0] = fe_mgr.get_cell_transformation ( it->index ( ) )->x ( tmp0 );
                tmp[1] = fe_mgr.get_cell_transformation ( it->index ( ) )->y ( tmp0 );
                tmp[2] = fe_mgr.get_cell_transformation ( it->index ( ) )->z ( tmp0 );
                cv_bk[vtx * gdim] = tmp[0];
                cv_bk[vtx * gdim + 1] = tmp[1];
                cv_bk[vtx * gdim + 2] = tmp[2];
            }
        }

        std::cout << "==================================" << std::endl;
        std::cout << "Coordinates Back Transformation:" << std::endl;
        std::cout << "==================================" << std::endl;
        if ( gdim == 2 )
        {
            for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
                std::cout << "ID: " << vtx << "  Coord x: " << cv_bk[vtx * gdim] << "  Coord y: " << cv_bk[vtx * gdim + 1] << std::endl;
        }
        else
        {
            for ( int vtx = 0; vtx < it->num_vertices ( ); ++vtx )
                std::cout << "ID: " << vtx << "  Coord x: " << cv_bk[vtx * gdim] << "  Coord y: " << cv_bk[vtx * gdim + 1] << "  Coord z: " << cv_bk[vtx * gdim + 2] << std::endl;
        }
        std::cout << "==================================" << std::endl;

    }
    return 0;

}
