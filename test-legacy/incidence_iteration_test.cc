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

/// \author Thomas Gengenbach

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const int DEBUG_LEVEL = 3;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESH_DATADIR;

void read_connectivities ( const char* filename, vector< pair<TDim, TDim> >& dims,
                           vector< vector< vector<Id> > >& incidence_lists )
{
    ifstream file ( filename );
    assert ( file.good ( ) );

    std::string line;

    while ( getline ( file, line ) )
    {
        if ( line.find ( "BEGIN" ) != string::npos ) break;
    }

    while ( getline ( file, line ) )
    {
        if ( line[0] == '#' )
        {
            continue;
        }
        else if ( line.find ( "[" ) != string::npos )
        {
            istringstream dim_pair_stream ( line );
            int d1 = 0, d2 = 0;
            char s;
            dim_pair_stream >> s >> d1 >> d2 >> s;
            dims.push_back ( make_pair ( d1, d2 ) );
            incidence_lists.push_back ( vector< vector<Id> >( 0 ) );
        }
        else if ( line.find ( "END" ) != string::npos )
        {
            break;
        }
        else
        {
            istringstream stream ( line );
            istream_iterator<Id> end;
            vector<Id> id_list ( istream_iterator<Id>( stream ), end );
            incidence_lists.back ( ).push_back ( id_list );
        }
    }
}

void check_incidence_iteration ( TDim d1, TDim d2,
                                 const Mesh& mesh, const vector< vector<Id> >& correct_incidences )
{
    LOG_DEBUG ( 1, "Testing " << d1 << " -> " << d2 );
    TEST_EQUAL ( mesh.num_entities ( d1 ), static_cast < int > ( correct_incidences.size ( ) ) );
    int i = 0;
    for ( EntityIterator it = mesh.begin ( d1 ); it != mesh.end ( d1 ); ++it )
    {
        LOG_DEBUG ( 2, "Entities incident to entity " << it->id ( ) );
        const vector<Id>& correct_incidence = correct_incidences[i];
        unsigned int j = 0;

        for ( IncidentEntityIterator it2 = mesh.begin_incident ( *it, d2 );
              it2 != mesh.end_incident ( *it, d2 ); ++it2 )
        {
            LOG_DEBUG ( 2, *it2 );
        }

        for ( IncidentEntityIterator it2 = mesh.begin_incident ( *it, d2 );
              it2 != mesh.end_incident ( *it, d2 ); ++it2 )
        {

            TEST_LESS ( j, correct_incidence.size ( ) );
            TEST ( find ( correct_incidence.begin ( ), correct_incidence.end ( ), it2->id ( ) )
                   != correct_incidence.end ( ) );
            ++j;
        }
        TEST_EQUAL ( j, correct_incidence.size ( ) );
        ++i;
    }
}

int main ( int, char** )
{
    ofstream debug_file ( "incidence_test.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    string filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    vector< pair<int, int> > dims;
    vector< vector< vector<Id> > > incidence_lists;

    string check_filename = string ( datadir ) + string ( "hexa_cell.con" );
    read_connectivities ( check_filename.c_str ( ), dims, incidence_lists );

    for ( unsigned int i = 0; i < incidence_lists.size ( ); ++i )
    {
        check_incidence_iteration ( dims[i].first, dims[i].second, *mesh, incidence_lists[i] );
    }

    delete mb;

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );
    return 0;
}
