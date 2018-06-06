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

/// \author Martin Baumann, Thomas Gengenbach, Staffan Ronnas

#include <string>
#include <mpi.h>
#include <math.h>

#include "test.h"
#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

/// Extended DegreeOfFreedom Identification test
///
/// \brief For a list of files it is analysed, whether the DoF points given
///        by the initial numbering (discontinuous numbering) have common DoF
///        IDs for DoFs that are located at common points when a continuous
///        DoF numbering is calculated.
///
/// 1) Use discontinuous numbering and ...
/// 2) Iterate over all cells 'ca'
/// 3) -> Iterate over all DoF points 'pa' (local point index) in cell 'ca'
/// 4) -> -> Iterate over all cells 'cb'
/// 5) -> -> -> check wheter any point in cell 'cb' has same coordinates as
///             point 'ca', if so safe information
///                 (ca, pa, cb, pb)
/// 6) Use continous numbering, and check whether all DoF IDs that
///    correspond to (ca, pa) and (cb, pb) match.

static const char* datadir = MESH_DATADIR;

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );

    // Which files should be checked?

    std::vector<std::string> filenames;
    std::vector<TDim> tdims;
    std::vector<GDim> gdims;

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_triangles_2d.inp" ) );
    tdims.push_back ( 2 );
    gdims.push_back ( 2 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_quads_2d.inp" ) );
    tdims.push_back ( 2 );
    gdims.push_back ( 2 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_tetras_3d.inp" ) );
    tdims.push_back ( 3 );
    gdims.push_back ( 3 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_hexas_3d.inp" ) );
    tdims.push_back ( 3 );
    gdims.push_back ( 3 );

    for ( int test_number = 0; test_number<static_cast < int > ( filenames.size ( ) ); ++test_number )
    {

        std::string filename = filenames.at ( test_number );
        TDim tdim = tdims.at ( test_number );
        GDim gdim = gdims.at ( test_number );

        /////////////////////////////////////
        // mesh

        MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
        ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );
        MeshPtr org_mesh;
        reader->read ( filename.c_str ( ), org_mesh );

        MeshPtr refined_mesh;
        refined_mesh = org_mesh->refine ( );

        MeshPtr mesh = refined_mesh;

        /////////////////////////////////////
        // finite element manager

        int nvar = 1;
        int degree = 1;

        std::vector<int> data;
        data.push_back ( degree );

        // continuous fe manager
        FEManager<double>* fe_manager_c = new FEManager<double>( gdim, nvar );
        fe_manager_c->set_mesh ( *mesh );
        fe_manager_c->init_fe_tank ( 0, FEType<double>::LAGRANGE, data );

        // discontinuous fe manager
        FEManager<double>* fe_manager_dc = new FEManager<double>( gdim, nvar );
        fe_manager_dc->set_mesh ( *mesh );
        fe_manager_dc->init_fe_tank ( 0, FEType<double>::LAGRANGE, data );
        fe_manager_dc->set_ca ( 0, false ); // var 0 should be discontinuous

        /////////////////////////////////////
        // degrees of freedom

        // continuous dof
        DegreeOfFreedom<double>* dof_c = new DegreeOfFreedom<double>;
        dof_c->set_mesh ( mesh.get ( ) );
        dof_c->set_fe_manager ( fe_manager_c );
        dof_c->number ( );

        // discontinuous dof
        DegreeOfFreedom<double>* dof_dc = new DegreeOfFreedom<double>;
        dof_dc->set_mesh ( mesh.get ( ) );
        dof_dc->set_fe_manager ( fe_manager_dc );
        dof_dc->number ( );

        /////////////////////////////////////
        // testing

        // information to store:
        // for given index i, the following dofs correspond
        //      (cell_a, local_dof_a) <--> (cell_b, local_dof_b)

        std::vector<int> cell_a_vec;
        std::vector<int> cell_b_vec;
        std::vector<int> local_dof_a_vec;
        std::vector<int> local_dof_b_vec;

        // 1. calculate correspondence information based on a geometrical analysis
        //    using a numbering in discontinous mode

        // loop over cells

        for ( EntityIterator it = mesh->begin ( gdim );
              it != mesh->end ( gdim );
              ++it )
        {
            // get coordinates and DoF IDs

            std::vector<std::vector<double> > coords;
            dof_dc->get_coord_on_cell ( 0, it->index ( ), coords );

            std::vector<int> ids;
            dof_dc->get_dofs_on_cell ( 0, it->index ( ), ids );

            interminable_assert ( coords.size ( ) == ids.size ( ) );

            // iterate over DoF points

            for ( int dof_point = 0; dof_point<static_cast < int > ( coords.size ( ) ); ++dof_point )
            {

                std::vector<double> the_coords = coords.at ( dof_point );

                // iterate over other cells and look for corresponding DoF points

                for ( EntityIterator it2 = mesh->begin ( gdim );
                      it2 != mesh->end ( gdim );
                      ++it2 )
                {
                    // don't treat same cell

                    if ( it->index ( ) != it2->index ( ) )
                    {

                        // get coordinates and DoF IDs

                        std::vector<std::vector<double> > coords2;
                        dof_dc->get_coord_on_cell ( 0, it2->index ( ), coords2 );

                        std::vector<int> ids2;
                        dof_dc->get_dofs_on_cell ( 0, it2->index ( ), ids2 );

                        // look for corresponding DoF point in cell it2 as point 'the_coordinates'

                        int dof_point2 = 0;

                        while ( dof_point2 < static_cast < int > ( ids2.size ( ) ) )
                        {

                            double eps = 1.0e-12;

                            bool corresponding = true;
                            for ( int comp = 0; comp<static_cast < int > ( the_coords.size ( ) ); ++comp )
                                if ( fabs ( the_coords.at ( comp ) - coords2.at ( dof_point2 ).at ( comp ) ) > eps )
                                    corresponding = false;

                            if ( corresponding == true )
                            {
                                int cell_a = it->index ( );
                                int local_dof_a = dof_point;
                                int cell_b = it2->index ( );
                                int local_dof_b = dof_point2;

                                cell_a_vec.push_back ( cell_a );
                                cell_b_vec.push_back ( cell_b );
                                local_dof_a_vec.push_back ( local_dof_a );
                                local_dof_b_vec.push_back ( local_dof_b );

                                //std::cout << "FOUND IDENTICAL DOFs:" << std::endl;
                                //std::cout << "  (" << cell_a << ", " << local_dof_a << ") <->"
                                //          <<  " (" << cell_b << ", " << local_dof_b << ") \t //";
                                //for (int comp=0; comp<the_coords.size(); ++comp)
                                //  std::cout << " " << the_coords.at(comp) << " /";
                                //std::cout << "/" << std::endl;

                                break;
                            }

                            ++dof_point2;

                        } // while (dof_point2 < ids2.size())

                    } // if (it->index() != it2->index())

                } // for (MeshEntityIterator it2 = mesh->begin(gdim)

            } // for (int dof_point=0; dof_point<coords.size(); ++dof_point)

        } // for (EntityIterator it = mesh->begin(gdim);

        // 2. check that found identifyable DoFs have common DoF ID in
        //    the continuous DoF structure

        // cell_a_vec.size() must be even number
        interminable_assert ( cell_a_vec.size ( ) == 2 * ( cell_a_vec.size ( ) / 2 ) );

        std::cout << "# DoF tuples: " << cell_a_vec.size ( ) / 2 << std::endl;

        for ( int i = 0; i<static_cast < int > ( cell_a_vec.size ( ) ); ++i )
        {
            std::vector<int> dofs_a;
            dof_c->get_dofs_on_cell ( 0, cell_a_vec.at ( i ), dofs_a );
            int id_a = dofs_a.at ( local_dof_a_vec.at ( i ) );

            std::vector<int> dofs_b;
            dof_c->get_dofs_on_cell ( 0, cell_b_vec.at ( i ), dofs_b );
            int id_b = dofs_b.at ( local_dof_b_vec.at ( i ) );

            //interminable_assert(id_a == id_b);
            TEST_EQUAL ( id_a, id_b );
        }

        // TEST_EQUAL(4, n_dof);

        // free memory

        delete dof_c;
        delete dof_dc;

        delete fe_manager_c;
        delete fe_manager_dc;

        delete mb;
    } // for (int test_number=0; test_number<filenames.size(); ++test_number)

    MPI_Finalize ( );

    return 0;

}
