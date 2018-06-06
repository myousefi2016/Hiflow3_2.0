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

#include "fecellwisemanager.h"
#include <cassert>
#include "mesh/iterator.h"
#include "common/macros.h"

#include "bilinearquadtransformation.h"
#include "lineartriangletransformation.h"
#include "lineartetrahedrontransformation.h"
#include "trilinearhexahedrontransformation.h"
#include "linearpyramidtransformation.h"
#include "felagrange.h"
#include "felagrange_tri.h"
#include "felagrange_quad.h"
#include "felagrange_tet.h"
#include "felagrange_hex.h"
#include "felagrange_pyr.h"

/// \author Michael Schick<br>Martin Baumann

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FECellwiseManager<DataType>::FECellwiseManager ( int dim, int nb_var ) : FEManager<DataType>( dim, nb_var )
        {
        }

        template<class DataType>
        void FECellwiseManager<DataType>::reinit_fe_tank ( const std::vector<std::vector<int> >& fe_data )
        {
            assert ( 0 );
            //  int nb_awake_cells = mesh_->num_entities(dim_);

            //  fe_inst_.clear();
            //  fe_tank_.clear();
            //  fe_tank_.resize(nb_var_*nb_awake_cells);

            //  assert(fe_tank_.size() == fe_data.size());

            //  for(int var = 0; var < nb_var_; ++var)
            //  {
            //    for (MeshEntityIterator it = mesh_->begin(dim_); it != mesh_->end(dim_); ++it)
            //    {
            //      // 1st: Get Cell Type (Tet, Hex, Quad ...)
            //      CellType::Type cell_type = mesh::get_cell_type(dim_, it->num_vertices())->cell_type();

            //      // 2nd: Store Finite Element in FE Tank

            //      switch(cell_type) {

            //      case CellType::TRIANGLE : {
            //                if(fe_data[it->index() + var*nb_awake_cells][0] == FEType::LAGRANGE)
            //                  fe_tank_.at(it->index() + var*nb_awake_cells) =
            //                    fe_inst_.get_lagrange_tri_element(fe_data[it->index() + var*nb_awake_cells][1]);
            //                else
            //                  assert(0);
            //                break; }
            //
            //      case CellType::QUADRILATERAL : {
            //                if(fe_data[it->index() + var*nb_awake_cells][0] == FEType::LAGRANGE)
            //                  fe_tank_.at(it->index() + var*nb_awake_cells) =
            //                    fe_inst_.get_lagrange_quad_element(fe_data[it->index() + var*nb_awake_cells][1]);
            //                else
            //                  assert(0);
            //                break; }

            //      case CellType::TETRAHEDRON : {
            //                if(fe_data[it->index() + var*nb_awake_cells][0] == FEType::LAGRANGE)
            //                  fe_tank_.at(it->index() + var*nb_awake_cells) =
            //                    fe_inst_.get_lagrange_tet_element(fe_data[it->index() + var*nb_awake_cells][1]);
            //                else
            //                  assert(0);
            //                break; }

            //      case CellType::HEXAHEDRON : {
            //                if(fe_data[it->index() + var*nb_awake_cells][0] == FEType::LAGRANGE)
            //                {
            //                  fe_tank_.at(it->index() + var*nb_awake_cells) =
            //                    fe_inst_.get_lagrange_hex_element(fe_data[it->index() + var*nb_awake_cells][1]);
            //                }
            //                else
            //                  assert(0);
            //                break; }

            //        default: assert(0);
            //      }
            //    }
            //  }
        }

        template<class DataType>
        void FECellwiseManager<DataType>::init_fe_tank ( int var,
                                                         const typename FEType<DataType>::FEAnsatz &ansatz,
                                                         const std::vector<int>& param )
        {
            assert ( !( this->initialized_[var] ) );

            // Initialize pointers to existing entity types
            mesh::EntityNumber tri_representer = -1;
            mesh::EntityNumber quad_representer = -1;
            mesh::EntityNumber tet_representer = -1;
            mesh::EntityNumber hex_representer = -1;
            mesh::EntityNumber pyr_representer = -1;

            const size_t offset = var * this->num_entities_;

            for ( mesh::EntityIterator it = this->mesh_->begin ( this->dim_ ), e_it = this->mesh_->end ( this->dim_ ); it != e_it; ++it )
            {
                // 1st: Get Cell Type (Tet, Hex, Quad ...)
                mesh::CellType::Tag cell_type = it->cell_type ( ).tag ( );

                int cell_index = it->index ( );
                int lagrange_degree = param[cell_index];

                // 2nd: Store Finite Element in FE Tank
                switch ( cell_type )
                {

                    case mesh::CellType::TRIANGLE:
                    {
                        if ( ansatz == FEType<DataType>::LAGRANGE )
                        {
                            tri_representer = it->index ( );
                            this->fe_tank_[it->index ( ) + offset] =
                                    this->fe_inst_.get_lagrange_tri_element ( lagrange_degree );
                        }
                        else
                            assert ( 0 );
                        break;
                    }

                    case mesh::CellType::QUADRILATERAL:
                    {
                        if ( ansatz == FEType<DataType>::LAGRANGE )
                        {
                            quad_representer = it->index ( );
                            this->fe_tank_[it->index ( ) + offset] =
                                    this->fe_inst_.get_lagrange_quad_element ( lagrange_degree );
                        }
                        else
                            assert ( 0 );
                        break;
                    }

                    case mesh::CellType::TETRAHEDRON:
                    {
                        if ( ansatz == FEType<DataType>::LAGRANGE )
                        {
                            tet_representer = it->index ( );
                            this->fe_tank_[it->index ( ) + offset] =
                                    this->fe_inst_.get_lagrange_tet_element ( lagrange_degree );
                        }
                        else
                            assert ( 0 );
                        break;
                    }

                    case mesh::CellType::HEXAHEDRON:
                    {
                        if ( ansatz == FEType<DataType>::LAGRANGE )
                        {
                            hex_representer = it->index ( );
                            this->fe_tank_[it->index ( ) + offset] =
                                    this->fe_inst_.get_lagrange_hex_element ( lagrange_degree );
                        }
                        else
                            assert ( 0 );
                        break;
                    }

                    case mesh::CellType::PYRAMID:
                    {
                        if ( ansatz == FEType<DataType>::LAGRANGE )
                        {
                            pyr_representer = it->index ( );
                            this->fe_tank_[it->index ( ) + offset] =
                                    this->fe_inst_.get_lagrange_pyr_element ( lagrange_degree );
                        }
                        else
                            assert ( 0 );
                        break;
                    }

                    default: assert ( 0 );
                }
            }

            // At least one cell should have been found
            assert ( !( tri_representer == -1 && quad_representer == -1
                     && tet_representer == -1 && hex_representer == -1
                     && pyr_representer == -1 ) );

            // Last:  setting up numbering of reference cell via transformation of an arbitrary
            //        (here the last) mesh entity back to reference cell and iteration
            //        through the vertex ids

            if ( tri_representer != -1 )
            {
                std::vector<DataType> coordinates;
                mesh::Entity entity = this->mesh_->get_entity ( this->dim_, tri_representer );
                entity.get_coordinates ( coordinates );

                LinearTriangleTransformation<DataType> ltt ( 2 );
                ltt.reinit ( coordinates );

                std::vector<Coord> ref_coords ( 3 );
                for ( int i = 0; i < 3; ++i )
                    ref_coords[i].resize ( 2 );

                for ( int i = 0; i < 3; ++i )
                    ltt.inverse ( coordinates[2 * i], coordinates[2 * i + 1], ref_coords[i][0], ref_coords[i][1] );

                this->fe_inst_.init_triangle_elements ( ref_coords );
            }

            if ( quad_representer != -1 )
            {
                std::vector<DataType> coordinates;
                mesh::Entity entity = this->mesh_->get_entity ( this->dim_, quad_representer );
                entity.get_coordinates ( coordinates );

                BiLinearQuadTransformation<DataType> blqt ( 2 );
                blqt.reinit ( coordinates );

                std::vector<Coord> ref_coords ( 4 );
                for ( int i = 0; i < 4; ++i )
                    ref_coords[i].resize ( 2 );

                for ( int i = 0; i < 4; ++i )
                    blqt.inverse ( coordinates[2 * i], coordinates[2 * i + 1], ref_coords[i][0], ref_coords[i][1] );

                this->fe_inst_.init_quadrilateral_elements ( ref_coords );
            }

            if ( tet_representer != -1 )
            {
                std::vector<DataType> coordinates;
                mesh::Entity entity = this->mesh_->get_entity ( this->dim_, tet_representer );
                entity.get_coordinates ( coordinates );

                LinearTetrahedronTransformation<DataType> ltt ( 3 );
                ltt.reinit ( coordinates );

                std::vector<Coord> ref_coords ( 4 );
                for ( int i = 0; i < 4; ++i )
                    ref_coords[i].resize ( 3 );

                for ( int i = 0; i < 4; ++i )
                    ltt.inverse ( coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2],
                                  ref_coords[i][0], ref_coords[i][1], ref_coords[i][2] );

                this->fe_inst_.init_tetrahedron_elements ( ref_coords );
            }

            if ( hex_representer != -1 )
            {
                std::vector<DataType> coordinates;
                mesh::Entity entity = this->mesh_->get_entity ( this->dim_, hex_representer );
                entity.get_coordinates ( coordinates );

                TriLinearHexahedronTransformation<DataType> tlht ( 3 );
                tlht.reinit ( coordinates );

                std::vector<Coord> ref_coords ( 8 );
                for ( int i = 0; i < 8; ++i )
                    ref_coords[i].resize ( 3 );

                for ( int i = 0; i < 8; ++i )
                    tlht.inverse ( coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2],
                                   ref_coords[i][0], ref_coords[i][1], ref_coords[i][2] );

                this->fe_inst_.init_hexahedron_elements ( ref_coords );
            }

            if ( pyr_representer != -1 )
            {
                std::vector<DataType> coordinates;
                mesh::Entity entity = this->mesh_->get_entity ( this->dim_, pyr_representer );
                entity.get_coordinates ( coordinates );

                LinearPyramidTransformation<DataType> lpr ( 3 );
                lpr.reinit ( coordinates );

                std::vector<Coord> ref_coords ( 5 );
                for ( int i = 0; i < 5; ++i )
                    ref_coords[i].resize ( 3 );

                for ( int i = 0; i < 5; ++i )
                    lpr.inverse ( coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2],
                                  ref_coords[i][0], ref_coords[i][1], ref_coords[i][2] );

                this->fe_inst_.init_pyramid_elements ( ref_coords );
            }

            this->initialized_[var] = true;
        }

        template class FECellwiseManager<double>;
        template class FECellwiseManager<float>;

    }
} // namespace hiflow
