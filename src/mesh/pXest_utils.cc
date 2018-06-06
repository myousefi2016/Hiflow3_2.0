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

#include "pXest_utils.h"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include "common/log.h"

#include "communication.h"
#include "iterator.h"
#include "mpi_communication.h"
#include "interface.h"

#ifdef WITH_P4EST
#    include "p4est_algorithms.h"
#    include "p8est_algorithms.h"
#endif

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        void print_pack ( EntityPackage& pack )
        {
            std::cout << " EntPack " << pack.tdim << " " << pack.gdim << std::endl;
            std::cout << " connections " << std::endl;
            for ( int l = 0; l < pack.connections.size ( ); ++l )
            {
                std::cout << pack.connections[l] << " ";
            }
            std::cout << std::endl << " coords " << std::endl;
            for ( int l = 0; l < pack.coords.size ( ); ++l )
            {
                std::cout << pack.coords[l] << " ";
            }
            std::cout << std::endl << " -------- " << std::endl;
        }

        void pack_ghost_data ( GhostCommData* ret, QuadData& quad_data )
        {
            ret->cell_index = quad_data.cell_index;
        }

        void unpack_ghost_data ( GhostCommData& ghost_data, QuadData& quad_data )
        {
            quad_data.cell_index = ghost_data.cell_index;
        }

        void print_ghost_data ( GhostCommData22* ret )
        {
            std::cout << " cell id ";
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                std::cout << ret->cell_id[l] << " ";
            }
            std::cout << std::endl;
            std::cout << " cell id remote ";
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                std::cout << ret->cell_id_remote[l] << " ";
            }
            std::cout << std::endl;
            std::cout << " cell index " << ret->cell_index << std::endl;
            std::cout << " cell index remote " << ret->cell_index_remote << std::endl;
            std::cout << " tree id " << ret->tree_id << std::endl;
            std::cout << " tree id remote " << ret->tree_id_remote << std::endl;
            std::cout << " owner " << ret->owner_rank << std::endl;
            std::cout << " active " << ret->active << std::endl;
            std::cout << " refine " << ret->refine << std::endl;
            std::cout << " coarsen " << ret->coarsen << std::endl;
            std::cout << " tdim " << ret->tdim << std::endl;
            std::cout << " gdim " << ret->gdim << std::endl;

            std::cout << std::endl << " offsets " << std::endl;
            for ( int l = 0; l < 10; ++l )
            {
                std::cout << ret->offsets[l] << " ";
            }

            std::cout << std::endl << " conn_offsets " << std::endl;
            for ( int l = 0; l < 6; ++l )
            {
                std::cout << ret->conn_offsets[l] << " ";
            }

            std::cout << std::endl << " conn_offsets " << std::endl;
            for ( int l = 0; l < 6; ++l )
            {
                std::cout << ret->coord_offsets[l] << " ";
            }

            std::cout << std::endl << " connections " << std::endl;
            for ( int l = 0; l < 12; ++l )
            {
                std::cout << ret->connections[l] << " ";
            }

            std::cout << std::endl << " coords " << std::endl;
            for ( int l = 0; l < 24; ++l )
            {
                std::cout << ret->coords[l] << " ";
            }
            std::cout << std::endl << " ###### " << std::endl;
            ;
        }

        void print_ghost_data ( GhostCommData33* ret )
        {
            std::cout << " cell id ";
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                std::cout << ret->cell_id[l] << " ";
            }
            std::cout << std::endl;
            std::cout << " cell id remote ";
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                std::cout << ret->cell_id_remote[l] << " ";
            }
            std::cout << std::endl;
            std::cout << " cell index " << ret->cell_index << std::endl;
            std::cout << " cell index remote " << ret->cell_index_remote << std::endl;
            std::cout << " tree id " << ret->tree_id << std::endl;
            std::cout << " tree id remote " << ret->tree_id_remote << std::endl;
            std::cout << " owner " << ret->owner_rank << std::endl;
            std::cout << " active " << ret->active << std::endl;
            std::cout << " refine " << ret->refine << std::endl;
            std::cout << " coarsen " << ret->coarsen << std::endl;
            std::cout << " tdim " << ret->tdim << std::endl;
            std::cout << " gdim " << ret->gdim << std::endl;

            std::cout << std::endl << " offsets " << std::endl;
            for ( int l = 0; l < 14; ++l )
            {
                std::cout << ret->offsets[l] << " ";
            }

            std::cout << std::endl << " conn_offsets " << std::endl;
            for ( int l = 0; l < 8; ++l )
            {
                std::cout << ret->conn_offsets[l] << " ";
            }

            std::cout << std::endl << " coord_offsets " << std::endl;
            for ( int l = 0; l < 8; ++l )
            {
                std::cout << ret->coord_offsets[l] << " ";
            }

            std::cout << std::endl << " connections " << std::endl;
            for ( int l = 0; l < 32; ++l )
            {
                std::cout << ret->connections[l] << " ";
            }

            std::cout << std::endl << " coords " << std::endl;
            for ( int l = 0; l < 96; ++l )
            {
                std::cout << ret->coords[l] << " ";
            }
            std::cout << std::endl << " ###### " << std::endl;
            ;

        }

        void pack_ghost_data ( GhostCommData22* ret, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs )
        {
            const int num_vert_cell = 4;
            const int num_vert_facet = 2;
            const int gdim = 2;

            const int max_num_cell = P4EST_MAXLEVEL;
            const int max_num_facet = P4EST_MAXLEVEL * 4;
            assert ( cell_packs.size ( ) <= max_num_cell );
            assert ( facet_packs.size ( ) <= max_num_facet );

            ret->num_cell = cell_packs.size ( );
            ret->num_facet = facet_packs.size ( );
            ret->num_vert_cell = num_vert_cell;
            ret->num_vert_facet = num_vert_facet;
            ret->cell_index = quad_data.cell_index;
            ret->cell_index_remote = quad_data.cell_index_remote;
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                ret->cell_id[l] = quad_data.cell_id[l];
                ret->cell_id_remote[l] = quad_data.cell_id_remote[l];
            }
            ret->tree_id = quad_data.tree_id;
            ret->tree_id_remote = quad_data.tree_id_remote;
            ret->owner_rank = quad_data.owner_rank;
            ret->tdim = quad_data.tdim;
            ret->active = quad_data.active;
            ret->refine = quad_data.refine;
            ret->coarsen = quad_data.coarsen;
            ret->gdim = gdim;
            ret->conn_offsets[0] = 0;
            ret->coord_offsets[0] = 0;

            for ( int jc = 0; jc < ret->num_cell; ++jc )
            {
                for ( int jv = 0; jv < num_vert_cell; ++jv )
                {
                    ret->connections[ ret->conn_offsets[jc] + jv] = cell_packs[jc].connections[jv];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        ret->coords[ ret->coord_offsets[jc] + jv * gdim + jd ] = cell_packs[jc].coords[jv * gdim + jd];
                    }
                }
                ret->tdims[jc] = cell_packs[jc].tdim;
                ret->offsets[jc * 2] = cell_packs[jc].offsets[0];
                ret->offsets[jc * 2 + 1] = cell_packs[jc].offsets[1];
                ret->conn_offsets[jc + 1] = ret->conn_offsets[jc] + num_vert_cell;
                ret->coord_offsets[jc + 1] = ret->coord_offsets[jc] + num_vert_cell*gdim;
                ret->material_numbers[jc] = cell_packs[jc].material_numbers[0];
            }
            for ( int jc = 0; jc < ret->num_facet; ++jc )
            {
                for ( int jv = 0; jv < num_vert_facet; ++jv )
                {
                    ret->connections[ ret->conn_offsets[ret->num_cell + jc] + jv ] = facet_packs[jc].connections[jv];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        ret->coords[ ret->coord_offsets[ret->num_cell + jc] + jv * gdim + jd ] = facet_packs[jc].coords[jv * gdim + jd];
                    }
                }
                ret->tdims[ret->num_cell + jc] = facet_packs[jc].tdim;
                ret->offsets[ret->num_cell * 2 + jc * 2] = facet_packs[jc].offsets[0];
                ret->offsets[ret->num_cell * 2 + jc * 2 + 1] = facet_packs[jc].offsets[1];
                ret->conn_offsets[ret->num_cell + jc + 1] = ret->conn_offsets[ret->num_cell + jc] + num_vert_facet;
                ret->coord_offsets[ret->num_cell + jc + 1] = ret->coord_offsets[ret->num_cell + jc] + num_vert_facet*gdim;
                ret->material_numbers[ret->num_cell + jc] = facet_packs[jc].material_numbers[0];
            }
        }

        void pack_ghost_data ( GhostCommData33* ret, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs )
        {
            const int num_vert_cell = 8;
            const int num_vert_facet = 4;
            const int gdim = 3;

            const int max_num_cell = P4EST_MAXLEVEL;
            const int max_num_facet = P4EST_MAXLEVEL * 6;
            assert ( cell_packs.size ( ) <= max_num_cell );
            assert ( facet_packs.size ( ) <= max_num_facet );

            ret->num_cell = cell_packs.size ( );
            ret->num_facet = facet_packs.size ( );
            ret->num_vert_cell = num_vert_cell;
            ret->num_vert_facet = num_vert_facet;
            ret->cell_index = quad_data.cell_index;
            ret->cell_index_remote = quad_data.cell_index_remote;
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                ret->cell_id[l] = quad_data.cell_id[l];
                ret->cell_id_remote[l] = quad_data.cell_id_remote[l];
            }
            ret->tree_id = quad_data.tree_id;
            ret->tree_id_remote = quad_data.tree_id_remote;
            ret->owner_rank = quad_data.owner_rank;
            ret->tdim = quad_data.tdim;
            ret->active = quad_data.active;
            ret->refine = quad_data.refine;
            ret->coarsen = quad_data.coarsen;
            ret->gdim = gdim;
            ret->conn_offsets[0] = 0;
            ret->coord_offsets[0] = 0;

            for ( int jc = 0; jc < ret->num_cell; ++jc )
            {
                for ( int jv = 0; jv < num_vert_cell; ++jv )
                {
                    ret->connections[ ret->conn_offsets[jc] + jv] = cell_packs[jc].connections[jv];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        ret->coords[ ret->coord_offsets[jc] + jv * gdim + jd ] = cell_packs[jc].coords[jv * gdim + jd];
                    }
                }
                ret->tdims[jc] = cell_packs[jc].tdim;
                ret->offsets[jc * 2] = cell_packs[jc].offsets[0];
                ret->offsets[jc * 2 + 1] = cell_packs[jc].offsets[1];
                ret->conn_offsets[jc + 1] = ret->conn_offsets[jc] + num_vert_cell;
                ret->coord_offsets[jc + 1] = ret->coord_offsets[jc] + num_vert_cell*gdim;
                ret->material_numbers[jc] = cell_packs[jc].material_numbers[0];
            }
            for ( int jc = 0; jc < ret->num_facet; ++jc )
            {
                for ( int jv = 0; jv < num_vert_facet; ++jv )
                {
                    ret->connections[ ret->conn_offsets[ret->num_cell + jc] + jv ] = facet_packs[jc].connections[jv];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        ret->coords[ ret->coord_offsets[ret->num_cell + jc] + jv * gdim + jd ] = facet_packs[jc].coords[jv * gdim + jd];
                    }
                }
                ret->tdims[ret->num_cell + jc] = facet_packs[jc].tdim;
                ret->offsets[ret->num_cell * 2 + jc * 2] = facet_packs[jc].offsets[0];
                ret->offsets[ret->num_cell * 2 + jc * 2 + 1] = facet_packs[jc].offsets[1];
                ret->conn_offsets[ret->num_cell + jc + 1] = ret->conn_offsets[ret->num_cell + jc] + num_vert_facet;
                ret->coord_offsets[ret->num_cell + jc + 1] = ret->coord_offsets[ret->num_cell + jc] + num_vert_facet*gdim;
                ret->material_numbers[ret->num_cell + jc] = facet_packs[jc].material_numbers[0];
            }
        }

        // BUG_TODO material numbers

        void pack_ghost_data ( GhostCommData22* ret, QuadData& quad_data )
        {
            const int num_vert_cell = 4;
            const int num_vert_facet = 2;
            const int num_cell = 0;
            const int num_facet = 0;
            const int gdim = 2;

            ret->num_cell = num_cell;
            ret->num_facet = num_facet;
            ret->num_vert_cell = num_vert_cell;
            ret->num_vert_facet = num_vert_facet;
            ret->cell_index = quad_data.cell_index;
            ret->cell_index_remote = quad_data.cell_index_remote;
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                ret->cell_id[l] = quad_data.cell_id[l];
                ret->cell_id_remote[l] = quad_data.cell_id_remote[l];
            }
            ret->tree_id = quad_data.tree_id;
            ret->tree_id_remote = quad_data.tree_id_remote;
            ret->owner_rank = quad_data.owner_rank;
            ret->tdim = quad_data.tdim;
            ret->active = quad_data.active;
            ret->refine = quad_data.refine;
            ret->coarsen = quad_data.coarsen;
            ret->gdim = gdim;
        }

        void pack_ghost_data ( GhostCommData33* ret, QuadData& quad_data )
        {
            const int num_vert_cell = 8;
            const int num_vert_facet = 4;
            const int num_cell = 0;
            const int num_facet = 0;
            const int gdim = 3;

            ret->num_cell = num_cell;
            ret->num_facet = num_facet;
            ret->num_vert_cell = num_vert_cell;
            ret->num_vert_facet = num_vert_facet;
            ret->cell_index = quad_data.cell_index;
            ret->cell_index_remote = quad_data.cell_index_remote;
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                ret->cell_id[l] = quad_data.cell_id[l];
                ret->cell_id_remote[l] = quad_data.cell_id_remote[l];
            }
            ret->tree_id = quad_data.tree_id;
            ret->tree_id_remote = quad_data.tree_id_remote;
            ret->owner_rank = quad_data.owner_rank;
            ret->tdim = quad_data.tdim;
            ret->active = quad_data.active;
            ret->refine = quad_data.refine;
            ret->coarsen = quad_data.coarsen;
            ret->gdim = gdim;
        }

        void unpack_ghost_data ( GhostCommData22& ghost_data, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs )
        {
            const int num_cell = ghost_data.num_cell;
            const int num_facet = ghost_data.num_facet;
            const int num_vert_cell = ghost_data.num_vert_cell;
            const int num_vert_facet = ghost_data.num_vert_facet;
            const int gdim = ghost_data.gdim;

            assert ( gdim == 2 );

            cell_packs.resize ( num_cell );
            facet_packs.resize ( num_facet );

            quad_data.cell_index = ghost_data.cell_index;
            quad_data.cell_index_remote = ghost_data.cell_index_remote;
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                quad_data.cell_id[l] = ghost_data.cell_id[l];
                quad_data.cell_id_remote[l] = ghost_data.cell_id_remote[l];
            }
            quad_data.tree_id = ghost_data.tree_id;
            quad_data.tree_id_remote = ghost_data.tree_id_remote;
            quad_data.owner_rank = ghost_data.owner_rank;
            quad_data.tdim = ghost_data.tdim;
            quad_data.active = ghost_data.active;
            quad_data.refine = ghost_data.refine;
            quad_data.coarsen = ghost_data.coarsen;

            for ( int jc = 0; jc < num_cell; ++jc )
            {
                cell_packs[jc].connections.resize ( num_vert_cell, 0 );
                cell_packs[jc].coords.resize ( num_vert_cell*gdim, 0. );
                cell_packs[jc].offsets.resize ( 2, 0 );
                cell_packs[jc].material_numbers.resize ( 1, 0 );
                cell_packs[jc].offsets[0] = ghost_data.offsets[jc * 2 + 0];
                cell_packs[jc].offsets[1] = ghost_data.offsets[jc * 2 + 1];
                cell_packs[jc].tdim = ghost_data.tdims[jc];
                cell_packs[jc].gdim = ghost_data.gdim;
                cell_packs[jc].material_numbers[0] = ghost_data.material_numbers[jc];

                for ( int jv = 0; jv < num_vert_cell; ++jv )
                {
                    cell_packs[jc].connections[jv] = ghost_data.connections[ ghost_data.conn_offsets[jc] + jv ];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        cell_packs[jc].coords[jv * gdim + jd] = ghost_data.coords[ ghost_data.coord_offsets[jc] + jv * gdim + jd];
                    }
                }
            }
            for ( int jc = 0; jc < num_facet; ++jc )
            {
                facet_packs[jc].connections.resize ( num_vert_facet, 0 );
                facet_packs[jc].coords.resize ( num_vert_facet*gdim, 0. );
                facet_packs[jc].offsets.resize ( 2, 0 );
                facet_packs[jc].material_numbers.resize ( 1, 0 );
                facet_packs[jc].offsets[0] = ghost_data.offsets[num_cell * 2 + jc * 2 + 0];
                facet_packs[jc].offsets[1] = ghost_data.offsets[num_cell * 2 + jc * 2 + 1];
                facet_packs[jc].tdim = ghost_data.tdims[num_cell + jc];
                facet_packs[jc].gdim = ghost_data.gdim;
                facet_packs[jc].material_numbers[0] = ghost_data.material_numbers[num_cell + jc];
                for ( int jv = 0; jv < num_vert_facet; ++jv )
                {
                    facet_packs[jc].connections[jv] = ghost_data.connections[ ghost_data.conn_offsets[num_cell + jc] + jv ];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        facet_packs[jc].coords[jv * gdim + jd] = ghost_data.coords[ ghost_data.coord_offsets[num_cell + jc] + jv * gdim + jd];
                    }
                }
            }
        }

        void unpack_ghost_data ( GhostCommData33& ghost_data, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs )
        {
            const int num_cell = ghost_data.num_cell;
            const int num_facet = ghost_data.num_facet;
            const int num_vert_cell = ghost_data.num_vert_cell;
            const int num_vert_facet = ghost_data.num_vert_facet;
            const int gdim = ghost_data.gdim;

            assert ( gdim == 3 );

            cell_packs.resize ( num_cell );
            facet_packs.resize ( num_facet );

            quad_data.cell_index = ghost_data.cell_index;
            quad_data.cell_index_remote = ghost_data.cell_index_remote;
            for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
            {
                quad_data.cell_id[l] = ghost_data.cell_id[l];
                quad_data.cell_id_remote[l] = ghost_data.cell_id_remote[l];
            }
            quad_data.tree_id = ghost_data.tree_id;
            quad_data.tree_id_remote = ghost_data.tree_id_remote;
            quad_data.owner_rank = ghost_data.owner_rank;
            quad_data.tdim = ghost_data.tdim;
            quad_data.active = ghost_data.active;
            quad_data.refine = ghost_data.refine;
            quad_data.coarsen = ghost_data.coarsen;

            for ( int jc = 0; jc < num_cell; ++jc )
            {
                cell_packs[jc].connections.resize ( num_vert_cell, 0 );
                cell_packs[jc].coords.resize ( num_vert_cell*gdim, 0. );
                cell_packs[jc].offsets.resize ( 2, 0 );
                cell_packs[jc].material_numbers.resize ( 1, 0 );
                cell_packs[jc].offsets[0] = ghost_data.offsets[jc * 2 + 0];
                cell_packs[jc].offsets[1] = ghost_data.offsets[jc * 2 + 1];
                cell_packs[jc].tdim = ghost_data.tdims[jc];
                cell_packs[jc].gdim = ghost_data.gdim;
                cell_packs[jc].material_numbers[0] = ghost_data.material_numbers[jc];

                for ( int jv = 0; jv < num_vert_cell; ++jv )
                {
                    cell_packs[jc].connections[jv] = ghost_data.connections[ ghost_data.conn_offsets[jc] + jv ];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        cell_packs[jc].coords[jv * gdim + jd] = ghost_data.coords[ ghost_data.coord_offsets[jc] + jv * gdim + jd];
                    }
                }
            }
            for ( int jc = 0; jc < num_facet; ++jc )
            {
                facet_packs[jc].connections.resize ( num_vert_facet, 0 );
                facet_packs[jc].coords.resize ( num_vert_facet*gdim, 0. );
                facet_packs[jc].offsets.resize ( 2, 0 );
                facet_packs[jc].material_numbers.resize ( 1, 0 );
                facet_packs[jc].offsets[0] = ghost_data.offsets[num_cell * 2 + jc * 2 + 0];
                facet_packs[jc].offsets[1] = ghost_data.offsets[num_cell * 2 + jc * 2 + 1];
                facet_packs[jc].tdim = ghost_data.tdims[num_cell + jc];
                facet_packs[jc].gdim = ghost_data.gdim;
                facet_packs[jc].material_numbers[0] = ghost_data.material_numbers[num_cell + jc];
                for ( int jv = 0; jv < num_vert_facet; ++jv )
                {
                    facet_packs[jc].connections[jv] = ghost_data.connections[ ghost_data.conn_offsets[num_cell + jc] + jv ];
                    for ( int jd = 0; jd < gdim; ++jd )
                    {
                        facet_packs[jc].coords[jv * gdim + jd] = ghost_data.coords[ ghost_data.coord_offsets[num_cell + jc] + jv * gdim + jd];
                    }
                }
            }
        }

#ifdef WITH_P4EST

        void pXest_copy_ghost ( p4est_ghost_t* input, p4est_ghost_t*& ghost )
        {
            ghost = new p4est_ghost_t;
            int num_mirrors = input->mirrors.elem_count;
            int mpisize = input->mpisize;
            int num_trees = input->num_trees;

            // scalars
            ghost->mpisize = input->mpisize;
            ghost->num_trees = input->num_trees;
            ghost->btype = input->btype;

            // fields
            ghost->tree_offsets = new p4est_locidx_t [num_trees + 1];
            for ( int i = 0; i < num_trees + 1; ++i )
            {
                ghost->tree_offsets[i] = input->tree_offsets[i];
            }

            ghost->proc_offsets = new p4est_locidx_t [mpisize + 1];
            for ( int i = 0; i < mpisize + 1; ++i )
            {
                ghost->proc_offsets[i] = input->proc_offsets[i];
            }

            ghost->mirror_tree_offsets = new p4est_locidx_t [num_trees + 1];
            for ( int i = 0; i < num_trees + 1; ++i )
            {
                ghost->mirror_tree_offsets[i] = input->mirror_tree_offsets[i];
            }

            ghost->mirror_proc_offsets = new p4est_locidx_t [mpisize + 1];
            for ( int i = 0; i < mpisize + 1; ++i )
            {
                ghost->mirror_proc_offsets[i] = input->mirror_proc_offsets[i];
            }

            ghost->mirror_proc_mirrors = new p4est_locidx_t [ ghost->mirror_proc_offsets[mpisize] ];
            for ( int i = 0; i < ghost->mirror_proc_offsets[mpisize]; ++i )
            {
                ghost->mirror_proc_mirrors[i] = input->mirror_proc_mirrors[i];
            }

            ghost->mirror_proc_front_offsets = new p4est_locidx_t [mpisize + 1];
            for ( int i = 0; i < mpisize + 1; ++i )
            {
                ghost->mirror_proc_front_offsets[i] = input->mirror_proc_front_offsets[i];
            }

            ghost->mirror_proc_fronts = new p4est_locidx_t [ ghost->mirror_proc_front_offsets[mpisize] ];
            for ( int i = 0; i < ghost->mirror_proc_front_offsets[mpisize]; ++i )
            {
                ghost->mirror_proc_fronts[i] = input->mirror_proc_fronts[i];
            }

            // quadrants
            size_t elem_size = input->ghosts.elem_size;
            sc_array_init ( &ghost->ghosts, elem_size );
            sc_array_copy ( &ghost->ghosts, &input->ghosts );
            for ( int i = 0; i < ghost->ghosts.elem_count; ++i )
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( p4est_quadrant_array_index ( &input->ghosts, i ) );
                pXest_set_quad_data ( p4est_quadrant_array_index ( &ghost->ghosts, i ), *q_ptr );
            }

            elem_size = input->mirrors.elem_size;
            sc_array_init ( &ghost->mirrors, elem_size );
            sc_array_copy ( &ghost->mirrors, &input->mirrors );
            for ( int i = 0; i < ghost->mirrors.elem_count; ++i )
            {
                p4est_quadrant_t* new_quad = p4est_quadrant_array_index ( &ghost->mirrors, i );
                p4est_quadrant_t* old_quad = p4est_quadrant_array_index ( &input->mirrors, i );

                new_quad->p.piggy3.which_tree = old_quad->p.piggy3.which_tree;
                new_quad->p.piggy3.local_num = old_quad->p.piggy3.local_num;
            }

            /*
            std::cout << input->ghosts.elem_count << std::endl;
            p4est_quadrant_t* quad = p4est_quadrant_array_index(&input->ghosts, 1);
            QuadData* q_ptr = pXest_get_quad_data_ptr(quad);

            q_ptr->print();
            std::cout << (int) quad->level << std::endl;
            std::cout << quad->x << std::endl;
            std::cout << quad->p.which_tree << std::endl;
            std::cout << quad->p.piggy1.which_tree << " " << quad->p.piggy1.owner_rank << std::endl;
            std::cout << quad->p.piggy3.which_tree << " " << quad->p.piggy3.local_num << std::endl;
            std::cout << " ------------------------- " << std::endl;
             */
        }

        void pXest_copy_ghost ( p8est_ghost_t* input, p8est_ghost_t*& ghost )
        {
            ghost = new p8est_ghost_t;

            int num_mirrors = input->mirrors.elem_count;
            int mpisize = input->mpisize;
            int num_trees = input->num_trees;

            // scalars
            ghost->mpisize = input->mpisize;
            ghost->num_trees = input->num_trees;
            ghost->btype = input->btype;

            // fields
            ghost->tree_offsets = new p4est_locidx_t [num_trees + 1];
            for ( int i = 0; i < num_trees + 1; ++i )
            {
                ghost->tree_offsets[i] = input->tree_offsets[i];
            }

            ghost->proc_offsets = new p4est_locidx_t [mpisize + 1];
            for ( int i = 0; i < mpisize + 1; ++i )
            {
                ghost->proc_offsets[i] = input->proc_offsets[i];
            }

            ghost->mirror_tree_offsets = new p4est_locidx_t [num_trees + 1];
            for ( int i = 0; i < num_trees + 1; ++i )
            {
                ghost->mirror_tree_offsets[i] = input->mirror_tree_offsets[i];
            }

            ghost->mirror_proc_offsets = new p4est_locidx_t [mpisize + 1];
            for ( int i = 0; i < mpisize + 1; ++i )
            {
                ghost->mirror_proc_offsets[i] = input->mirror_proc_offsets[i];
            }

            ghost->mirror_proc_mirrors = new p4est_locidx_t [ ghost->mirror_proc_offsets[mpisize] ];
            for ( int i = 0; i < ghost->mirror_proc_offsets[mpisize]; ++i )
            {
                ghost->mirror_proc_mirrors[i] = input->mirror_proc_mirrors[i];
            }

            ghost->mirror_proc_front_offsets = new p4est_locidx_t [mpisize + 1];
            for ( int i = 0; i < mpisize + 1; ++i )
            {
                ghost->mirror_proc_front_offsets[i] = input->mirror_proc_front_offsets[i];
            }

            ghost->mirror_proc_fronts = new p4est_locidx_t [ ghost->mirror_proc_front_offsets[mpisize] ];
            for ( int i = 0; i < ghost->mirror_proc_front_offsets[mpisize]; ++i )
            {
                ghost->mirror_proc_fronts[i] = input->mirror_proc_fronts[i];
            }

            // quadrants
            size_t elem_size = input->ghosts.elem_size;
            sc_array_init ( &ghost->ghosts, elem_size );
            sc_array_copy ( &ghost->ghosts, &input->ghosts );
            for ( int i = 0; i < ghost->ghosts.elem_count; ++i )
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( p8est_quadrant_array_index ( &input->ghosts, i ) );
                pXest_set_quad_data ( p8est_quadrant_array_index ( &ghost->ghosts, i ), *q_ptr );
            }

            elem_size = input->mirrors.elem_size;
            sc_array_init ( &ghost->mirrors, elem_size );
            sc_array_copy ( &ghost->mirrors, &input->mirrors );
            for ( int i = 0; i < ghost->mirrors.elem_count; ++i )
            {
                p8est_quadrant_t* new_quad = p8est_quadrant_array_index ( &ghost->mirrors, i );
                p8est_quadrant_t* old_quad = p8est_quadrant_array_index ( &input->mirrors, i );

                new_quad->p.piggy3.which_tree = old_quad->p.piggy3.which_tree;
                new_quad->p.piggy3.local_num = old_quad->p.piggy3.local_num;
            }
        }

        void pXest_build_conn ( int num_vertices, int num_cells, const std::vector<double>& vertices, const std::vector<int>& tree_to_vertices,
                                p4est_connectivity_t*& conn )
        {
            LOG_DEBUG ( 1, "num_vertices: " << num_vertices << " num_cells: " << num_cells << " vertices.size(): " << vertices.size ( )
                        << " tree_to_vertices.size(): " << tree_to_vertices.size ( ) );

            assert ( vertices.size ( ) > 0 );
            assert ( tree_to_vertices.size ( ) > 0 );
            assert ( num_cells * 4 == tree_to_vertices.size ( ) );

            if ( DEBUG_LEVEL > 0 )
            {
                int min_id = 1e8;
                int max_id = -1;
                for ( int l = 0; l < tree_to_vertices.size ( ); ++l )
                {
                    if ( tree_to_vertices[l] < min_id )
                    {
                        min_id = tree_to_vertices[l];
                    }
                    if ( tree_to_vertices[l] > max_id )
                    {
                        max_id = tree_to_vertices[l];
                    }
                }
                LOG_DEBUG ( 1, "Max(tree_to_vertices) = " << max_id << " Min(tree_to_vertices) = " << min_id );
            }

            // create empty connectivity
            conn = p4est_connectivity_new ( num_vertices, num_cells, 0, 0 );

            // copy data
            std::copy ( vertices.begin ( ), vertices.end ( ), conn->vertices );
            std::copy ( tree_to_vertices.begin ( ), tree_to_vertices.end ( ), conn->tree_to_vertex );

            for ( int tree = 0; tree < conn->num_trees; ++tree )
            {
                for ( int face = 0; face < P4EST_FACES; ++face )
                {
                    conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
                    conn->tree_to_face[P4EST_FACES * tree + face] = face;

                    LOG_DEBUG ( 3, " tree " << tree << " , face " << face << " , index1 " << P4EST_FACES * tree + face << " index2 " << P4EST_FACES * tree + face );
                }
            }

            // check if everything is correct
            assert ( p4est_connectivity_is_valid ( conn ) );
            p4est_connectivity_complete ( conn );
        }

        void pXest_build_conn ( int num_vertices, int num_cells, const std::vector<double>& vertices, const std::vector<int>& tree_to_vertices,
                                p8est_connectivity_t*& conn )
        {
            LOG_DEBUG ( 1, "num_vertices: " << num_vertices << " num_cells: " << num_cells << " vertices.size(): " << vertices.size ( )
                        << " tree_to_vertices.size(): " << tree_to_vertices.size ( ) );

            assert ( vertices.size ( ) > 0 );
            assert ( tree_to_vertices.size ( ) > 0 );
            assert ( num_cells * 8 == tree_to_vertices.size ( ) );

            if ( DEBUG_LEVEL > 0 )
            {
                int min_id = 1e8;
                int max_id = -1;
                for ( int l = 0; l < tree_to_vertices.size ( ); ++l )
                {
                    if ( tree_to_vertices[l] < min_id )
                    {
                        min_id = tree_to_vertices[l];
                    }
                    if ( tree_to_vertices[l] > max_id )
                    {
                        max_id = tree_to_vertices[l];
                    }
                }
                LOG_DEBUG ( 1, "Max(tree_to_vertices) = " << max_id << " Min(tree_to_vertices) = " << min_id );
            }

            // create empty connectivity
            conn = p8est_connectivity_new ( num_vertices, num_cells, 0, 0, 0, 0 );

            // copy data
            std::copy ( vertices.begin ( ), vertices.end ( ), conn->vertices );
            std::copy ( tree_to_vertices.begin ( ), tree_to_vertices.end ( ), conn->tree_to_vertex );

            for ( int tree = 0; tree < conn->num_trees; ++tree )
            {
                for ( int face = 0; face < P8EST_FACES; ++face )
                {
                    conn->tree_to_tree[P8EST_FACES * tree + face] = tree;
                    conn->tree_to_face[P8EST_FACES * tree + face] = face;

                    LOG_DEBUG ( 3, " tree " << tree << " , face " << face << " , index1 " << P8EST_FACES * tree + face << " index2 " << P8EST_FACES * tree + face );
                }
            }

            // check if everything is correct
            assert ( p8est_connectivity_is_valid ( conn ) );
            p8est_connectivity_complete ( conn );
        }

        static int reorder_comp ( const void *a, const void *b )
        {
            const int *A = ( const int * ) a;
            const int *B = ( const int * ) b;

            if ( A[0] < B[0] )
            {
                return -1;
            }
            else if ( B[0] < A[0] )
            {
                return 1;
            }
            else
            {
                return (A[1] - B[1] );
            }
        }

        void pXest_reorder_conn ( const std::vector<int>& part, p4est_connectivity_t * conn, std::vector<size_t>& perm )
        {
            int n = part.size ( );
            int i;
            size_t zz;
            sc_array_t *newid;
            size_t *zp;
            sc_array_t *sorter;
            int *ip;
            perm.resize ( n );

            /* now that everyone has part, each process computes the renumbering
             * for itself*/
            newid = sc_array_new_size ( sizeof (size_t ), ( size_t ) n );
            sorter = sc_array_new_size ( 2 * sizeof (int ), ( size_t ) n );
            for ( i = 0; i < n; i++ )
            {
                ip = ( int * ) sc_array_index ( sorter, i );
                ip[0] = part[i];
                ip[1] = i;
            }

            /* sort current index by partition given */
            /* this will be the same on every process because the comparison operation
             * does not allow equality between different trees */
            sc_array_sort ( sorter, reorder_comp );
            for ( i = 0; i < n; i++ )
            {
                ip = ( int * ) sc_array_index ( sorter, i );
                zp = ( size_t * ) sc_array_index ( newid, ip[1] );
                *zp = i;

                perm[ip[1]] = i;
            }
            sc_array_destroy ( sorter );

            p4est_connectivity_permute ( conn, newid, 1 );

            sc_array_destroy ( newid );
        }

        void pXest_reorder_conn ( const std::vector<int>& part, p8est_connectivity_t * conn, std::vector<size_t>& perm )
        {
            int n = part.size ( );
            int i;
            size_t zz;
            sc_array_t *newid;
            size_t *zp;
            sc_array_t *sorter;
            int *ip;

            perm.resize ( n );

            /* now that everyone has part, each process computes the renumbering
             * for itself*/
            newid = sc_array_new_size ( sizeof (size_t ), ( size_t ) n );
            sorter = sc_array_new_size ( 2 * sizeof (int ), ( size_t ) n );
            for ( i = 0; i < n; i++ )
            {
                ip = ( int * ) sc_array_index ( sorter, i );
                ip[0] = part[i];
                ip[1] = i;
            }

            /* sort current index by partition given */
            /* this will be the same on every process because the comparison operation
             * does not allow equality between different trees */
            sc_array_sort ( sorter, reorder_comp );
            for ( i = 0; i < n; i++ )
            {
                ip = ( int * ) sc_array_index ( sorter, i );
                zp = ( size_t * ) sc_array_index ( newid, ip[1] );
                *zp = i;

                perm[ip[1]] = i;
            }
            sc_array_destroy ( sorter );

            p8est_connectivity_permute ( conn, newid, 1 );

            sc_array_destroy ( newid );
        }

        p4est_gloidx_t pXest_partition_forest ( p4est_t * forest, int partition_for_coarsening, std::vector<p4est_locidx_t>& num_quads_per_proc )
        {
            p4est_gloidx_t global_shipped = 0;
            const p4est_gloidx_t global_num_quadrants = forest->global_num_quadrants;

            P4EST_ASSERT ( p4est_is_valid ( forest ) );
            P4EST_GLOBAL_PRODUCTIONF ( "Into " P4EST_STRING "_partition with %lld total quadrants\n", ( long long ) forest->global_num_quadrants );

            /* this function does nothing in a serial setup */
            if ( forest->mpisize == 1 )
            {
                P4EST_GLOBAL_PRODUCTION ( "Done " P4EST_STRING "_partition no shipping\n" );
                return global_shipped;
            }

            p4est_log_indent_push ( );

            // FEATURE_TODO: remove comments, problem: p4est_partition_for_coarsening is not contained in p4est.h, although it is defined
            // in p4est.c. Either, one has modify p4est.h accordingly, or, one has to include p4est_partition_for_coarsening in hiflow
            /* correct partition */
            /*
                        if (partition_for_coarsening) 
                        {
                            p4est_gloidx_t num_corrected = p4est_partition_for_coarsening (forest, &num_quads_per_proc[0]);
                            P4EST_GLOBAL_INFOF ("Designated partition for coarsening %lld quadrants moved\n", (long long) num_corrected);
                        }
             */
            /* run the partition algorithm with proper quadrant counts */
            global_shipped = p4est_partition_given ( forest, &num_quads_per_proc[0] );

            /* check validity of the p4est */
            P4EST_ASSERT ( p4est_is_valid ( forest ) );

            // update owner ranks in QuadData
            int my_rank = forest->mpirank;

            treeId first_treeId = pXest_get_first_local_treeId ( forest );
            treeId last_treeId = pXest_get_last_local_treeId ( forest );

            for ( treeId jt = first_treeId; jt <= last_treeId; ++jt )
            {
                p4est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                int num_quads = tree->quadrants.elem_count;

                // Loop over all local quads in tree
                for ( int jq = 0; jq < num_quads; ++jq )
                {
                    p4est_quadrant_t* quad = p4est_quadrant_array_index ( &tree->quadrants, jq );

                    QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
                    q_ptr->set_owner_rank ( my_rank );
                }
            }

            p4est_log_indent_pop ( );
            P4EST_GLOBAL_PRODUCTIONF ( "Done " P4EST_STRING "_partition shipped %lld quadrants %.3g%%\n",
                                       ( long long ) global_shipped, global_shipped * 100. / global_num_quadrants );

            return global_shipped;
        }

        p4est_gloidx_t pXest_partition_forest ( p8est_t * forest, int partition_for_coarsening, std::vector<p4est_locidx_t>& num_quads_per_proc )
        {
            p4est_gloidx_t global_shipped = 0;
            const p4est_gloidx_t global_num_quadrants = forest->global_num_quadrants;

            P4EST_ASSERT ( p8est_is_valid ( forest ) );
            P4EST_GLOBAL_PRODUCTIONF ( "Into " P4EST_STRING "_partition with %lld total quadrants\n", ( long long ) forest->global_num_quadrants );

            /* this function does nothing in a serial setup */
            if ( forest->mpisize == 1 )
            {
                P4EST_GLOBAL_PRODUCTION ( "Done " P4EST_STRING "_partition no shipping\n" );
                return global_shipped;
            }

            p4est_log_indent_push ( );

            // FEATURE_TODO: remove comments, problem: p8est_partition_for_coarsening is not contained in p8est.h, although it is defined
            // in p8est.c. Either, one has modify p8est.h accordingly, or, one has to include p8est_partition_for_coarsening in hiflow
            /* correct partition */
            /*
                        if (partition_for_coarsening) 
                        {
                            p4est_gloidx_t num_corrected = p8est_partition_for_coarsening (forest, &num_quads_per_proc[0]);
                            P4EST_GLOBAL_INFOF ("Designated partition for coarsening %lld quadrants moved\n", (long long) num_corrected);
                        }
             */
            /* run the partition algorithm with proper quadrant counts */
            global_shipped = p8est_partition_given ( forest, &num_quads_per_proc[0] );

            /* check validity of the p4est */
            P4EST_ASSERT ( p8est_is_valid ( forest ) );

            // update owner ranks in QuadData
            int my_rank = forest->mpirank;

            treeId first_treeId = pXest_get_first_local_treeId ( forest );
            treeId last_treeId = pXest_get_last_local_treeId ( forest );

            for ( treeId jt = first_treeId; jt <= last_treeId; ++jt )
            {
                p8est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                int num_quads = tree->quadrants.elem_count;

                // Loop over all local quads in tree
                for ( int jq = 0; jq < num_quads; ++jq )
                {
                    p8est_quadrant_t* quad = p8est_quadrant_array_index ( &tree->quadrants, jq );

                    QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
                    assert ( q_ptr != NULL );
                    q_ptr->set_owner_rank ( my_rank );
                }
            }

            p4est_log_indent_pop ( );
            P4EST_GLOBAL_PRODUCTIONF ( "Done " P4EST_STRING "_partition shipped %lld quadrants %.3g%%\n",
                                       ( long long ) global_shipped, global_shipped * 100. / global_num_quadrants );

            return global_shipped;
        }

        std::vector<int> pXest_get_partition_from_initial_forest ( p4est_t* forest )
        {
            int first_local_tree = forest->first_local_tree;
            int last_local_tree = forest->last_local_tree;
            int local_num_quad = last_local_tree - first_local_tree + 1;

            std::vector<int> first_remote_tree ( forest->mpisize );
            std::vector<int> last_remote_tree ( forest->mpisize );

            MPI_Allgather ( &first_local_tree, 1, MPI_INT, vec2ptr ( first_remote_tree ), 1, MPI_INT, forest->mpicomm );
            MPI_Allgather ( &last_local_tree, 1, MPI_INT, vec2ptr ( last_remote_tree ), 1, MPI_INT, forest->mpicomm );

            std::vector<int> partition ( forest->trees->elem_count );

            for ( int p = 0; p < forest->mpisize; ++p )
            {
                for ( int l = first_remote_tree[p]; l <= last_remote_tree[p]; ++l )
                {
                    partition[l] = p;
                }
            }
            return partition;
        }

        std::vector<int> pXest_get_partition_from_initial_forest ( p8est_t* forest )
        {
            int first_local_tree = forest->first_local_tree;
            int last_local_tree = forest->last_local_tree;
            int local_num_quad = last_local_tree - first_local_tree + 1;

            std::vector<int> first_remote_tree ( forest->mpisize );
            std::vector<int> last_remote_tree ( forest->mpisize );

            MPI_Allgather ( &first_local_tree, 1, MPI_INT, vec2ptr ( first_remote_tree ), 1, MPI_INT, forest->mpicomm );
            MPI_Allgather ( &last_local_tree, 1, MPI_INT, vec2ptr ( last_remote_tree ), 1, MPI_INT, forest->mpicomm );

            std::vector<int> partition ( forest->trees->elem_count );

            for ( int p = 0; p < forest->mpisize; ++p )
            {
                for ( int l = first_remote_tree[p]; l <= last_remote_tree[p]; ++l )
                {
                    partition[l] = p;
                }
            }
            return partition;
        }

        // PERF_TODO check efficiency, maybe use MPI_PACK, check whether deadlock can occure -> use nonblocking communication modes

        void pXest_distribute_coarse_maps ( MPI_Comm comm, int master_rank,
                                            std::vector<int>& partitioning, std::vector<int>& coarse_cell_ids, std::vector<int>& local_cell_ids,
                                            EntityToQuadMap* coarse_map, EntityToQuadMap& map )
        {
            int rank;
            int num_ranks;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &num_ranks );

            int size_map = 10;
            int local_ids = 0;
            int num_cells = partitioning.size ( );
            for ( std::vector<int>::iterator it = partitioning.begin ( ); it != partitioning.end ( ); ++it )
            {
                if ( *it == rank )
                {
                    ++local_ids;
                }
            }
            map.clear ( );

            if ( rank == master_rank )
            {
                LOG_DEBUG ( 3, " i am master " );
                // Loop over all cell indices in coarse mesh
                for ( EntityNumber jc = 0; jc < num_cells; ++jc )
                {
                    int proc = partitioning[jc];
                    Id coarse_id = coarse_cell_ids[jc];
                    QuadCoord coord = ( *coarse_map )[coarse_id];

                    LOG_DEBUG ( 3, " considered id " << coarse_id );
                    if ( proc != master_rank )
                    {
                        LOG_DEBUG ( 3, " send coarse id " << coarse_id << " to proc " << proc << " with base tag " << jc * size_map );

                        MPI_Send ( &coarse_id, 1, MPI_INT32_T, proc, jc * size_map + 0, comm );
                        MPI_Send ( &coord.tree, 1, MPI_INT32_T, proc, jc * size_map + 1, comm );
                        MPI_Send ( &coord.level, 1, MPI_INT, proc, jc * size_map + 2, comm );
                        MPI_Send ( &coord.x, 1, MPI_INT32_T, proc, jc * size_map + 3, comm );
                        MPI_Send ( &coord.y, 1, MPI_INT32_T, proc, jc * size_map + 4, comm );
                        MPI_Send ( &coord.z, 1, MPI_INT32_T, proc, jc * size_map + 5, comm );
                        MPI_Send ( &coord.localId, 1, MPI_INT, proc, jc * size_map + 6, comm );
                        MPI_Send ( &coord.tdim, 1, MPI_INT, proc, jc * size_map + 7, comm );
                        MPI_Send ( &coord.coarseId, 1, MPI_INT, proc, jc * size_map + 8, comm );
                    }
                    else
                    {
                        LOG_DEBUG ( 3, " keep id " << coarse_id << " as local id " << local_cell_ids[jc] );
                        assert ( local_cell_ids[jc] >= 0 );
                        map.insert ( std::pair<Id, QuadCoord>( local_cell_ids[jc], coord ) );
                    }
                }
            }
            else
            {
                for ( EntityNumber jc = 0; jc < num_cells; ++jc )
                {
                    int proc = partitioning[jc];
                    if ( proc == rank )
                    {
                        int tag_id = jc;
                        Id coarse_id;
                        QuadCoord coord ( 0, 0, 0, 0, 0, 0, 0, 0 );
                        LOG_DEBUG ( 3, " proc " << rank << " receives data from " << master_rank << " with base tag " << tag_id * size_map );

                        std::vector<MPI_Status> status ( size_map );
                        MPI_Recv ( &coarse_id, 1, MPI_INT32_T, master_rank, tag_id * size_map + 0, comm, &status[0] );
                        MPI_Recv ( &coord.tree, 1, MPI_INT32_T, master_rank, tag_id * size_map + 1, comm, &status[1] );
                        MPI_Recv ( &coord.level, 1, MPI_INT, master_rank, tag_id * size_map + 2, comm, &status[2] );
                        MPI_Recv ( &coord.x, 1, MPI_INT32_T, master_rank, tag_id * size_map + 3, comm, &status[3] );
                        MPI_Recv ( &coord.y, 1, MPI_INT32_T, master_rank, tag_id * size_map + 4, comm, &status[4] );
                        MPI_Recv ( &coord.z, 1, MPI_INT32_T, master_rank, tag_id * size_map + 5, comm, &status[5] );
                        MPI_Recv ( &coord.localId, 1, MPI_INT, master_rank, tag_id * size_map + 6, comm, &status[6] );
                        MPI_Recv ( &coord.tdim, 1, MPI_INT, master_rank, tag_id * size_map + 7, comm, &status[7] );
                        MPI_Recv ( &coord.coarseId, 1, MPI_INT, master_rank, tag_id * size_map + 8, comm, &status[8] );

                        LOG_DEBUG ( 3, " proc " << rank << " received coarse id " << coarse_id << " which is pushed back as local id "
                                    << local_cell_ids[jc] );

                        map.insert ( std::pair<Id, QuadCoord>( local_cell_ids[jc], coord ) );
                        assert ( local_cell_ids[jc] >= 0 );
                    }
                }
            }
        }

        void pXest_extract_local_maps ( MPI_Comm comm, int master_rank,
                                        std::vector<int>& partitioning, std::vector<int>& coarse_cell_ids, std::vector<int>& local_cell_ids,
                                        EntityToQuadMap* coarse_map, EntityToQuadMap& map )
        {
            int rank;
            int num_ranks;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &num_ranks );

            int size_map = 10;
            int local_ids = 0;
            int num_cells = partitioning.size ( );
            map.clear ( );

            // Loop over all cell indices in coarse mesh
            for ( EntityNumber jc = 0; jc < num_cells; ++jc )
            {
                int proc = partitioning[jc];
                Id coarse_id = coarse_cell_ids[jc];
                QuadCoord coord = ( *coarse_map )[coarse_id];

                LOG_DEBUG ( 3, " considered id " << coarse_id );
                if ( proc == rank )
                {
                    LOG_DEBUG ( 3, " extract id " << coarse_id << " as local id " << local_cell_ids[jc] );
                    assert ( local_cell_ids[jc] >= 0 );
                    map.insert ( std::pair<Id, QuadCoord>( local_cell_ids[jc], coord ) );
                }
            }
        }

        mortonId pXest_get_mortonId ( p4est_quadrant_t* quad, int level )
        {
            if ( quad->level >= level )
            {
                return p4est_quadrant_linear_id ( quad, level );
            }
            else
            {
                return -1;
            }
        }

        mortonId pXest_get_mortonId ( p8est_quadrant_t* quad, int level )
        {
            if ( quad->level >= level )
            {
                return p8est_quadrant_linear_id ( quad, level );
            }
            else
            {
                return -1;
            }
        }

        mortonId pXest_get_last_mortonId_in_tree ( p4est_tree_t* tree, int level )
        {
            int num_quads = tree->quadrants.elem_count;
            if ( num_quads == 0 )
            {
                return -1;
            }
            if ( num_quads == 1 )
            {
                return p4est_quadrant_linear_id ( &tree->first_desc, level );
            }
            else
            {
                return p4est_quadrant_linear_id ( &tree->last_desc, level );
            }
        }

        mortonId pXest_get_last_mortonId_in_tree ( p8est_tree_t* tree, int level )
        {
            int num_quads = tree->quadrants.elem_count;
            if ( num_quads == 0 )
            {
                return -1;
            }
            if ( num_quads == 1 )
            {
                return p8est_quadrant_linear_id ( &tree->first_desc, level );
            }
            else
            {
                return p8est_quadrant_linear_id ( &tree->last_desc, level );
            }
        }

        bool pXest_quad_arrays_in_forest_sorted ( p4est_t* forest )
        {
            for ( int jt = 0; jt < forest->trees->elem_count; ++jt )
            {
                p4est_tree_t* cur_tree = pXest_get_tree_in_forest ( forest, jt );
                if ( !sc_array_is_sorted ( &cur_tree->quadrants, pXest_compare_quad4_in_tree ) )
                {
                    return false;
                }
            }
            return true;
        }

        bool pXest_quad_arrays_in_forest_sorted ( p8est_t* forest )
        {
            for ( int jt = 0; jt < forest->trees->elem_count; ++jt )
            {
                p8est_tree_t* cur_tree = pXest_get_tree_in_forest ( forest, jt );
                if ( !sc_array_is_sorted ( &cur_tree->quadrants, pXest_compare_quad8_in_tree ) )
                {
                    return false;
                }
            }
            return true;
        }

        void pXest_sort_quad_arrays_in_forest ( p4est_t* forest )
        {
            for ( int jt = 0; jt < forest->trees->elem_count; ++jt )
            {
                p4est_tree_t* cur_tree = pXest_get_tree_in_forest ( forest, jt );
                if ( !sc_array_is_sorted ( &cur_tree->quadrants, pXest_compare_quad4_in_tree ) )
                {
                    sc_array_sort ( &cur_tree->quadrants, pXest_compare_quad4_in_tree );
                }
            }
        }

        void pXest_sort_quad_arrays_in_forest ( p8est_t* forest )
        {
            for ( int jt = 0; jt < forest->trees->elem_count; ++jt )
            {
                p8est_tree_t* cur_tree = pXest_get_tree_in_forest ( forest, jt );
                if ( !sc_array_is_sorted ( &cur_tree->quadrants, pXest_compare_quad8_in_tree ) )
                {
                    sc_array_sort ( &cur_tree->quadrants, pXest_compare_quad8_in_tree );
                }
            }
        }

        int64_t pXest_get_quad_pos_in_tree_array ( p4est_tree_t* tree, p4est_quadrant_t* quad )
        {
            // binary search
            int64_t pos = sc_array_bsearch ( &tree->quadrants, quad, pXest_compare_quad4_in_tree );
            return pos;

            // Note: the following code can be used for testing sc_array_bsearch. In this case, simply remove the above return pos.
            // In principle, the result of this function should not be changed
            if ( pos >= 0 )
            {
                LOG_DEBUG ( 2, "Binary search in tree quadrants SUCCESS" );
                return pos;
            }

            // naive search
            const int num_quads = tree->quadrants.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, jq );
                if ( p4est_quadrant_is_equal ( cur_quad, quad ) )
                {
                    LOG_DEBUG ( 2, "Binary search in tree quadrants FAILED" );
                    assert ( false );
                    return jq;
                }
            }
            LOG_DEBUG ( 2, "Search in tree quadrants FAILED" );
            return -1;
        }

        int64_t pXest_get_quad_pos_in_tree_array ( p8est_tree_t* tree, p8est_quadrant_t* quad )
        {
            // binary search
            int64_t pos = sc_array_bsearch ( &tree->quadrants, quad, pXest_compare_quad8_in_tree );
            return pos;

            // Note: the following code can be used for testing sc_array_bsearch. In this case, simply remove the above return pos.
            // In principle, the result of this function should not be changed
            if ( pos >= 0 )
            {
                LOG_DEBUG ( 2, "Binary search in tree quadrants SUCCESS" );
                return pos;
            }

            // naive search
            const int num_quads = tree->quadrants.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, jq );
                if ( p8est_quadrant_is_equal ( cur_quad, quad ) )
                {
                    LOG_DEBUG ( 2, "Binary search in tree quadrants FAILED" );
                    assert ( false );
                    return jq;
                }
            }
            LOG_DEBUG ( 2, "Search in tree quadrants FAILED" );
            return -1;
        }

        int64_t pXest_get_quad_pos_in_ghost_array ( p4est_ghost_t* ghost, p4est_quadrant_t* quad, treeId tree_id )
        {
            // binary search for quadrant in ghost layer
            int64_t pos = p4est_ghost_bsearch ( ghost, -1, tree_id, quad );
            return pos;

            // Note: the following code can be used for testing sc_array_bsearch. In this case, simply remove the above return pos.
            // In principle, the result of this function should not be changed
            if ( pos >= 0 )
            {
                LOG_DEBUG ( 2, " Binary ghost layer search SUCCESS " );
                return pos;
            }
            LOG_DEBUG ( 2, " Binary ghost layer search failed " );

            // naive search
            const int num_quads = ghost->ghosts.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_ghost ( ghost, jq );
                LOG_DEBUG ( 3, " QuadData " << pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) );

                if ( pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) == tree_id && p4est_quadrant_is_equal ( cur_quad, quad ) )
                {
                    return jq;
                }
            }
            return -1;
        }

        int64_t pXest_get_quad_pos_in_ghost_array ( p8est_ghost_t* ghost, p8est_quadrant_t* quad, treeId tree_id )
        {
            // binary search for quadrant in ghost layer
            int64_t pos = p8est_ghost_bsearch ( ghost, -1, tree_id, quad );
            return pos;

            // Note: the following code can be used for testing sc_array_bsearch. In this case, simply remove the above return pos.
            // In principle, the result of this function should not be changed
            if ( pos >= 0 )
            {
                LOG_DEBUG ( 2, " Binary ghost layer search SUCCESS " );
                return pos;
            }
            LOG_DEBUG ( 2, " Binary ghost layer search failed " );

            // naive search
            const int num_quads = ghost->ghosts.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_ghost ( ghost, jq );
                LOG_DEBUG ( 3, " QuadData " << pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) );

                if ( pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) == tree_id && p8est_quadrant_is_equal ( cur_quad, quad ) )
                {
                    return jq;
                }
            }
            return -1;
        }

        // PERF_TODO make more efficient by using range of possible mortonIds (computed from input quad) 
        // and binary search
        // PERF_TODO: BOTTLENECK

        std::vector<int64_t> pXest_get_leaf_desc_in_tree_array ( p4est_tree_t* tree, p4est_quadrant_t* quad )
        {
            std::vector<int64_t> leaf_pos;
            const int num_quads = tree->quadrants.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, jq );
                if ( p4est_quadrant_is_ancestor ( quad, cur_quad ) )
                {
                    leaf_pos.push_back ( jq );
                }
            }
            return leaf_pos;
        }

        // PERF_TODO: BOTTLENECK

        std::vector<int64_t> pXest_get_leaf_desc_in_tree_array ( p8est_tree_t* tree, p8est_quadrant_t* quad )
        {
            std::vector<int64_t> leaf_pos;
            const int num_quads = tree->quadrants.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, jq );
                if ( p8est_quadrant_is_ancestor ( quad, cur_quad ) )
                {
                    leaf_pos.push_back ( jq );
                }
            }
            return leaf_pos;
        }

        // PERF_TODO make more efficient by using range of possible mortonIds (computed from input quad) 
        // and binary search
        // PERF_TODO: BOTTLENECK

        std::vector<int64_t> pXest_get_leaf_desc_in_ghost_array ( p4est_ghost_t* ghost, p4est_quadrant_t* quad, treeId tree_id )
        {
            std::vector<int64_t> leaf_pos;
            const int num_quads = ghost->ghosts.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p4est_quadrant_t* cur_quad = p4est_quadrant_array_index ( &ghost->ghosts, jq );
                LOG_DEBUG ( 3, " QuadData " << pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) );

                if ( pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) == tree_id && p4est_quadrant_is_ancestor ( quad, cur_quad ) )
                {
                    leaf_pos.push_back ( jq );
                }
            }
            return leaf_pos;
        }

        // PERF_TODO: BOTTLENECK

        std::vector<int64_t> pXest_get_leaf_desc_in_ghost_array ( p8est_ghost_t* ghost, p8est_quadrant_t* quad, treeId tree_id )
        {
            std::vector<int64_t> leaf_pos;
            const int num_quads = ghost->ghosts.elem_count;
            for ( int jq = 0; jq < num_quads; ++jq )
            {
                p8est_quadrant_t* cur_quad = p8est_quadrant_array_index ( &ghost->ghosts, jq );
                LOG_DEBUG ( 3, " QuadData " << pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) );

                if ( pXest_get_quad_data_ptr ( cur_quad )->get_tree_id ( ) == tree_id && p8est_quadrant_is_ancestor ( quad, cur_quad ) )
                {
                    leaf_pos.push_back ( jq );
                }
            }
            return leaf_pos;
        }

        QuadData* pXest_get_quad_data_ptr ( p4est_tree_t* tree, p4est_quadrant_t* quad )
        {
            const int64_t pos = pXest_get_quad_pos_in_tree_array ( tree, quad );
            if ( pos >= 0 )
            {
                p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, pos );
                return pXest_get_quad_data_ptr ( cur_quad );
            }
            else
            {
                return NULL;
            }
        }

        QuadData* pXest_get_quad_data_ptr ( p8est_tree_t* tree, p8est_quadrant_t* quad )
        {
            const int64_t pos = pXest_get_quad_pos_in_tree_array ( tree, quad );
            if ( pos >= 0 )
            {
                p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, pos );
                return pXest_get_quad_data_ptr ( cur_quad );
            }
            else
            {
                return NULL;
            }
        }

        QuadData* pXest_get_quad_data_ptr ( p4est_ghost_t* ghost, p4est_quadrant_t* quad, treeId tree_id )
        {
            const int64_t pos = pXest_get_quad_pos_in_ghost_array ( ghost, quad, tree_id );
            if ( pos >= 0 )
            {
                p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_ghost ( ghost, pos );
                return pXest_get_quad_data_ptr ( cur_quad );
            }
            else
            {
                return NULL;
            }
        }

        QuadData* pXest_get_quad_data_ptr ( p8est_ghost_t* ghost, p8est_quadrant_t* quad, treeId tree_id )
        {
            const int64_t pos = pXest_get_quad_pos_in_ghost_array ( ghost, quad, tree_id );
            if ( pos >= 0 )
            {
                p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_ghost ( ghost, pos );
                return pXest_get_quad_data_ptr ( cur_quad );
            }
            else
            {
                return NULL;
            }
        }

        void pXest_set_quad_data ( p4est_quadrant_t* quad, QuadData& data )
        {
            QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
            q_ptr->set ( data );
        }

        void pXest_set_quad_data ( p8est_quadrant_t* quad, QuadData& data )
        {
            QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
            q_ptr->set ( data );
        }

        void pXest_get_quad_data ( p4est_quadrant_t* quad, QuadData& data )
        {
            QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
            data.set ( *q_ptr );
        }

        void pXest_get_quad_data ( p8est_quadrant_t* quad, QuadData& data )
        {
            QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
            data.set ( *q_ptr );
        }

        QuadCoord pXest_get_quad_coord ( p4est_quadrant_t* quad, treeId tree_id )
        {
            QuadCoord coord ( tree_id, quad->level, quad->x, quad->y, 0, -1, 2, -1 );
            return coord;
        }

        QuadCoord pXest_get_quad_coord ( p8est_quadrant_t* quad, treeId tree_id )
        {
            QuadCoord coord ( tree_id, quad->level, quad->x, quad->y, quad->z, -1, 3, -1 );
            return coord;
        }

        void pXest_build_ref_tree ( std::vector<p4est_quadrant_t>& quads, RefinementTree* tree )
        {
            int num_quads = quads.size ( );
            p4est_quadrant_t ancestor = quads[0];

            LOG_DEBUG ( 2, " Build RefTree with " << num_quads << " quadrants " );
            LOG_DEBUG ( 2, " Root cell is on level " << ( int ) ancestor.level );

            // sort descendants according to their level and their family
            std::vector< std::map< int, std::vector<int> > > quads_index;

            // BUG_TODO check resize
            quads_index.reserve ( 5 * num_quads );

            for ( int jq = 0; jq < num_quads; ++jq )
            {
                // level index
                int diff_level = ( int ) quads[jq].level - ( int ) ancestor.level;
                if ( diff_level >= quads_index.size ( ) )
                {
                    quads_index.resize ( diff_level + 1 );
                }

                // family index on current level 
                int family = 0;
                for ( int jl = 1; jl < diff_level; ++jl )
                {
                    // maximum number of descendant families on level quad[jq].level a quad on level ancestor.level+jl can have
                    int max_num_desc_family = static_cast < int > ( std::pow ( static_cast < double > ( 4 ), diff_level - jl - 1 ) );
                    int ancestor_number = p4est_quadrant_ancestor_id ( &quads[jq], ancestor.level + jl );
                    LOG_DEBUG ( 2, "Intermediate level " << jl << " max_num_desc_families " << max_num_desc_family << " ancestor_number " << ancestor_number );

                    family += ancestor_number * max_num_desc_family;
                }

                LOG_DEBUG ( 2, "Current quad " << jq << " is on level " << ( int ) quads[jq].level << " and has family index " << family
                            << " and x " << quads[jq].x << " and y " << quads[jq].y );

                // put quadrant number into map
                std::map<int, std::vector<int> >::iterator it = quads_index[diff_level].find ( family );
                if ( it == quads_index[diff_level].end ( ) )
                {
                    // First occurence of family
                    if ( jq > 0 )
                    {
                        std::vector<int> tmp ( 4, -1 );
                        tmp[ pXest_z_order_to_ccw_order ( pXest_get_sibling_nr ( &quads[jq] ) )] = jq;
                        quads_index[diff_level].insert ( std::pair<int, std::vector<int> > ( family, tmp ) );
                    }
                    else
                    {
                        std::vector<int> tmp ( 1, jq );
                        quads_index[diff_level].insert ( std::pair<int, std::vector<int> > ( family, tmp ) );
                    }
                }
                else
                {
                    // Repeated occurence of family
                    it->second.at ( pXest_z_order_to_ccw_order ( pXest_get_sibling_nr ( &quads[jq] ) ) ) = jq;
                }
            }
            int num_levels = quads_index.size ( );

            if ( DEBUG_LEVEL >= 2 )
            {
                for ( int jl = 0; jl < num_levels; ++jl )
                {
                    std::cout << " Level: " << jl << std::endl;
                    for ( std::map<int, std::vector<int> >::iterator jf = quads_index[jl].begin ( ); jf != quads_index[jl].end ( ); ++jf )
                    {
                        std::cout << "    Family: " << jf->first << std::endl;

                    }
                }
            }

            std::map<int, int> quadrant_to_cell_index;

            // Loop over levels
            for ( int jl = 1; jl < num_levels; ++jl )
            {
                LOG_DEBUG ( 2, "Level: " << jl );

                int quad_level = ancestor.level + jl;
                int max_num_family_on_level = static_cast < int > ( std::pow ( static_cast < double > ( 4 ), jl - 1 ) );

                // Loop over families in level
                for ( std::map<int, std::vector<int> >::iterator jf = quads_index[jl].begin ( ); jf != quads_index[jl].end ( ); ++jf )
                {
                    std::vector<int> quadrant_numbers = jf->second;
                    int tmp[] = { 1, 2, 3, 4 };
                    std::vector<int> sub_cell_numbers ( tmp, tmp + sizeof (tmp ) / sizeof (int ) );

                    LOG_DEBUG ( 2, "   Family: " << jf->first << " consists of quadrant numbers " << quadrant_numbers[0] << " "
                                << quadrant_numbers[1] << " " << quadrant_numbers[2] << " " << quadrant_numbers[3] );

                    assert ( quadrant_numbers.size ( ) > 0 );

                    int parent_family = jf->first / 4;
                    int parent_sibling = jf->first - parent_family * 4;
                    int parent_quadrant_number = quads_index[jl - 1][parent_family][ pXest_z_order_to_ccw_order ( parent_sibling )];

                    int parent_cell_index = 0;
                    if ( jl > 1 )
                    {
                        parent_cell_index = quadrant_to_cell_index[parent_quadrant_number];
                    }

                    LOG_DEBUG ( 2, "   Family: " << jf->first << " has parent family: " << parent_family << ", parent sibling: " << parent_sibling
                                << ", parent_quadrant_number: " << parent_quadrant_number << " parent_cell_index: " << parent_cell_index );

                    // add subtree 
                    int num_cells = tree->num_cells ( );
                    tree->add_subtree ( parent_cell_index, sub_cell_numbers, quadrant_numbers );

                    for ( int jc = 0; jc < quadrant_numbers.size ( ); ++jc )
                    {
                        int cell_index = num_cells + jc;
                        quadrant_to_cell_index.insert ( std::pair<int, int> ( quadrant_numbers[jc], cell_index ) );
                    }
                }
            }
        }

        void pXest_build_ref_tree ( std::vector<p8est_quadrant_t>& quads, RefinementTree* tree )
        {
            int num_quads = quads.size ( );
            p8est_quadrant_t ancestor = quads[0];

            // sort descendants according to their level and their family
            std::vector< std::map< int, std::vector<int> > > quads_index;
            quads_index.reserve ( 5 * num_quads ); // BUG_TODO check resize

            for ( int jq = 0; jq < num_quads; ++jq )
            {
                // level index
                int diff_level = ( int ) quads[jq].level - ( int ) ancestor.level;
                if ( diff_level >= quads_index.size ( ) )
                {
                    quads_index.resize ( diff_level + 1 );
                }

                // family index on current level 
                int family = 0;
                for ( int jl = 1; jl < diff_level; ++jl )
                {
                    // maximum number of descendant families on level quad[jq].level a quad on level ancestor.level+jl can have
                    int max_num_desc_family = static_cast < int > ( std::pow ( static_cast < double > ( 8 ), diff_level - jl - 1 ) );
                    int ancestor_number = p8est_quadrant_ancestor_id ( &quads[jq], ancestor.level + jl );
                    family += ancestor_number * max_num_desc_family;
                }

                LOG_DEBUG ( 2, "Current quad " << jq << " is on level " << ( int ) quads[jq].level << " and has family index " << family );

                std::map<int, std::vector<int> >::iterator it = quads_index[diff_level].find ( family );
                if ( it == quads_index[diff_level].end ( ) )
                {// First occurence of family
                    if ( jq > 0 )
                    {
                        std::vector<int> tmp ( 8, -1 );
                        tmp[ pXest_z_order_to_ccw_order ( pXest_get_sibling_nr ( &quads[jq] ) )] = jq;
                        quads_index[diff_level].insert ( std::pair<int, std::vector<int> > ( family, tmp ) );
                    }
                    else
                    {
                        std::vector<int> tmp ( 1, jq );
                        quads_index[diff_level].insert ( std::pair<int, std::vector<int> > ( family, tmp ) );
                    }
                }
                else
                {// Repeated occurence of family
                    it->second.at ( pXest_z_order_to_ccw_order ( pXest_get_sibling_nr ( &quads[jq] ) ) ) = jq;
                }

            }
            int num_levels = quads_index.size ( );

            if ( DEBUG_LEVEL >= 2 )
            {
                for ( int jl = 0; jl < num_levels; ++jl )
                {
                    std::cout << " Level: " << jl << std::endl;
                    for ( std::map<int, std::vector<int> >::iterator jf = quads_index[jl].begin ( ); jf != quads_index[jl].end ( ); ++jf )
                    {
                        std::cout << "    Family: " << jf->first << std::endl;

                    }
                }
            }
            std::map<int, int> quadrant_to_cell_index;

            // Loop over levels
            for ( int jl = 1; jl < num_levels; ++jl )
            {
                LOG_DEBUG ( 2, "Level: " << jl );

                int quad_level = ancestor.level + jl;
                int max_num_family_on_level = static_cast < int > ( std::pow ( static_cast < double > ( 8 ), jl - 1 ) );

                // Loop over families in level
                for ( std::map<int, std::vector<int> >::iterator jf = quads_index[jl].begin ( ); jf != quads_index[jl].end ( ); ++jf )
                {
                    std::vector<int> quadrant_numbers = jf->second;
                    int tmp[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
                    std::vector<int> sub_cell_numbers ( tmp, tmp + sizeof (tmp ) / sizeof (int ) );

                    LOG_DEBUG ( 2, "   Family: " << jf->first
                                << " consists of quadrant numbers " << string_from_range ( quadrant_numbers.begin ( ), quadrant_numbers.end ( ) ) );

                    assert ( quadrant_numbers.size ( ) > 0 );

                    int parent_family = jf->first / 8;
                    int parent_sibling = jf->first - parent_family * 8;
                    int parent_quadrant_number = quads_index[jl - 1][parent_family][pXest_z_order_to_ccw_order ( parent_sibling )];

                    int parent_cell_index = 0;
                    if ( jl > 1 )
                    {
                        parent_cell_index = quadrant_to_cell_index[parent_quadrant_number];
                    }

                    LOG_DEBUG ( 2, "   Family: " << jf->first << " has parent family: " << parent_family << ", parent sibling: " << parent_sibling
                                << ", parent_quadrant_number: " << parent_quadrant_number << " parent_cell_index: " << parent_cell_index );

                    // add subtree 
                    int num_cells = tree->num_cells ( );
                    tree->add_subtree ( parent_cell_index, sub_cell_numbers, quadrant_numbers );

                    for ( int jc = 0; jc < quadrant_numbers.size ( ); ++jc )
                    {
                        int cell_index = num_cells + jc;
                        quadrant_to_cell_index.insert ( std::pair<int, int> ( quadrant_numbers[jc], cell_index ) );
                    }
                }
            }
        }

        // PERF_TODO make more efficient, avoid for loop

        int pXest_get_db_desc_nr ( p4est_quadrant_t* desc, p4est_quadrant_t* ancestor )
        {
            int diff_level = desc->level - ancestor->level;
            int desc_nr = 0;

            for ( int r = ancestor->level + 1; r <= desc->level; ++r )
            {
                int co_num_cells_on_level = static_cast < int > ( std::pow ( static_cast < double > ( 4 ), diff_level - ( r - ancestor->level ) ) );
                int sibling_nr_on_level = pXest_z_order_to_ccw_order ( p4est_quadrant_ancestor_id ( desc, r ) );

                desc_nr += sibling_nr_on_level * co_num_cells_on_level;
            }
            return desc_nr;
        }

        int pXest_get_db_desc_nr ( p8est_quadrant_t* desc, p8est_quadrant_t* ancestor )
        {
            int diff_level = desc->level - ancestor->level;
            int desc_nr = 0;

            for ( int r = ancestor->level + 1; r <= desc->level; ++r )
            {
                int co_num_cells_on_level = static_cast < int > ( std::pow ( static_cast < double > ( 8 ), diff_level - ( r - ancestor->level ) ) );
                int sibling_nr_on_level = pXest_z_order_to_ccw_order ( p8est_quadrant_ancestor_id ( desc, r ) );

                desc_nr += sibling_nr_on_level * co_num_cells_on_level;
            }
            return desc_nr;
        }

        // 3D_TODO check 

        int pXest_z_order_to_ccw_order ( int sibling_id )
        {
            switch ( sibling_id )
            {
                case 0:
                {
                    return 0;
                    break;
                }
                case 1:
                {
                    return 1;
                    break;
                }
                case 2:
                {
                    return 3;
                    break;
                }
                case 3:
                {
                    return 2;
                    break;
                }
                case 4:
                {
                    return 4;
                    break;
                }
                case 5:
                {
                    return 5;
                    break;
                }
                case 6:
                {
                    return 7;
                    break;
                }
                case 7:
                {
                    return 6;
                    break;
                }
            }
            return -1;
        }

        // 3D_TODO check

        int pXest_facet_number_pXest_to_db ( int facet )
        {
            switch ( facet )
            {
                case 0:
                {
                    return 3;
                    break;
                }
                case 1:
                {
                    return 1;
                    break;
                }
                case 2:
                {
                    return 0;
                    break;
                }
                case 3:
                {
                    return 2;
                    break;
                }
                case 4:
                {
                    return 4;
                    break;
                }
                case 5:
                {
                    return 5;
                    break;
                }
            }
        }

        // 3D_TODO check

        int pXest_facet_number_db_to_pXest ( int facet )
        {
            switch ( facet )
            {
                case 0:
                {
                    return 2;
                    break;
                }
                case 1:
                {
                    return 1;
                    break;
                }
                case 2:
                {
                    return 3;
                    break;
                }
                case 3:
                {
                    return 0;
                    break;
                }
                case 4:
                {
                    return 4;
                    break;
                }
                case 5:
                {
                    return 5;
                    break;
                }
            }
        }

        /// function that is called for each quadrant when initializing the forest, 
        /// allocates memory

        void pXest_init_fn ( p4est_t* forest, p4est_topidx_t tree, p4est_quadrant_t* quad )
        {
            // get pointer to external user data
            ForestInitData* i_ptr = ( ForestInitData * ) forest->user_pointer;

            // allocate memory for user data in quadrant
            QuadData* q_ptr = new QuadData ( tree, i_ptr->rank, 2 );

            // set correct pointer
            quad->p.user_data = q_ptr;
        }

        void pXest_init_fn ( p8est_t* forest, p4est_topidx_t tree, p8est_quadrant_t* quad )
        {
            // get pointer to external user data
            ForestInitData* i_ptr = ( ForestInitData * ) forest->user_pointer;

            // allocate memory for user data in quadrant
            QuadData* q_ptr = new QuadData ( tree, i_ptr->rank, 3 );

            // set correct pointer
            quad->p.user_data = q_ptr;
        }

#endif        
    }
}
