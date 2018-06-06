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

/// \author Philipp Gerstner

#ifndef PXEST_STRUCTS_H
#    define PXEST_STRUCTS_H

#    include <vector>
#    include <map>
#    include <iostream>
#    include <string>
#    include "common/sorted_array.h"
#    include "config.h"
#    include "mesh/types.h"
#    include <mpi.h>
#    include "communication.h"

#    ifdef WITH_P4EST
#        include "p4est.h"
#        include "p8est.h"
#    else
#        define P4EST_MAXLEVEL 1
#    endif

namespace hiflow
{
    namespace mesh
    {

        // Forward declaration
        class MeshPXestDatabase;
        typedef SharedPtr< MeshPXestDatabase >::Type MeshPXestDatabasePtr;

        struct QuadCoord;
        typedef int64_t mortonId;
        typedef int32_t treeId;
        typedef typename std::map < Id, QuadCoord > EntityToQuadMap;
        typedef typename std::vector< EntityToQuadMap > EntityToQuadMapping;

        /// \brief struct containing all information to find the corresponding p4est quadrant for a given hiflow entitty

        struct QuadCoord
        {
            // constructor

            QuadCoord ( )
            : tree ( -1 ), level ( -1 ), x ( -1 ), y ( -1 ), z ( -1 ), localId ( -1 ), tdim ( -1 ), coarseId ( -1 )
            {
            }

            QuadCoord ( int32_t t, int l, int32_t co_x, int32_t co_y, int32_t co_z, int Id, int dim, mesh::Id coarse_id )
            : tree ( t ), level ( l ), x ( co_x ), y ( co_y ), z ( co_z ), localId ( Id ), tdim ( dim ), coarseId ( coarse_id )
            {
            }

            // copy operator

            void copy ( QuadCoord& coord )
            {
                tree = coord.tree;
                level = coord.level;
                x = coord.x;
                y = coord.y;
                z = coord.z;
                localId = coord.localId;
                tdim = coord.tdim;
                coarseId = coord.coarseId;
            }

            // comparison operator, checks if two quad_maps point to the same quadrant in p4est forest

            bool is_same_quad ( QuadCoord& coord1, QuadCoord& coord2 )
            {
                if ( coord1.tree != coord2.tree )
                    return false;
                if ( coord1.level != coord2.level )
                    return false;
                if ( coord1.x != coord2.x )
                    return false;
                if ( coord1.y != coord2.y )
                    return false;
                if ( coord1.z != coord2.z )
                    return false;
                return true;
            }

            void print ( )
            {
                std::cout << "tree: " << tree << " level: " << level << " x: " << x << " y: " << y << " z: " << z
                        << " localId: " << localId << " tdim: " << tdim << " coarseId: " << coarseId << std::endl;
            }

            bool is_valid ( )
            {
                if ( tree < 0 )
                {
                    return false;
                }
                if ( level < 0 )
                {
                    return false;
                }
                if ( tdim < 0 )
                {
                    return false;
                }
                if ( x < 0 )
                {
                    return false;
                }
                if ( y < 0 )
                {
                    return false;
                }
                if ( z < 0 )
                {
                    return false;
                }
                return true;
            }

            // Id of tree
            int32_t tree;
            // Level of refinement
            int level;
            // Integer coordinates of corresponding p4est quadrant
            int32_t x;
            int32_t y;
            int32_t z;
            // Id of entity inside of corresponding quadrant. If entity is a cell: 0, facet:0-6, edge: 0-11, vertex:0-7
            int localId;
            // topological dimension of entity
            int tdim;
            // Id of cell in coarse mesh
            mesh::Id coarseId;

        };

        /// \brief struct containing external data used in p4est_init_fn

        struct ForestInitData
        {

            ForestInitData ( int mpi_rank )
            : rank ( mpi_rank )
            {
            }

            int rank;
        };

        /// \brief struct containing external data used in p4est_refine_init_fn

        struct ForestRefineData
        {

            ForestRefineData ( int mpi_rank, MeshPXestDatabasePtr db_ptr )
            : rank ( mpi_rank ), db ( db_ptr )
            {
            }

            void set_children_ids ( std::vector< std::vector<Id> > ids )
            {
                children_ids = ids;
            }

            int rank;
            std::vector< std::vector<Id> > children_ids;
            MeshPXestDatabasePtr db;
        };

        /// \brief struct containing external data used in p4est_coarsen_decide_fn

        struct ForestCoarsenData
        {
            // Ids of cells to be coarsened
            SortedArray<Id> coarsened_cells;
            // Ids of parent cells whose children were coarsened
            SortedArray<Id> parent_cells;
        };

        struct ForestPatchData
        {
            SortedArray<EntityNumber> coarsenable_cells;
        };

        /// \brief struct containing external data used in p4est_balance_replace_fn

        struct ForestBalanceData
        {
            SortedArray<EntityNumber> replaced_cell_indices;
            std::map<EntityNumber, int> ref_steps;
#    ifdef WITH_P4EST
            std::map<EntityNumber, std::vector<p4est_quadrant_t> > new_quads4;
            std::map<EntityNumber, std::vector<p8est_quadrant_t> > new_quads8;
#    endif
        };

        struct GhostCommData
        {
            EntityNumber cell_index;
        };

        /// \brief struct that contains that should be communicated with ghost cells. '22' -> tdim=2, gdim=2

        struct GhostCommData22
        {
            // tdim = 2, gdim = 2
            // #cell = 1, #facets = 4
            // num_vert(cell) = 4, num_vert(facet) = 2
            int num_cell;
            int num_facet;
            int num_vert_cell;
            int num_vert_facet;
            int connections [P4EST_MAXLEVEL * 12]; // size = #cells x num_vert(cell) + #facets x num_vert(facet) = 12
            double coords [P4EST_MAXLEVEL * 24]; // size = #cells x num_vert(cell)*gdim + #facetsx num_vert(facet)*gdim = 24
            int tdims [P4EST_MAXLEVEL * 5]; // size = #cells + #facets = 1 + 4 = 5
            int gdim;
            int offsets[P4EST_MAXLEVEL * 10]; // size = 2* (#cells + #facets) = 10
            int material_numbers[P4EST_MAXLEVEL * 5]; // size = #cells + #facets = 5

            int conn_offsets[P4EST_MAXLEVEL * 6]; // size = #cells + #facets + 1 = 6
            int coord_offsets[P4EST_MAXLEVEL * 6]; // size = #cells + #facets + 1 = 6

            EntityNumber cell_index;
            EntityNumber cell_index_remote;
            Id cell_id[P4EST_MAXLEVEL];
            Id cell_id_remote[P4EST_MAXLEVEL];
            int mesh_history[P4EST_MAXLEVEL];
            treeId tree_id;
            treeId tree_id_remote;
            int owner_rank;
            int tdim;
            bool active;
            bool refine;
            int coarsen;
        };

        struct GhostCommData33
        {
            // tdim = 3, gdim = 3
            // #cell = 1, #facets = 6
            // num_vert(cell) = 8, num_vert(facet) = 4
            int num_cell;
            int num_facet;
            int num_vert_cell;
            int num_vert_facet;
            int connections [P4EST_MAXLEVEL * 32]; // size = #cells x num_vert(cell) + #facets x num_vert(facet) = 32
            double coords [P4EST_MAXLEVEL * 96]; // size = #cells x num_vert(cell)*gdim + #facetsx num_vert(facet)*gdim = 24 + 72 = 96
            int tdims [P4EST_MAXLEVEL * 7]; // size = #cells + #facets = 1 + 6 = 7
            int gdim;
            int offsets[P4EST_MAXLEVEL * 14]; // size = 2* (#cells + #facets) = 14
            int material_numbers[P4EST_MAXLEVEL * 7]; // size = #cells + #facets = 7

            int conn_offsets[P4EST_MAXLEVEL * 8]; // size = #cells + #facets + 1 = 8
            int coord_offsets[P4EST_MAXLEVEL * 8]; // size = #cells + #facets + 1 = 8

            EntityNumber cell_index;
            EntityNumber cell_index_remote;
            Id cell_id[P4EST_MAXLEVEL];
            Id cell_id_remote[P4EST_MAXLEVEL];
            int mesh_history[P4EST_MAXLEVEL];
            treeId tree_id;
            treeId tree_id_remote;
            int owner_rank;
            int tdim;
            bool active;
            bool refine;
            int coarsen;
        };

        // PERF_TODO copy/set effizienter machen
        /// \brief struct that contains data which can be accessed by p4est_init_fn when initialzing the forest

        struct QuadData
        {

            QuadData ( )
            {
                init ( -1, -1, -1, -1, -1, -1, -1, -1, true, false, 0 );
            }

            QuadData ( treeId tree, int rank, int t_dim )
            {
                init ( -1, -1, -1, -1, tree, tree, rank, t_dim, true, false, 0 );
            }

            QuadData ( const QuadData& data )
            {
                set ( data );
            }

            void init ( Id id, Id id_remote, EntityNumber index, EntityNumber index_remote, treeId tree, treeId tree_remote, int rank, int t_dim,
                        bool status, bool ref, int coa )
            {
                cell_id[0] = id;
                cell_id_remote[0] = id_remote;
                for ( int l = 1; l < P4EST_MAXLEVEL; ++l )
                {
                    cell_id[l] = -1;
                    cell_id_remote[l] = -1;
                }
                cell_index = index;
                cell_index_remote = index_remote;
                tree_id = tree;
                tree_id_remote = tree_remote;
                owner_rank = rank;
                tdim = t_dim;
                active = status;
                refine = ref;
                coarsen = coa;
            }

            void set ( const QuadData& data )
            {
                for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
                {
                    cell_id[l] = data.cell_id[l];
                    cell_id_remote[l] = data.cell_id_remote[l];
                }
                cell_index_remote = data.cell_index_remote;
                cell_index = data.cell_index;
                tdim = data.tdim;
                owner_rank = data.owner_rank;
                tree_id = data.tree_id;
                tree_id_remote = data.tree_id_remote;
                active = data.active;
                refine = data.refine;
                coarsen = data.coarsen;
            }

            inline void set_cell_index ( EntityNumber index )
            {
                cell_index = index;
            }

            inline void set_remote_cell_index ( EntityNumber index )
            {
                cell_index_remote = index;
            }

            inline void set_cell_id ( Id id, int level )
            {
                cell_id[level] = id;
            }

            inline void set_remote_cell_id ( Id id, int level )
            {
                cell_id_remote[level] = id;
            }

            inline void set_active ( bool flag )
            {
                active = flag;
            }

            inline void set_refine ( bool flag )
            {
                refine = flag;
            }

            inline void set_coarsen ( int flag )
            {
                coarsen = flag;
            };

            inline void set_owner_rank ( int rank )
            {
                owner_rank = rank;
            };

            inline Id get_cell_id ( int level )
            {
                return this->cell_id[level];
            }

            inline Id get_remote_cell_id ( int level )
            {
                return this->cell_id_remote[level];
            }

            inline EntityNumber get_remote_cell_index ( )
            {
                return this->cell_index_remote;
            }

            inline EntityNumber get_cell_index ( )
            {
                return this->cell_index;
            }

            inline treeId get_tree_id ( )
            {
                return this->tree_id;
            }

            inline treeId get_remote_tree_id ( )
            {
                return this->tree_id;
            }

            inline int get_owner_rank ( )
            {
                return this->owner_rank;
            }

            inline bool is_active ( )
            {
                return this->active;
            }

            inline bool do_refine ( )
            {
                return this->refine;
            }

            inline int do_coarsen ( )
            {
                return this->coarsen;
            }

            void print ( )
            {
                std::cout << " owner rank: " << owner_rank << " is active: " << active << std::endl
                        << " tree_id: " << tree_id << " remote_tree_id: " << tree_id_remote << std::endl
                        << " cell index: " << cell_index << " remote cell_index: " << cell_index_remote << std::endl
                        << " tdim: " << tdim << std::endl
                        << " refine: " << refine << " coarsen: " << coarsen << std::endl;

                for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
                {
                    std::cout << cell_id[l] << " ";
                }
                std::cout << std::endl;
                for ( int l = 0; l < P4EST_MAXLEVEL; ++l )
                {
                    std::cout << cell_id_remote[l] << " ";
                }
                std::cout << std::endl << " ----------------------- " << std::endl;
            }

            // Index of cell corresponding to this quadrant, locally on process
            EntityNumber cell_index;

            // Id of cell corresponding to this quadrant, locally on process,
            // For interior quadrants: cell_id(proc) = cell_id_remote(proc)
            // For ghost quadrant: cell_id_remote(proc) = cell_id(owner proc)
            EntityNumber cell_index_remote;

            // Id of cell corresponding to this quadrant, locally on process
            Id cell_id[P4EST_MAXLEVEL];

            // Id of cell corresponding to this quadrant, locally on process,
            // For interior quadrants: cell_id(proc) = cell_id_remote(proc)
            // For ghost quadrant: cell_id_remote(proc) = cell_id(owner proc)
            Id cell_id_remote[P4EST_MAXLEVEL];

            // Id of tree
            treeId tree_id;
            treeId tree_id_remote;

            // Rank of owner process
            int owner_rank;

            // topological dimension of entity
            int tdim;

            // indicates whether quadrant is in current mesh, i.e. has no children
            bool active;

            // indicates whether quadrant should be refined
            bool refine;

            // indicates whether quadrant should be coarsened
            // if sum_{quad in family} quad.coarsen <= -1 * #{quad in family}, then family is coarsened
            int coarsen;
        };

        /// \brief return size of struct QuadData in bytes

        inline int QuadData_size ( int tdim )
        {
            int size = sizeof (QuadData );
            return size;
        }
    }
}
#endif
