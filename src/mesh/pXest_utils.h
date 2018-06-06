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

#ifndef PXEST_UTILS_H
#    define PXEST_UTILS_H

#    include <vector>
#    include <map>
#    include <iostream>
#    include <string>
#    include "common/sorted_array.h"
#    include "config.h"
#    include "mesh/types.h"
#    include <mpi.h>
#    include "communication.h"
#    include "pXest_structs.h"
#    include "refinement.h"

#    ifdef WITH_P4EST
#        include "p4est.h"
#        include "p8est.h"
#        include "p4est_ghost.h"
#        include "p8est_ghost.h"
#        include "p4est_connectivity.h"
#        include "p8est_connectivity.h"
#        include "p4est_bits.h"
#        include "p8est_bits.h"
#        include "p4est_iterate.h"
#        include "p8est_iterate.h"
#    endif

// TODO p4est namespace?
namespace hiflow
{
    namespace mesh
    {

        // BUG_TODO material number fehlt
        /// \brief functions for building ghost communication data structs
        /// @param[out] ret pointer to packed data strcut, ready for communicating
        /// @param[in] quad_data reference to data stored in p4est (mirror) quadrant
        /// @param[in] cell_packs cell entity packages to be communicated
        /// @param[in] facet_packs facet entity packages to be communicated
        void pack_ghost_data ( GhostCommData* ret, QuadData& quad_data );
        void pack_ghost_data ( GhostCommData22* ret, QuadData& quad_data );
        void pack_ghost_data ( GhostCommData22* ret, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs );
        void pack_ghost_data ( GhostCommData33* ret, QuadData& quad_data );
        void pack_ghost_data ( GhostCommData33* ret, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs );

        void unpack_ghost_data ( GhostCommData& ghost_data, QuadData& quad_data );
        void unpack_ghost_data ( GhostCommData22& ghost_data, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs );
        void unpack_ghost_data ( GhostCommData33& ghost_data, QuadData& quad_data, std::vector<EntityPackage>& cell_packs, std::vector<EntityPackage>& facet_packs );

        void print_ghost_data ( GhostCommData22* ret );
        void print_ghost_data ( GhostCommData33* ret );
        void print_pack ( EntityPackage& pack );

#    ifdef WITH_P4EST
        /// \brief function that makes a deep copy of a given p4est ghost struct
        /// @param[in] ghost struct to copy from
        /// @param[out]  pointer to newly allocated ghost
        void pXest_copy_ghost ( p4est_ghost_t* input, p4est_ghost_t*& output );
        void pXest_copy_ghost ( p8est_ghost_t* input, p8est_ghost_t*& output );

        /// \brief function that builds a p4est connectivity out of read in data structures
        /// @param[in] num_vertices number of vertices in coarse mesh
        /// @param[in] num_cells number of cells in coarse mesh
        /// @param[in] vertices coordinates of vertices
        /// @param[in] tree_to_vertices cell->vertex connectivity
        /// @param[out] conn pointer to newly created connectivity struct
        void pXest_build_conn ( int num_vertices, int num_cells, const std::vector<double>& vertices, const std::vector<int>& tree_to_vertices,
                                p4est_connectivity_t*& conn );
        void pXest_build_conn ( int num_vertices, int num_cells, const std::vector<double>& vertices, const std::vector<int>& tree_to_vertices,
                                p8est_connectivity_t*& conn );

        /// \brief reoder connectivity struct such that forest is distributed according to given partition
        /// @param[in] part partitioning of coarse mesh
        /// @param[in/out] conn connectivity to be permuted
        /// @param[out] perm permutation of trees
        void pXest_reorder_conn ( const std::vector<int>& part, p4est_connectivity_t * conn, std::vector<size_t>& perm );
        void pXest_reorder_conn ( const std::vector<int>& part, p8est_connectivity_t * conn, std::vector<size_t>& perm );

        p4est_gloidx_t pXest_partition_forest ( p4est_t * forest, int partition_for_coarsening, std::vector<p4est_locidx_t>& num_quads_per_proc );
        p4est_gloidx_t pXest_partition_forest ( p8est_t * forest, int partition_for_coarsening, std::vector<p4est_locidx_t>& num_quads_per_proc );

        /// \brief extract mesh partition from unrefined forest
        /// @param[in] forest pointer to forest object
        /// @param[in] master_rank
        /// @return map: coarse cell id -> MPI rank
        std::vector<int> pXest_get_partition_from_initial_forest ( p4est_t* forest );
        std::vector<int> pXest_get_partition_from_initial_forest ( p8est_t* forest );

        /// \brief extract mesh partition from refined forest -> consider childless quadrants (leafs) as current mesh
        /// @param[in] forest pointer to forest object
        /// @param[in] master_rank
        /// @return map: coarse cell id -> MPI rank
        std::vector<int> pXest_get_partition_from_forest ( p4est_t* forest );
        std::vector<int> pXest_get_partition_from_forest ( p8est_t* forest );

        /// \brief broadcast coarse entity_to_quad maps
        /// @param[in] comm MPI communicator
        /// @param[in] master_rank
        /// @param[in] partitioning previuously created partition map
        /// @param[in] coarse_cell_ids ids of cells in coarse mesh
        /// @param[in] local_cell_ids ids of distributed cells in local mesh
        /// @param[in] coarse_maps vector of coarse_entity_to_quad maps, must be nonempty iff on master process
        /// @param[out] maps distributed entity_to_quad maps with adjusted id (now local id instead of global coarse id)
        void pXest_distribute_coarse_maps ( MPI_Comm comm, int master_rank,
                                            std::vector<int>& partitioning, std::vector<int>& coarse_cell_ids, std::vector<int>& local_cell_ids,
                                            EntityToQuadMap* coarse_map, EntityToQuadMap& map );

        /// \brief extract local entity_to_quad maps from set of coarse maps
        /// @param[in] comm MPI communicator
        /// @param[in] master_rank
        /// @param[in] partitioning previuously created partition map
        /// @param[in] coarse_cell_ids ids of cells in coarse mesh
        /// @param[in] local_cell_ids ids of distributed cells in local mesh
        /// @param[in] coarse_maps vector of coarse_entity_to_quad maps, must be nonempty iff on master process
        /// @param[out] maps distributed entity_to_quad maps with adjusted id (now local id instead of global coarse id)
        void pXest_extract_local_maps ( MPI_Comm comm, int master_rank,
                                        std::vector<int>& partitioning, std::vector<int>& coarse_cell_ids, std::vector<int>& local_cell_ids,
                                        EntityToQuadMap* coarse_map, EntityToQuadMap& map );

        /// \brief get tree id of quadrant
        /// @param[in] quad pointer to some quadrant
        /// @return id of tree containing quadrant
        treeId pXest_get_treeId ( p4est_quadrant_t* quad );
        treeId pXest_get_treeId ( p8est_quadrant_t* quad );

        /// \brief creates a quadrant object out of coordinates describing its positin in the forest
        /// @param[in] QuadCoord coordinates describing position of quadrant in forest
        /// @return quadrant object (without user data)
        p4est_quadrant_t pXest_create_quad4_from_coord ( QuadCoord );
        p8est_quadrant_t pXest_create_quad8_from_coord ( QuadCoord );

        /// \brief get Morton index of quadrant
        /// @param[in] quad pointer to some quadrant
        /// @param[in] level level of refinement in tree
        /// @return morton index of quadrant on given level
        mortonId pXest_get_mortonId ( p4est_quadrant_t* quad, int level );
        mortonId pXest_get_mortonId ( p8est_quadrant_t* quad, int level );

        /// \brief compute morton index of last quadrant stored on given tree
        /// @param[in] tree pointer to tree object
        /// @param[in] level
        /// @return morton index of last quadrant on given level
        mortonId pXest_get_last_mortonId_in_tree ( p4est_tree_t* tree, int level );
        mortonId pXest_get_last_mortonId_in_tree ( p8est_tree_t* tree, int level );

        /// \brief compute the position in quadrant array of tree object for some given quad
        /// @param[in] pointer to tree object
        /// @param[in] quad pointer to some quadrant
        /// @return index of quadrant in tree array
        int64_t pXest_get_quad_pos_in_tree_array ( p4est_tree_t* tree, p4est_quadrant_t* quad );
        int64_t pXest_get_quad_pos_in_tree_array ( p8est_tree_t* tree, p8est_quadrant_t* quad );

        /// \brief compute the position in quadrant array of ghost object for some given quad
        /// @param[in] ghost pointer to ghost object
        /// @param[in] quad pointer to some quadrant
        /// @param[in] tree_id global id of tree of quadrant
        /// @return index of quadrant in ghost array
        int64_t pXest_get_quad_pos_in_ghost_array ( p4est_ghost_t* ghost, p4est_quadrant_t* quad, treeId tree_id );
        int64_t pXest_get_quad_pos_in_ghost_array ( p8est_ghost_t* ghost, p8est_quadrant_t* quad, treeId tree_id );

        /// \brief  compute positions in quadrant array of some tree for those quadrants that do not have descendants
        /// @param[in] pointer to tree object
        /// @param[in] quad pointer to some quadrant
        /// @return indices of leaf quadrants in tree array
        std::vector<int64_t> pXest_get_leaf_desc_in_tree_array ( p4est_tree_t* tree, p4est_quadrant_t* quad );
        std::vector<int64_t> pXest_get_leaf_desc_in_tree_array ( p8est_tree_t* tree, p8est_quadrant_t* quad );

        /// \brief compute positions in quadrant array of ghost for those quadrants that do not have descendants
        /// @param[in] ghost pointer to ghost object
        /// @param[in] quad pointer to some quadrant
        /// @param[in] tree_id global id of tree
        /// @return indices of leaf quadrants in tree array
        std::vector<int64_t> pXest_get_leaf_desc_in_ghost_array ( p4est_ghost_t* ghost, p4est_quadrant_t* quad, treeId tree_id );
        std::vector<int64_t> pXest_get_leaf_desc_in_ghost_array ( p8est_ghost_t* ghost, p8est_quadrant_t* quad, treeId tree_id );

        /// \brief check if quadrant arrays in trees of given forest are sorted
        /// @param[in] forest pointer to underlying forest object
        /// @return true if quadrant arrays are sorted
        bool pXest_quad_arrays_in_forest_sorted ( p4est_t* forest );
        bool pXest_quad_arrays_in_forest_sorted ( p8est_t* forest );

        /// \brief sort quadrant arrays in trees of given forest
        /// @param[in] forest pointer to underlying forest object
        void pXest_sort_quad_arrays_in_forest ( p4est_t* forest );
        void pXest_sort_quad_arrays_in_forest ( p8est_t* forest );

        /// \brief set QuadData in p4est quadrant
        /// @param[in] quad pointer to some quadrant
        /// @param[in] data reference to user quadrant data struct
        void pXest_set_quad_data ( p4est_quadrant_t* quad, QuadData& data );
        void pXest_set_quad_data ( p8est_quadrant_t* quad, QuadData& data );

        /// \brief get QuadData in p4est quadrant. Note: not all quadrants contain user data
        /// @param[in] quad pointer to some quadrant
        /// @param[in] data reference to user quadrant data struct
        void pXest_get_quad_data ( p4est_quadrant_t* quad, QuadData& data );
        void pXest_get_quad_data ( p8est_quadrant_t* quad, QuadData& data );

        /// \brief get pointer to data of quadrant in given tree. Here, quadrant denotes an empty quadrant, i.e. without user data
        /// @param[in] tree pointer to tree object
        /// @param[in] quad pointer to some (empty) quadrant
        /// @return pointer to quadrant user data
        QuadData* pXest_get_quad_data_ptr ( p4est_tree_t* tree, p4est_quadrant_t* quad );
        QuadData* pXest_get_quad_data_ptr ( p8est_tree_t* tree, p8est_quadrant_t* quad );

        /// \brief get pointer to data of quadrant in given ghost layer. Here, quadrant denotes an empty quadrant, i.e. without user data
        /// @param[in] ghost pointer to ghost object
        /// @param[in] quad pointer to some (empty) quadrant
        /// @param[in] tree_id (optional) id of tree to which wuadrant belongs.
        /// @return pointer to quadrant user data
        QuadData* pXest_get_quad_data_ptr ( p4est_ghost_t* ghost, p4est_quadrant_t* quad, treeId tree_id );
        QuadData* pXest_get_quad_data_ptr ( p8est_ghost_t* ghost, p8est_quadrant_t* quad, treeId tree_id );

        /// \brief compute QuadCoord from quadrant
        /// @param[in] quad pointer to some quadrant
        /// @param[in] tree_id global id of tree
        /// @return quadrant coordinates object
        QuadCoord pXest_get_quad_coord ( p4est_quadrant_t* quad, treeId tree_id );
        QuadCoord pXest_get_quad_coord ( p8est_quadrant_t* quad, treeId tree_id );

        /// \brief build a RefinementTree object whose structure matches the relationship of the given quadrants
        /// @param[in] ancestor quadrant corresponding to root cell
        /// @param[in] descendants quadrants corresponding to descendants of root cell
        /// @param[out] tree initialized RefinementTree. sub_cell_numbers in tree correspond to array indices+1
        void pXest_build_ref_tree ( std::vector<p4est_quadrant_t>& quads, RefinementTree* tree );
        void pXest_build_ref_tree ( std::vector<p8est_quadrant_t>& quads, RefinementTree* tree );

        // PERF_TODO more efficient solution possible?
        /// \brief convert p4est family index to hiflow family index for a refined cell.
        /// \@param[in] p4est sibling id
        /// @return index of cell in counter clockwise ordering
        int pXest_z_order_to_ccw_order ( int sibling_id );

        /// \brief get sibling id of quad, i.e. the position in the z curve within a family created by refining one single cell
        /// @param[in] quad pointer to some quadrant
        /// @return sibling number

        inline int pXest_get_sibling_nr ( const p4est_quadrant_t* quad )
        {
            return p4est_quadrant_child_id ( quad );
        }

        inline int pXest_get_sibling_nr ( const p8est_quadrant* quad )
        {
            return p8est_quadrant_child_id ( quad );
        }

        /// \brief
        /// @param[in] desc pointer to descendant quadrant
        /// @param[in] ancestor pointer to ancestor quadrant
        /// @return
        int pXest_get_db_desc_nr ( p4est_quadrant_t* desc, p4est_quadrant_t* ancestor );
        int pXest_get_db_desc_nr ( p8est_quadrant_t* desc, p8est_quadrant_t* ancestor );

        /// \brief convert local facet number from p4est numbering to mesh database numbering
        /// @param[in] p4est number of facet
        /// @return MeshDBView number of facet
        int pXest_facet_number_pXest_to_db ( int facet );

        /// \brief convert local facet number from mesh database numbering to p4est numbering
        /// @param[in] MeshDBView number of facet
        /// @return p4est number of facet
        int pXest_facet_number_db_to_pXest ( int facet );

        /// \brief function which is called for each quadrant when creating a new p4est forest object
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of tree
        /// @param[in] quad pointer to some quadrant
        void pXest_init_fn ( p4est_t* forest, p4est_topidx_t tree, p4est_quadrant_t* quad );
        void pXest_init_fn ( p8est_t* forest, p4est_topidx_t tree, p8est_quadrant_t* quad );

        /// \brief function which is called for each quadrant when partitioning a forest
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of tree
        /// @param[in] quad pointer to some quadrant
        /// @return weight of quadrant
        int pXest_weight_fn ( p4est_t* forest, p4est_topidx_t tree, p4est_quadrant_t* quad );
        int pXest_weight_fn ( p8est_t* forest, p4est_topidx_t tree, p8est_quadrant_t* quad );

        /// \brief get id of non-empty first tree on processor
        /// @param[in] forest pointer to underlying forest object
        /// @return id of tree

        inline treeId pXest_get_first_local_treeId ( p4est_t* forest )
        {
            return forest->first_local_tree;
        }

        inline treeId pXest_get_first_local_treeId ( p8est_t* forest )
        {
            return forest->first_local_tree;
        }

        /// \brief get id of last non-empty tree on processor
        /// @param[in] forest pointer to underlying forest object
        /// @return id of tree

        inline treeId pXest_get_last_local_treeId ( p4est_t* forest )
        {
            return forest->last_local_tree;
        }

        inline treeId pXest_get_last_local_treeId ( p8est_t* forest )
        {
            return forest->last_local_tree;
        }

        /// \brief get pointer to specific tree in forest
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of tree
        /// @return pointer to tree object

        inline p4est_tree_t* pXest_get_tree_in_forest ( p4est_t* forest, treeId tree_id )
        {
            return p4est_tree_array_index ( forest->trees, tree_id );
        }

        inline p8est_tree_t* pXest_get_tree_in_forest ( p8est_t* forest, treeId tree_id )
        {
            return p8est_tree_array_index ( forest->trees, tree_id );
        }

        /// \brief check if some quadrant is ancestor of another one
        /// @param[in] test pointer to potential ancestor quadrant
        /// @param[in] desc pointer to descendant quadrant
        /// @return true if test is ancestor of desc

        inline bool pXest_quad_is_ancestor ( p4est_quadrant_t* test, p4est_quadrant_t* desc )
        {
            return p4est_quadrant_is_ancestor ( test, desc );
        }

        inline bool pXest_quad_is_ancestor ( p8est_quadrant_t* test, p8est_quadrant_t* desc )
        {
            return p8est_quadrant_is_ancestor ( test, desc );
        }

        /// \brief set user data pointer in forest
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] ptr pointer to user data struct

        inline void pXest_set_forest_user_ptr ( p4est_t* forest, void * ptr )
        {
            forest->user_pointer = ptr;
        }

        inline void pXest_set_forest_user_ptr ( p8est_t* forest, void * ptr )
        {
            forest->user_pointer = ptr;
        }

        /// \brief get pointer to QuadData in p4est quadrant
        /// Note: not all quadrants contain user data
        /// @param[in] quad pointer to some quadrant
        /// @return pointer to quadrant user data

        inline QuadData* pXest_get_quad_data_ptr ( p4est_quadrant_t* quad )
        {
            return ( QuadData* ) quad->p.user_data;
        }

        inline QuadData* pXest_get_quad_data_ptr ( p8est_quadrant_t* quad )
        {
            return ( QuadData* ) quad->p.user_data;
        }

        /// \brief set user data pointer in quadrant
        /// @param[in] quad pointer to some quadrant
        /// @param[in] q_ptr pointer to user quadrant data struct

        inline void pXest_set_quad_data_ptr ( p4est_quadrant_t* quad, QuadData* q_ptr )
        {
            quad->p.user_data = q_ptr;
        }

        inline void pXest_set_quad_data_ptr ( p8est_quadrant_t* quad, QuadData* q_ptr )
        {
            quad->p.user_data = q_ptr;
        }

        /// \brief access a quadrant in the tree quadrant array
        /// @param[in] tree pointer to tree object
        /// @param[in] pos index of quadrant in tree array
        /// @return pointer to pos-th quadrant in given tree

        inline p4est_quadrant_t* pXest_get_local_quad_in_tree ( p4est_tree_t* tree, int64_t pos )
        {
            return p4est_quadrant_array_index ( &tree->quadrants, pos );
        }

        inline p8est_quadrant_t* pXest_get_local_quad_in_tree ( p8est_tree_t* tree, int64_t pos )
        {
            return p8est_quadrant_array_index ( &tree->quadrants, pos );
        }

        /// \brief access a quadrant in the ghost quadrant array
        /// @param[in] ghost pointer to ghost object
        /// @param[in] index of quadrant in ghost array
        /// @return pointer to pos-th quadrant in given ghost layer

        inline p4est_quadrant_t* pXest_get_local_quad_in_ghost ( p4est_ghost_t* ghost, int64_t pos )
        {
            return p4est_quadrant_array_index ( &ghost->ghosts, pos );
        }

        inline p8est_quadrant_t* pXest_get_local_quad_in_ghost ( p8est_ghost_t* ghost, int64_t pos )
        {
            return p8est_quadrant_array_index ( &ghost->ghosts, pos );
        }

        /// \brief compare two quadrants w.r.t. their morton index
        /// @param[in] v1 pointer to some quadrant
        /// @param[in] v2 pointer to some quadrant
        /// @return -1 if v1 < v2, 0 if v1 = v2, 1 if v1 > v2

        inline int pXest_compare_quad4_in_tree ( const void* q1, const void* q2 )
        {
            const p4est_quadrant_t *v1 = ( const p4est_quadrant_t * ) q1;
            const p4est_quadrant_t *v2 = ( const p4est_quadrant_t * ) q2;
            return p4est_quadrant_compare ( v1, v2 );
        }

        inline int pXest_compare_quad8_in_tree ( const void* q1, const void* q2 )
        {
            const p8est_quadrant_t *v1 = ( const p8est_quadrant_t * ) q1;
            const p8est_quadrant_t *v2 = ( const p8est_quadrant_t * ) q2;
            return p8est_quadrant_compare ( v1, v2 );
        }

        /// \brief compute morton index of some quadrant w.r.t. the current level of the quadrant
        /// @param[in] quad pointer to some quadrant
        /// @return morton index

        inline mortonId pXest_get_mortonId ( p4est_quadrant_t* quad )
        {
            return p4est_quadrant_linear_id ( quad, quad->level );
        }

        inline mortonId pXest_get_mortonId ( p8est_quadrant_t* quad )
        {
            return p8est_quadrant_linear_id ( quad, quad->level );
        }

        /// \brief get Morton Index of first and last local quadrant on tree for given level of refinement
        /// @param[in] tree pointer to tree object
        /// @param[in] level level of refinement
        /// @return morton index

        inline mortonId pXest_get_first_mortonId_in_tree ( p4est_tree_t* tree, int level )
        {
            return p4est_quadrant_linear_id ( &tree->first_desc, level );
        }

        inline mortonId pXest_get_first_mortonId_in_tree ( p8est_tree_t* tree, int level )
        {
            return p8est_quadrant_linear_id ( &tree->first_desc, level );
        }
#    endif
    }
}
#endif
