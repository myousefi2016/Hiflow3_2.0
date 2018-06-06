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

#include "mesh_pXest.h"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include "common/log.h"

#include "cell_type.h"
#include "communication.h"
#include "connectivity.h"
#include "entity.h"
#include "iterator.h"
#include "mpi_communication.h"
#include "refinement.h"
#include "interface.h"

#ifdef WITH_P4EST
#    include "p4est_iterate.h"
#    include "p8est_iterate.h"
#    include "p4est_extended.h"
#    include "p8est_extended.h"
#    include "sc.h"
#endif

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    namespace mesh
    {
        //////////// Helper ////////////////////////////////////////////////////
#ifdef WITH_P4EST

        /// \brief function which is called for each facet when iterating through a given forest
        /// @param[in] info pointer to interface object
        /// @param[in] user_data pointer to user data

        void pXest_interface_list_fn ( p4est_iter_face_info_t *info, void *user_data )
        {
            // get pointer to interface list
            InterfaceList* list = ( InterfaceList* ) user_data;

            // get information of current facet
            int8_t orientation = info->orientation;
            int8_t is_tree_boundary = info->tree_boundary;
            sc_array_t *sides = &( info->sides );
            int num_sides = sides->elem_count;
            int num_hanging_cells = 2;

            assert ( num_sides <= 2 );

            /* TODO
            if (tdim == 3)
            {
                num_hanging_cells = 4;
            }
             * */
            LOG_DEBUG ( 1, " interface: num_sides " << num_sides << " orientation " << static_cast < int16_t > ( orientation )
                        << " is_tree_boundary " << static_cast < int16_t > ( is_tree_boundary ) );

            // collect data from adjacent quadrants
            std::vector<int> is_hanging;
            std::vector<int> pXest_facet_nr;
            std::vector< std::vector<p4est_quadrant_t*> > quads;
            std::vector< std::vector<int> > is_ghost;
            std::vector<p4est_iter_face_side_t*> side;

            side.resize ( num_sides );
            pXest_facet_nr.resize ( num_sides, -1 );
            is_hanging.resize ( num_sides, 0 );
            quads.resize ( num_sides );
            is_ghost.resize ( num_sides );

            // Loop over all sides corresponding to current facet
            for ( int js = 0; js < num_sides; ++js )
            {
                side[js] = p4est_iter_fside_array_index_int ( sides, js );
                pXest_facet_nr[js] = static_cast < int16_t > ( side[js]->face );
                is_hanging[js] = static_cast < int16_t > ( side[js]->is_hanging );

                // Get pointers to quadrant adjacent to current side          
                if ( is_hanging[js] )
                {
                    // facet contains hanging node 
                    quads[js].resize ( num_hanging_cells );
                    is_ghost[js].resize ( num_hanging_cells, 0 );
                    for ( int jq = 0; jq < num_hanging_cells; jq++ )
                    {
                        quads[js][jq] = side[js]->is.hanging.quad[jq];
                        is_ghost[js][jq] = static_cast < int16_t > ( side[js]->is.hanging.is_ghost[jq] );

                        // ****
                        /*
                        p4est_locidx_t quadid = side[js]->is.hanging.quadid[jq];
                        QuadData* q_ptr = pXest_get_quad_data_ptr(quad);
                        std::cout << "   Quad " << jq << " in current side has local id " << quadid 
                                  << " and ghost flag " << is_ghost[js][jq] << std::endl;
                        std::cout << "   ";
                        q_ptr->print();
                        std::cout << " -- " << std::endl;
                         * */
                        // ****
                    }
                }
                else
                {
                    quads[js].resize ( 1 );
                    is_ghost[js].resize ( 1, 0 );

                    quads[js][0] = side[js]->is.full.quad;
                    is_ghost[js][0] = static_cast < int16_t > ( side[js]->is.full.is_ghost );

                    // ****    
                    /*
                    p4est_locidx_t quadid = side[js]->is.full.quadid;
                    QuadData* q_ptr = pXest_get_quad_data_ptr(quad);

                    std::cout << "   Quad  in current side has local id " << quadid 
                              << " and ghost flag " << is_ghost[js][0] << std::endl;
                    std::cout << "   ";
                    q_ptr->print();
                    std::cout << " -- " << std::endl;
                     * */
                    // ****
                }
            }

            // data to be stored in interface list
            EntityNumber master_index;
            std::vector<EntityNumber> slave_indices;
            int master_facet_number;
            std::vector<int> slave_facet_numbers;
            //std::tr1::unordered_map<EntityNumber, int> irregular_interfaces;

            // check type of facet: 
            // num_sides = 1 -> boundary
            // num_sides = 2:
            //     is_hanging = 0: standard interior facet
            //     is_hanging = 1: interior facet with hanging node    

            if ( num_sides == 1 ) // Add boundary facet
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[0][0] );
                master_index = q_ptr->get_cell_index ( );
                int master_facet_number = pXest_facet_number_pXest_to_db ( pXest_facet_nr[0] );

                list->add_interface ( master_index, master_facet_number );
                LOG_DEBUG ( 1, "side " << 0 << " tree id " << side[0]->treeid
                            << " face (db)" << master_facet_number
                            << " is hanging " << is_hanging[0] );
                LOG_DEBUG ( 1, "(Boundary facet) Master cell index: " << master_index );
            }
            else
            {
                assert ( is_hanging[0] + is_hanging[1] < 2 );

                if ( is_hanging[0] + is_hanging[1] == 0 ) // standard interior facet
                {
                    std::vector<QuadData*> q_ptr ( 2 );
                    std::vector<EntityNumber> cell_indices ( 2, -1 );
                    std::vector<int> db_facet_nr ( 2, -1 );

                    for ( int js = 0; js < 2; ++js )
                    {
                        q_ptr[js] = pXest_get_quad_data_ptr ( quads[js][0] );
                        cell_indices[js] = q_ptr[js]->get_cell_index ( );
                        db_facet_nr[js] = pXest_facet_number_pXest_to_db ( pXest_facet_nr[js] );

                        LOG_DEBUG ( 1, "side " << js << " tree id " << side[js]->treeid
                                    << " face (db) " << db_facet_nr[js]
                                    << " is hanging " << is_hanging[js] );
                    }

                    const int master = ( db_facet_nr[0] < db_facet_nr[1] ) ? 0 : 1;
                    const int slave = 1 - master;
                    Interface& interface = list->add_interface ( cell_indices[master], db_facet_nr[master] );
                    interface.add_slave ( cell_indices[slave], db_facet_nr[slave] );

                    LOG_DEBUG ( 1, "(Regular facet) Master cell index: " << cell_indices[master] );
                    LOG_DEBUG ( 1, "(Regular facet) Slave cell index: " << cell_indices[slave] );
                }
                else // facet with hanging node
                {
                    const int master_side = ( is_hanging[0] == 0 ) ? 0 : 1;
                    const int slave_side = 1 - master_side;

                    // add master cell
                    QuadData* m_ptr = pXest_get_quad_data_ptr ( quads[master_side][0] );
                    EntityNumber master_index = m_ptr->get_cell_index ( );
                    int master_facet_number = pXest_facet_number_pXest_to_db ( pXest_facet_nr[master_side] );

                    Interface& interface = list->add_interface ( master_index, master_facet_number );
                    LOG_DEBUG ( 1, "(Hanging facet) Master cell index: " << master_index );

                    // add slave cells
                    for ( int jq = 0; jq < num_hanging_cells; ++jq )
                    {
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[slave_side][jq] );
                        EntityNumber slave_index = q_ptr->get_cell_index ( );
                        int slave_facet_number = pXest_facet_number_pXest_to_db ( pXest_facet_nr[slave_side] );
                        interface.add_slave ( slave_index, slave_facet_number );
                        LOG_DEBUG ( 1, "(Hanging facet) Slave cell index: " << slave_index );
                    }
                }
            }
        }

        /// \brief function which is called for each facet when iterating through a given forest
        /// @param[in] info pointer to interface object
        /// @param[in] user_data pointer to user data

        void pXest_interface_list_fn ( p8est_iter_face_info_t *info, void *user_data )
        {
            // get pointer to interface list
            InterfaceList* list = ( InterfaceList* ) user_data;

            // get information of current facet
            int8_t orientation = info->orientation;
            int8_t is_tree_boundary = info->tree_boundary;
            sc_array_t *sides = &( info->sides );
            int num_sides = sides->elem_count;
            int num_hanging_cells = 4;

            assert ( num_sides <= 2 );

            LOG_DEBUG ( 1, " interface: num_sides " << num_sides << " orientation " << static_cast < int16_t > ( orientation )
                        << " is_tree_boundary " << static_cast < int16_t > ( is_tree_boundary ) );

            // collect data from adjacent quadrants
            std::vector<int> is_hanging;
            std::vector<int> pXest_facet_nr;
            std::vector< std::vector<p8est_quadrant_t*> > quads;
            std::vector< std::vector<int> > is_ghost;
            std::vector<p8est_iter_face_side_t*> side;

            side.resize ( num_sides );
            pXest_facet_nr.resize ( num_sides, -1 );
            is_hanging.resize ( num_sides, 0 );
            quads.resize ( num_sides );
            is_ghost.resize ( num_sides );

            // Loop over all sides corresponding to current facet
            for ( int js = 0; js < num_sides; ++js )
            {
                side[js] = p8est_iter_fside_array_index_int ( sides, js );
                pXest_facet_nr[js] = static_cast < int16_t > ( side[js]->face );
                is_hanging[js] = static_cast < int16_t > ( side[js]->is_hanging );

                // Get pointers to quadrant adjacent to current side          
                if ( is_hanging[js] )
                {
                    // facet contains hanging node 
                    quads[js].resize ( num_hanging_cells );
                    is_ghost[js].resize ( num_hanging_cells, 0 );
                    for ( int jq = 0; jq < num_hanging_cells; jq++ )
                    {
                        quads[js][jq] = side[js]->is.hanging.quad[jq];
                        is_ghost[js][jq] = static_cast < int16_t > ( side[js]->is.hanging.is_ghost[jq] );
                    }
                }
                else
                {
                    quads[js].resize ( 1 );
                    is_ghost[js].resize ( 1, 0 );

                    quads[js][0] = side[js]->is.full.quad;
                    is_ghost[js][0] = static_cast < int16_t > ( side[js]->is.full.is_ghost );
                }
            }

            // data to be stored in interface list
            EntityNumber master_index;
            std::vector<EntityNumber> slave_indices;
            int master_facet_number;
            std::vector<int> slave_facet_numbers;
            //std::tr1::unordered_map<EntityNumber, int> irregular_interfaces;

            // check type of facet: 
            // num_sides = 1 -> boundary
            // num_sides = 2:
            //     is_hanging = 0: standard interior facet
            //     is_hanging = 1: interior facet with hanging node    

            if ( num_sides == 1 ) // Add boundary facet
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[0][0] );
                master_index = q_ptr->get_cell_index ( );
                int master_facet_number = pXest_facet_number_pXest_to_db ( pXest_facet_nr[0] );

                list->add_interface ( master_index, master_facet_number );
                LOG_DEBUG ( 1, "side " << 0 << " tree id " << side[0]->treeid
                            << " face (db)" << master_facet_number
                            << " is hanging " << is_hanging[0] );
                LOG_DEBUG ( 1, "(Boundary facet) Master cell index: " << master_index );
            }
            else
            {
                assert ( is_hanging[0] + is_hanging[1] < 2 );

                if ( is_hanging[0] + is_hanging[1] == 0 ) // standard interior facet
                {
                    std::vector<QuadData*> q_ptr ( 2 );
                    std::vector<EntityNumber> cell_indices ( 2, -1 );
                    std::vector<int> db_facet_nr ( 2, -1 );

                    for ( int js = 0; js < 2; ++js )
                    {
                        q_ptr[js] = pXest_get_quad_data_ptr ( quads[js][0] );
                        cell_indices[js] = q_ptr[js]->get_cell_index ( );
                        db_facet_nr[js] = pXest_facet_number_pXest_to_db ( pXest_facet_nr[js] );

                        LOG_DEBUG ( 1, "side " << js << " tree id " << side[js]->treeid
                                    << " face (db) " << db_facet_nr[js]
                                    << " is hanging " << is_hanging[js] );
                    }

                    const int master = ( db_facet_nr[0] < db_facet_nr[1] ) ? 0 : 1;
                    const int slave = 1 - master;
                    Interface& interface = list->add_interface ( cell_indices[master], db_facet_nr[master] );
                    interface.add_slave ( cell_indices[slave], db_facet_nr[slave] );

                    LOG_DEBUG ( 1, "(Regular facet) Master cell index: " << cell_indices[master] );
                    LOG_DEBUG ( 1, "(Regular facet) Slave cell index: " << cell_indices[slave] );
                }
                else // facet with hanging node
                {
                    const int master_side = ( is_hanging[0] == 0 ) ? 0 : 1;
                    const int slave_side = 1 - master_side;

                    // add master cell
                    QuadData* m_ptr = pXest_get_quad_data_ptr ( quads[master_side][0] );
                    EntityNumber master_index = m_ptr->get_cell_index ( );
                    int master_facet_number = pXest_facet_number_pXest_to_db ( pXest_facet_nr[master_side] );

                    Interface& interface = list->add_interface ( master_index, master_facet_number );
                    LOG_DEBUG ( 1, "(Hanging facet) Master cell index: " << master_index );

                    // add slave cells
                    for ( int jq = 0; jq < num_hanging_cells; ++jq )
                    {
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[slave_side][jq] );
                        EntityNumber slave_index = q_ptr->get_cell_index ( );
                        int slave_facet_number = pXest_facet_number_pXest_to_db ( pXest_facet_nr[slave_side] );
                        interface.add_slave ( slave_index, slave_facet_number );
                        LOG_DEBUG ( 1, "(Hanging facet) Slave cell index: " << slave_index );
                    }
                }
            }
        }

        /// \brief function which is called on each local quadrant when calling p4est_refine().
        /// The return value of this function decides whether the quadrant should be refined
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of current tree object
        /// @param[in] quad pointer to current quadrant object
        /// @return true if quad should be refined

        int pXest_refine_decide_fn ( p4est_t* forest, p4est_topidx_t tree, p4est_quadrant_t* quad )
        {
            // return flag indicating whether to refine or not
            QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );

            bool refine = q_ptr->do_refine ( );
            q_ptr->set_refine ( false );
            return refine;
        }

        /// \brief function which is called on each local quadrant when calling p4est_refine().
        /// The return value of this function decides whether the quadrant should be refined
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of current tree object
        /// @param[in] quad pointer to current quadrant object
        /// @return true if quad should be refined

        int pXest_refine_decide_fn ( p8est_t* forest, p4est_topidx_t tree, p8est_quadrant_t* quad )
        {
            // return flag indicating whether to refine or not
            QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );

            bool refine = q_ptr->do_refine ( );
            q_ptr->set_refine ( false );
            return refine;
        }

        /// \brief function that is called for each children quadrant which is generated when calling 
        /// p4est_refine(). The purpose of this function is to set initialize the QuadData structure for 
        /// the newly created quadrants.
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of current tree object
        /// @param[in] quad pointer to current quadrant object

        void pXest_refine_init_fn ( p4est_t* forest, p4est_topidx_t tree_id, p4est_quadrant_t* quad )
        {
            // initialize quad data for new quadrant       
            QuadData* q_ptr = new QuadData ( );
            quad->p.user_data = q_ptr;
        }

        /// \brief function that is called for each children quadrant which is generated when calling 
        /// p4est_refine(). The purpose of this function is to set initialize the QuadData structure for 
        /// the newly created quadrants.
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of current tree object
        /// @param[in] quad pointer to current quadrant object

        void pXest_refine_init_fn ( p8est_t* forest, p4est_topidx_t tree_id, p8est_quadrant_t* quad )
        {
            // initialize quad data for new quadrant       
            QuadData* q_ptr = new QuadData ( );
            quad->p.user_data = q_ptr;
        }

        /// \brief function which is called on each parent and all its children quadrants in p4est_refine after 
        /// p4est_refine_init_fn() is called for the children quadrants.
        /// The purpose of this function is to copy the data from parent quadrant to its children quadrants
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] which_tree global id of current tree object
        /// @param[in] num_outgoing number of quadrants to be removed when refining a given quadrant
        /// @param[in] outgoing pointer to quadrant to be removed (i.e. refined)
        /// @param[in] num_incoming number of newly created quadrants
        /// @param[in] incoming pointers to newly created quadrants

        void pXest_refine_replace_fn ( p4est_t* forest, p4est_topidx_t which_tree,
                                       int num_outgoing, p4est_quadrant_t *outgoing[],
                                       int num_incoming, p4est_quadrant_t *incoming[] )
        {
            assert ( num_outgoing == 1 );
            assert ( num_incoming == 4 );

            // get pointer to external user data
            ForestRefineData* i_ptr = ( ForestRefineData * ) forest->user_pointer;
            MeshPXestDatabasePtr db_ptr = i_ptr->db;

            // get pointer to parent data
            QuadData* parent_ptr = pXest_get_quad_data_ptr ( outgoing[0] );
            assert ( parent_ptr != NULL );
            EntityNumber parent_index = parent_ptr->get_cell_index ( );

            // Loop over children quadrants
            for ( int jq = 0; jq < num_incoming; ++jq )
            {
                // get pointer to data
                p4est_quadrant_t* cur_quad = incoming[jq];
                QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );

                // copy from parent data
                q_ptr->set ( *parent_ptr );

                // get local quad number within family
                int p4est_sibling = pXest_get_sibling_nr ( cur_quad );
                int hiflow_sibling = pXest_z_order_to_ccw_order ( p4est_sibling );

                // set Id and index for child
                Id own_id = i_ptr->children_ids[parent_index][hiflow_sibling];
                q_ptr->set_cell_id ( own_id, cur_quad->level );
                q_ptr->set_remote_cell_id ( own_id, cur_quad->level );
                q_ptr->set_cell_index ( -1 );
                q_ptr->set_remote_cell_index ( -1 );

                // add new entity_to_quad mapping in mesh database
                QuadCoord coord = pXest_get_quad_coord ( cur_quad, which_tree );
                db_ptr->add_entity_to_quad_coord ( 2, own_id, coord, 1 );

                LOG_DEBUG ( 3, " parent index " << parent_index
                            << " parent x " << outgoing[0]->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " parent y " << outgoing[0]->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own x " << cur_quad->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own y " << cur_quad->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " p4est sibling id " << p4est_sibling << " hiflow sibling id " << hiflow_sibling
                            << " cell id " << own_id
                            << " owner rank: " << forest->mpirank );
            }
        }

        /// \brief function which is called on each parent and all its children quadrants in p4est_refine after 
        /// p4est_refine_init_fn() is called for the children quadrants.
        /// The purpose of this function is to copy the data from parent quadrant to its children quadrants
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] which_tree global id of current tree object
        /// @param[in] num_outgoing number of quadrants to be removed when refining a given quadrant
        /// @param[in] outgoing pointer to quadrant to be removed (i.e. refined)
        /// @param[in] num_incoming number of newly created quadrants
        /// @param[in] incoming pointers to newly created quadrants

        void pXest_refine_replace_fn ( p8est_t* forest, p4est_topidx_t which_tree,
                                       int num_outgoing, p8est_quadrant_t *outgoing[],
                                       int num_incoming, p8est_quadrant_t *incoming[] )
        {
            assert ( num_outgoing == 1 );
            assert ( num_incoming == 8 );

            // get pointer to external user data
            ForestRefineData* i_ptr = ( ForestRefineData * ) forest->user_pointer;
            MeshPXestDatabasePtr db_ptr = i_ptr->db;

            // get pointer to parent data
            QuadData* parent_ptr = pXest_get_quad_data_ptr ( outgoing[0] );
            assert ( parent_ptr != NULL );
            EntityNumber parent_index = parent_ptr->get_cell_index ( );

            // Loop over children quadrants
            for ( int jq = 0; jq < num_incoming; ++jq )
            {
                // get pointer to data
                p8est_quadrant_t* cur_quad = incoming[jq];
                QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );

                // copy from parent data
                q_ptr->set ( *parent_ptr );

                // get local quad number within family
                int p4est_sibling = pXest_get_sibling_nr ( cur_quad );
                int hiflow_sibling = pXest_z_order_to_ccw_order ( p4est_sibling );

                // set own Id
                Id own_id = i_ptr->children_ids[parent_index][hiflow_sibling];
                q_ptr->set_cell_id ( own_id, cur_quad->level );
                q_ptr->set_remote_cell_id ( own_id, cur_quad->level );
                q_ptr->set_cell_index ( -1 );
                q_ptr->set_remote_cell_index ( -1 );

                // add new entity_to_quad mapping in mesh database
                QuadCoord coord = pXest_get_quad_coord ( cur_quad, which_tree );
                db_ptr->add_entity_to_quad_coord ( 3, own_id, coord, 1 );

                LOG_DEBUG ( 3, " parent index " << parent_index
                            << " parent x " << outgoing[0]->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " parent y " << outgoing[0]->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own x " << cur_quad->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own y " << cur_quad->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own z " << cur_quad->z * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " p4est sibling id " << p4est_sibling << " hiflow sibling id " << hiflow_sibling
                            << " cell id " << own_id
                            << " owner rank: " << forest->mpirank );
            }
        }

        /// \brief Function which is called for each coarsenbale family of quadrants. This function decides whether the given 
        /// family should be coarsened
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of current tree object
        /// @param[in] quads pointers to quadrants forming a family
        /// @return true if family of quads shouzld be coarsened

        int pXest_coarsen_decide_fn ( p4est_t* forest, p4est_topidx_t tree, p4est_quadrant_t* quads[] )
        {
            // get pointer to external user data
            ForestCoarsenData* i_ptr = ( ForestCoarsenData * ) forest->user_pointer;

            // count number of set coarsen markers in family
            int num_markers = 0;
            for ( int jq = 0; jq < 4; ++jq )
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[jq] );
                assert ( q_ptr != NULL );
                num_markers += q_ptr->do_coarsen ( );
                q_ptr->set_coarsen ( 0 );
            }

            // Check if sufficiently many cells are marked for coarsening
            if ( num_markers <= -4 )
            {// yes
                LOG_DEBUG ( 2, "Tree " << tree << " num_coarsen_markers " << num_markers << " -> coarsen" );
                // loop over all quads in family
                for ( int jq = 0; jq < 4; ++jq )
                {
                    QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[jq] );

                    // save Ids of replaced cells in user data 
                    Id cell_id = q_ptr->get_cell_id ( quads[jq]->level );
                    i_ptr->coarsened_cells.find_insert ( cell_id );

                    // reset Id of finest level
                    q_ptr->set_cell_id ( -1, quads[jq]->level );
                    q_ptr->set_remote_cell_id ( -1, quads[jq]->level );

                    // add Id of parent cell to user data
                    i_ptr->parent_cells.find_insert ( q_ptr->get_cell_id ( quads[jq]->level - 1 ) );
                }
                return true;
            }
            else
            {// no
                LOG_DEBUG ( 2, "Tree " << tree << " num_coarsen_markers " << num_markers << " -> keep" );
                return false;
            }
        }

        /// \brief Function which is called for each coarsenbale family of quadrants. This function decides whether the given 
        /// family should be coarsened
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of current tree object
        /// @param[in] quads pointers to quadrants forming a family
        /// @return true if family of quads shouzld be coarsened

        int pXest_coarsen_decide_fn ( p8est_t* forest, p4est_topidx_t tree, p8est_quadrant_t* quads[] )
        {
            // get pointer to external user data
            ForestCoarsenData* i_ptr = ( ForestCoarsenData * ) forest->user_pointer;

            // count number of set coarsen markers in family
            int num_markers = 0;
            for ( int jq = 0; jq < 8; ++jq )
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[jq] );
                assert ( q_ptr != NULL );
                num_markers += q_ptr->do_coarsen ( );
                q_ptr->set_coarsen ( 0 );
            }

            // Check if sufficiently many cells are marked for coarsening
            if ( num_markers <= -8 )
            {// yes
                LOG_DEBUG ( 2, "Tree " << tree << " num_coarsen_markers " << num_markers << " -> coarsen" );
                // loop over all quads in family
                for ( int jq = 0; jq < 8; ++jq )
                {
                    QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[jq] );

                    // save Ids of replaced cells in user data 
                    Id cell_id = q_ptr->get_cell_id ( quads[jq]->level );
                    i_ptr->coarsened_cells.find_insert ( cell_id );

                    // reset Id of finest level
                    q_ptr->set_cell_id ( -1, quads[jq]->level );
                    q_ptr->set_remote_cell_id ( -1, quads[jq]->level );

                    // add Id of parent cell to user data
                    i_ptr->parent_cells.find_insert ( q_ptr->get_cell_id ( quads[jq]->level - 1 ) );
                }
                return true;
            }
            else
            {// no
                LOG_DEBUG ( 2, "Tree " << tree << " num_coarsen_markers " << num_markers << " -> keep" );
                return false;
            }
        }

        /// \brief Function which is called for each quadrant that is created when coarsening a faily of quadrants. 
        /// Here, the user data structure is initialized
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of current tree object
        /// @param[in] quad pointer to current quadrant object

        void pXest_coarsen_init_fn ( p4est_t* forest, p4est_topidx_t tree_id, p4est_quadrant_t* quad )
        {
            // init user data for new quadrant
            QuadData* q_ptr = new QuadData ( );
            quad->p.user_data = q_ptr;
        }

        /// \brief Function which is called for each quadrant that is created when coarsening a faily of quadrants. 
        /// Here, the user data structure is initialized
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of current tree object
        /// @param[in] quad pointer to current quadrant object

        void pXest_coarsen_init_fn ( p8est_t* forest, p4est_topidx_t tree_id, p8est_quadrant_t* quad )
        {
            // init user data for new quadrant
            QuadData* q_ptr = new QuadData ( );
            quad->p.user_data = q_ptr;
        }

        /// \brief function which is called on each parent and all its children quadrants in p4est_coarsen after 
        /// p4est_coarsen_init_fn() is called for the parent quadrants.
        /// The purpose of this function is to copy the data from chidlren quadrants to their parent quadrant
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] which_tree global id of current tree object
        /// @param[in] num_outgoing number of quadrants to be removed when coarsening a given family
        /// @param[in] outgoing pointers to quadrants to be removed (i.e. coarsened)
        /// @param[in] num_incoming number of newly created quadrant (parent of removed quadrants)
        /// @param[in] incoming pointer to newly created quadrant

        void pXest_coarsen_replace_fn ( p4est_t *forest, p4est_topidx_t which_tree,
                                        int num_outgoing, p4est_quadrant_t *outgoing[],
                                        int num_incoming, p4est_quadrant_t *incoming[] )
        {
            assert ( num_outgoing == 4 );
            assert ( num_incoming == 1 );

            // get data of some child 
            QuadData* child_ptr = pXest_get_quad_data_ptr ( outgoing[0] );

            // move data from child to parent
            QuadData* q_ptr = pXest_get_quad_data_ptr ( incoming[0] );
            q_ptr->set ( *child_ptr );
            q_ptr->set_cell_index ( -1 );
            q_ptr->set_remote_cell_index ( -1 );
        }

        /// \brief function which is called on each parent and all its children quadrants in p4est_coarsen after 
        /// p4est_coarsen_init_fn() is called for the parent quadrants.
        /// The purpose of this function is to copy the data from chidlren quadrants to their parent quadrant
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] which_tree global id of current tree object
        /// @param[in] num_outgoing number of quadrants to be removed when coarsening a given family
        /// @param[in] outgoing pointers to quadrants to be removed (i.e. coarsened)
        /// @param[in] num_incoming number of newly created quadrant (parent of removed quadrants)
        /// @param[in] incoming pointer to newly created quadrant

        void pXest_coarsen_replace_fn ( p8est_t *forest, p4est_topidx_t which_tree,
                                        int num_outgoing, p8est_quadrant_t *outgoing[],
                                        int num_incoming, p8est_quadrant_t *incoming[] )
        {
            assert ( num_outgoing == 8 );
            assert ( num_incoming == 1 );

            // get data of some child 
            QuadData* child_ptr = pXest_get_quad_data_ptr ( outgoing[0] );

            // move data from child to parent
            QuadData* q_ptr = pXest_get_quad_data_ptr ( incoming[0] );
            q_ptr->set ( *child_ptr );
            q_ptr->set_cell_index ( -1 );
            q_ptr->set_remote_cell_index ( -1 );
        }

        /// \brief Function which is called for each coarsenbale family of quadrants. This function decides whether the given 
        /// family should be coarsened. Here the return value is always false, since this function is used to track all coarsenable families of quadrants
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of current tree object
        /// @param[in] quads pointers to quadrants forming a family
        /// @return always false

        int pXest_coarsen_find_fn ( p4est_t* forest, p4est_topidx_t tree, p4est_quadrant_t* quads[] )
        {
            // get pointer to external user data
            ForestPatchData* i_ptr = ( ForestPatchData * ) forest->user_pointer;

            // Loop over coarsenable family and store ids of members
            for ( int jq = 0; jq < 4; ++jq )
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[jq] );
                assert ( q_ptr != NULL );

                EntityNumber index = q_ptr->get_cell_index ( );
                assert ( index >= 0 );

                i_ptr->coarsenable_cells.find_insert ( index );
            }
            return false;
        }

        /// \brief Function which is called for each coarsenbale family of quadrants. This function decides whether the given 
        /// family should be coarsened. Here the return value is always false, since this function is used to track all coarsenable families of quadrants
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree global id of current tree object
        /// @param[in] quads pointers to quadrants forming a family
        /// @return always false

        int pXest_coarsen_find_fn ( p8est_t* forest, p4est_topidx_t tree, p8est_quadrant_t* quads[] )
        {
            // get pointer to external user data
            ForestPatchData* i_ptr = ( ForestPatchData * ) forest->user_pointer;

            // Loop over coarsenable family and store ids of members
            for ( int jq = 0; jq < 8; ++jq )
            {
                QuadData* q_ptr = pXest_get_quad_data_ptr ( quads[jq] );
                assert ( q_ptr != NULL );

                EntityNumber index = q_ptr->get_cell_index ( );
                assert ( index >= 0 );

                i_ptr->coarsenable_cells.find_insert ( index );
            }
            return false;
        }

        /// \brief function that is called for each children quadrant which is generated when calling 
        /// p4est_balance(). The purpose of this function is to set initialize the QuadData structure for 
        /// the newly created quadrants.
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of current tree object
        /// @param[in] quad pointer to current quadrant object

        void pXest_balance_init_fn ( p4est_t* forest, p4est_topidx_t tree_id, p4est_quadrant_t* quad )
        {
            // initialize quad data for new quadrant       
            QuadData* q_ptr = new QuadData ( );
            quad->p.user_data = q_ptr;
        }

        /// \brief function that is called for each children quadrant which is generated when calling 
        /// p4est_balance(). The purpose of this function is to set initialize the QuadData structure for 
        /// the newly created quadrants.
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] tree_id global id of current tree object
        /// @param[in] quad pointer to current quadrant object

        void pXest_balance_init_fn ( p8est_t* forest, p4est_topidx_t tree_id, p8est_quadrant_t* quad )
        {
            // initialize quad data for new quadrant       
            QuadData* q_ptr = new QuadData ( );
            quad->p.user_data = q_ptr;
        }

        /// \brief function which is called on each parent and all its children quadrants in p4est_refine after 
        /// p4est_refine_init_fn() is called for the children quadrants.
        /// The purpose of this function is to copy the data from parent quadrant to its children quadrants
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] which_tree global id of current tree object
        /// @param[in] num_outgoing number of quadrants to be removed when refining a given quadrant
        /// @param[in] outgoing pointer to quadrant to be removed (i.e. refined)
        /// @param[in] num_incoming number of newly created quadrants
        /// @param[in] incoming pointers to newly created quadrants

        void pXest_balance_replace_fn ( p4est_t* forest, p4est_topidx_t which_tree,
                                        int num_outgoing, p4est_quadrant_t *outgoing[],
                                        int num_incoming, p4est_quadrant_t *incoming[] )
        {
            assert ( num_outgoing == 1 );
            assert ( num_incoming == 4 );

            // get pointer to external user data
            ForestBalanceData* i_ptr = ( ForestBalanceData * ) forest->user_pointer;

            // get data and index of parent 
            // Note: in case of multiple refinement, this data always corresponds to the initial quadrant to be refined
            // in the following called 'ancestor'
            QuadData* ancestor_ptr = pXest_get_quad_data_ptr ( outgoing[0] );
            assert ( ancestor_ptr != NULL );

            if ( DEBUG_LEVEL >= 2 )
            {
                ancestor_ptr->print ( );
            }

            // add ancestor_index as replaced cell to user data 
            EntityNumber ancestor_index = ancestor_ptr->get_cell_index ( );
            i_ptr->replaced_cell_indices.find_insert ( ancestor_index );

            // Loop over children quadrants
            for ( int jq = 0; jq < num_incoming; ++jq )
            {
                // get pointer to data
                p4est_quadrant_t* cur_quad = incoming[jq];
                QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                assert ( q_ptr != NULL );

                // copy data from parent to child  
                q_ptr->set ( *ancestor_ptr );

                // add newly created quadrants with ancestor given by ancestor_index to user data
                std::vector<p4est_quadrant_t> new_quad ( 1 );
                new_quad[0].level = cur_quad->level;
                new_quad[0].x = cur_quad->x;
                new_quad[0].y = cur_quad->y;
                new_quad[0].p.which_tree = q_ptr->get_tree_id ( );

                int number_desc_quads;
                std::map<EntityNumber, std::vector<p4est_quadrant_t> >::iterator it_q = i_ptr->new_quads4.find ( ancestor_index );
                if ( it_q == i_ptr->new_quads4.end ( ) )
                {
                    // First refinement of cell with ancestor_index
                    i_ptr->new_quads4.insert ( std::pair<EntityNumber, std::vector<p4est_quadrant_t> > ( ancestor_index, new_quad ) );
                    number_desc_quads = 1;
                }
                else
                {
                    // Repeated refinement of cell with ancestor_index
                    it_q->second.push_back ( new_quad[0] );
                    number_desc_quads = it_q->second.size ( );
                }

                LOG_DEBUG ( 2, " ancestor index " << ancestor_index
                            << " on process " << forest->mpirank
                            << " parent x " << outgoing[0]->x /** static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )*/
                            << " parent y " << outgoing[0]->y /** static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )*/
                            << " parent level " << ( int ) outgoing[0]->level
                            //<< " own x " << cur_quad->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            //<< " own y " << cur_quad->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " x " << cur_quad->x
                            << " y " << cur_quad->y
                            << " own level " << ( int ) cur_quad->level
                            << " own level " << ( int ) incoming[jq]->level
                            << " number of desc quads: " << number_desc_quads );
            }
        }

        /// \brief function which is called on each parent and all its children quadrants in p4est_refine after 
        /// p4est_refine_init_fn() is called for the children quadrants.
        /// The purpose of this function is to copy the data from parent quadrant to its children quadrants
        /// @param[in] forest pointer to underlying forest object
        /// @param[in] which_tree global id of current tree object
        /// @param[in] num_outgoing number of quadrants to be removed when refining a given quadrant
        /// @param[in] outgoing pointer to quadrant to be removed (i.e. refined)
        /// @param[in] num_incoming number of newly created quadrants
        /// @param[in] incoming pointers to newly created quadrants

        void pXest_balance_replace_fn ( p8est_t* forest, p4est_topidx_t which_tree,
                                        int num_outgoing, p8est_quadrant_t *outgoing[],
                                        int num_incoming, p8est_quadrant_t *incoming[] )
        {
            assert ( num_outgoing == 1 );
            assert ( num_incoming == 8 );

            // get pointer to external user data
            ForestBalanceData* i_ptr = ( ForestBalanceData * ) forest->user_pointer;

            // get data and index of parent 
            // Note: in case of multiple refinement, this data always corresponds to the initial quadrant to be refined
            // in the following called 'ancestor'
            QuadData* ancestor_ptr = pXest_get_quad_data_ptr ( outgoing[0] );
            assert ( ancestor_ptr != NULL );
            EntityNumber ancestor_index = ancestor_ptr->get_cell_index ( );

            // Loop over children quadrants
            for ( int jq = 0; jq < num_incoming; ++jq )
            {
                // get pointer to data
                p8est_quadrant_t* cur_quad = incoming[jq];
                QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                assert ( q_ptr != NULL );

                // copy from ancestor data
                q_ptr->set ( *ancestor_ptr );

                // add ancestor_index as replaced cell to user data 
                i_ptr->replaced_cell_indices.find_insert ( ancestor_index );

                // add newly created quadrants with ancestor given by ancestor_index to user data
                std::vector<p8est_quadrant_t> new_quad ( 1 );
                new_quad[0].level = cur_quad->level;
                new_quad[0].x = cur_quad->x;
                new_quad[0].y = cur_quad->y;
                new_quad[0].z = cur_quad->z;
                new_quad[0].p.which_tree = q_ptr->get_tree_id ( );

                int number_desc_quads;
                std::map<EntityNumber, std::vector<p8est_quadrant_t> >::iterator it_q = i_ptr->new_quads8.find ( ancestor_index );
                if ( it_q == i_ptr->new_quads8.end ( ) )
                {// First refinement of cell with ancestor_index
                    i_ptr->new_quads8.insert ( std::pair<EntityNumber, std::vector<p8est_quadrant_t> > ( ancestor_index, new_quad ) );
                    number_desc_quads = 1;
                }
                else
                {// Repeated refinement of cell with ancestor_index
                    it_q->second.push_back ( new_quad[0] );
                    number_desc_quads = it_q->second.size ( );
                }

                LOG_DEBUG ( 2, " ancestor index " << ancestor_index
                            << " on process " << forest->mpirank
                            << " parent x " << outgoing[0]->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " parent y " << outgoing[0]->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own x " << cur_quad->x * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own y " << cur_quad->y * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own z " << cur_quad->z * static_cast < int > ( std::pow ( static_cast < double > ( 2 ), -P4EST_MAXLEVEL ) )
                            << " own level: " << ( int ) cur_quad->level
                            << " number of desc quads: " << number_desc_quads );
            }
        }
#endif
        //////////////////////////////////////////////////////////////////////////////
        /////////// MeshPXest ////////////////////////////////////////////////////////

        /////////// public + mesh interface //////////////////////////////////////////

        MeshPXest::MeshPXest ( TDim tdim, GDim gdim,
                               MeshPXestDatabasePtr db,
                               const std::vector<Id>& cells,
                               int history_index,
                               const std::vector<bool>& cell_is_local,
                               std::vector<MasterSlave> period )
        :
        MeshDbView ( tdim, gdim, db, cells, period )
        {
            assert ( history_index >= 0 );
            assert ( cells.size ( ) == cell_is_local.size ( ) || cell_is_local.size ( ) == 0 );

            this->pXest_db_ = db;
            this->history_index_ = history_index;

            // If db contains forest: Set cell indices for local quadrants, without ghost cells
            this->update_cell_index_in_quads ( true );

            // Push back mesh history for cells
            std::vector<int> history_indices ( cells.size ( ), history_index );
            this->pXest_db_->add_cell_to_mesh_maps ( cells, history_indices );

            // if mesh has been refined at least once, setup sub cell numbers attribute
            if ( history_index >= 1 )
            {
                this->update_sub_cell_numbers ( );
            }

            // Set all cell locality flags
            if ( cell_is_local.size ( ) == 0 )
            {
                this->cell_is_local_.resize ( cells.size ( ), true );
            }
            else
            {
                this->cell_is_local_ = cell_is_local;
            }

            this->connection_mode_ = 0;
            this->patch_mode_ = false;
        }

        MeshPXest::MeshPXest ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        :
        MeshDbView ( tdim, gdim, period )
        {
            this->history_index_ = -1;
            this->connection_mode_ = 0;
            this->patch_mode_ = false;
        }

        MeshPtr MeshPXest::refine ( std::vector<EntityNumber>& refinements ) const
        {
#ifdef WITH_P4EST
            assert ( static_cast < int > ( refinements.size ( ) ) == this->num_entities ( tdim ( ) ) );

            treeId first_treeId = -1;
            treeId last_treeId = -2;

            const TDim tdim = this->tdim ( );
            const GDim gdim = this->gdim ( );
            const EntityCount num_cells = this->num_entities ( tdim );

            // Cells to be added to refined mesh
            SortedArray<Id> new_cells;
            std::map<Id, bool> new_cell_is_local;

            // Children ids per parent cell
            std::vector< std::vector<Id> > children_ids_per_parent ( num_cells );

            int rank;

            // **************************************************
            // Refine // coarsen mesh
            if ( tdim == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                rank = forest->mpirank;

                // ***************************************************
                // Collect data from refinement array and create new cell entities

                // Loop over all local trees, i.e. ignore refinement flags for ghost cells
                treeId first_local_tree = pXest_get_first_local_treeId ( forest );
                treeId last_local_tree = pXest_get_last_local_treeId ( forest );

                for ( treeId jt = first_local_tree; jt <= last_local_tree; ++jt )
                {
                    p4est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                    int num_quads = tree->quadrants.elem_count;

                    // Loop over all local quads in tree
                    for ( int jq = 0; jq < num_quads; ++jq )
                    {
                        p4est_quadrant* quad = pXest_get_local_quad_in_tree ( tree, jq );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
                        EntityNumber index = q_ptr->get_cell_index ( );

                        // check if quad corresponds to cell in current mesh
                        if ( index >= 0 )
                        {
                            // get refinement flag and entity
                            const int r = refinements[index];
                            const Entity current_cell = this->get_entity ( tdim, index );
                            const Id cell_id = current_cell.id ( );

                            if ( r > 0 )
                            {
                                // refinement (subdivision)
                                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Refining cell " << index << " with id " << cell_id << " with refinement " << 0 );

                                // refine cell entity
                                const RefinementTree* reftree = current_cell.cell_type ( ).refinement_tree ( 0 );
                                std::vector<int> local_sub_cell_numbers;
                                const std::vector<Id> children_cells = this->compute_refined_cells ( current_cell, reftree, local_sub_cell_numbers );
                                for ( int l = 0; l < children_cells.size ( ); ++l )
                                {
                                    new_cells.find_insert ( children_cells[l] );
                                    new_cell_is_local[children_cells[l]] = true;
                                }

                                // set refine flag in quad
                                q_ptr->set_refine ( true );
                                children_ids_per_parent[index].insert ( children_ids_per_parent[index].end ( ), children_cells.begin ( ), children_cells.end ( ) );
                            }
                            if ( r == 0 )
                            {
                                // keep cell
                                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Keep cell " << index << " with id " << cell_id );
                                new_cells.find_insert ( cell_id );
                                new_cell_is_local[cell_id] = true;
                                q_ptr->set_refine ( false );
                                children_ids_per_parent[index].push_back ( -1 );
                            }
                            if ( r < 0 )
                            {
                                // mark for coarsening
                                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Mark cell " << index << " with id " << cell_id << " for coarsening " << r );
                                new_cells.find_insert ( cell_id );
                                new_cell_is_local[cell_id] = true;
                                q_ptr->set_coarsen ( r );
                            }
                        }
                        else
                        {
                            q_ptr->set_refine ( false );
                            q_ptr->set_coarsen ( 0 );
                        }
                    }
                }

                // Loop over all ghost cells to ensure that current ghost cells are also present in refined mesh
                p4est_ghost_t* ghost = this->pXest_db_->get_p4est_ghost ( );
                if ( ghost != NULL )
                {
                    for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                    {
                        // get ghost quad
                        p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost->ghosts, jq );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                        assert ( q_ptr != NULL );

                        int level = ghost_quad->level;
                        const Id cell_id = q_ptr->get_cell_id ( level );

                        // keep ghost cell only if it is already contained in mesh
                        int pos;
                        if ( this->get_entities ( tdim ).find ( cell_id, &pos ) )
                        {
                            LOG_DEBUG ( 2, "[" << forest->mpirank << "] Keep ghost cell with " << cell_id );
                            new_cells.find_insert ( cell_id );
                            new_cell_is_local[cell_id] = false;
                            q_ptr->set_refine ( false );
                            q_ptr->set_coarsen ( 0 );
                        }
                    }
                }

                // ***************************************************
                // refine forest 
                ForestRefineData refine_data ( forest->mpirank, this->pXest_db_ );
                refine_data.set_children_ids ( children_ids_per_parent );
                pXest_set_forest_user_ptr ( forest, &refine_data );
                p4est_refine_ext ( forest, 0, -1, pXest_refine_decide_fn, pXest_refine_init_fn, pXest_refine_replace_fn );

                // resort quadrant arrays in forest
                pXest_sort_quad_arrays_in_forest ( forest );

                LOG_DEBUG ( 1, "[" << forest->mpirank << "] Number of cells after refinement: " << new_cells.size ( ) );
                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Cell ids after refinement: " << string_from_range ( new_cells.begin ( ), new_cells.end ( ) ) );

                // ***************************************************
                // coarsen forest 
                ForestCoarsenData coarsen_data;
                pXest_set_forest_user_ptr ( forest, &coarsen_data );
                p4est_coarsen_ext ( forest, 0, 0, pXest_coarsen_decide_fn, pXest_coarsen_init_fn, pXest_coarsen_replace_fn );

                // resort quadrant arrays in forest
                pXest_sort_quad_arrays_in_forest ( forest );

                // remove coarsened cells from list of marked cells
                for ( int l = 0; l < coarsen_data.coarsened_cells.size ( ); ++l )
                {
                    int pos;
                    if ( new_cells.find ( coarsen_data.coarsened_cells[l], &pos ) )
                    {
                        new_cells.erase ( pos );
                        new_cell_is_local.erase ( coarsen_data.coarsened_cells[l] );
                    }
                    // remove entity_to_quad maps
                    // this->pXest_db_->remove_entity_to_quad_maps(tdim, coarsen_data.coarsened_cells[l], 1);
                }

                // add parent cells of coarsened chidlren
                for ( int l = 0; l < coarsen_data.parent_cells.size ( ); ++l )
                {
                    new_cells.find_insert ( coarsen_data.parent_cells[l] );
                    new_cell_is_local[coarsen_data.parent_cells[l]] = true;
                }

                LOG_DEBUG ( 1, "[" << forest->mpirank << "] Number of cells after coarsening: " << new_cells.size ( ) );
                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Cell ids after coarsening: " << string_from_range ( new_cells.begin ( ), new_cells.end ( ) ) );
            }
            else if ( tdim == 3 )
            {
                // 3D_TODO new_cell_is_local ...

                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                rank = forest->mpirank;

                // ***************************************************
                // Collect data from refinement array and create new cell entities

                // Loop over all local trees, i.e. ignore refinement flags for ghost cells
                treeId first_local_tree = pXest_get_first_local_treeId ( forest );
                treeId last_local_tree = pXest_get_last_local_treeId ( forest );

                for ( treeId jt = first_local_tree; jt <= last_local_tree; ++jt )
                {
                    p8est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                    int num_quads = tree->quadrants.elem_count;

                    // Loop over all local quads in tree
                    for ( int jq = 0; jq < num_quads; ++jq )
                    {
                        p8est_quadrant* quad = pXest_get_local_quad_in_tree ( tree, jq );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
                        EntityNumber index = q_ptr->get_cell_index ( );

                        // check if quad corresponds to cell in current mesh
                        if ( index >= 0 )
                        {
                            // get refinement flag and entity
                            const int r = refinements[index];
                            const Entity current_cell = this->get_entity ( tdim, index );
                            const Id cell_id = current_cell.id ( );

                            if ( r > 0 )
                            {
                                // refinement (subdivision)
                                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Refining cell " << index << " with id " << cell_id << " with refinement " << 0 );

                                // refine cell entity
                                const RefinementTree* reftree = current_cell.cell_type ( ).refinement_tree ( 0 );
                                std::vector<int> local_sub_cell_numbers;
                                const std::vector<Id> children_cells = this->compute_refined_cells ( current_cell, reftree, local_sub_cell_numbers );
                                for ( int l = 0; l < children_cells.size ( ); ++l )
                                {
                                    new_cells.find_insert ( children_cells[l] );
                                    new_cell_is_local[children_cells[l]] = true;
                                }

                                // set refine flag in quad
                                q_ptr->set_refine ( true );
                                children_ids_per_parent[index].insert ( children_ids_per_parent[index].end ( ), children_cells.begin ( ), children_cells.end ( ) );
                            }
                            if ( r == 0 )
                            {
                                // keep cell
                                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Keep cell " << index << " with id " << cell_id );
                                new_cells.find_insert ( cell_id );
                                new_cell_is_local[cell_id] = true;
                                q_ptr->set_refine ( false );
                                children_ids_per_parent[index].push_back ( -1 );
                            }
                            if ( r < 0 )
                            {
                                // mark for coarsening
                                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Mark cell " << index << " with id " << cell_id << " for coarsening " << r );
                                new_cells.find_insert ( cell_id );
                                new_cell_is_local[cell_id] = true;
                                q_ptr->set_coarsen ( r );
                            }
                        }
                        else
                        {
                            q_ptr->set_refine ( false );
                            q_ptr->set_coarsen ( 0 );
                        }
                    }
                }

                // Loop over all ghost cells to ensure that current ghost cells are also present in refined mesh
                p8est_ghost_t* ghost = this->pXest_db_->get_p8est_ghost ( );
                if ( ghost != NULL )
                {
                    for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                    {
                        // get ghost quad
                        p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost->ghosts, jq );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                        int level = ghost_quad->level;
                        const Id cell_id = q_ptr->get_cell_id ( level );

                        // keep ghost cell only if it is already contained in mesh
                        int pos;
                        if ( this->get_entities ( tdim ).find ( cell_id, &pos ) )
                        {
                            LOG_DEBUG ( 2, "[" << forest->mpirank << "] Keep ghost cell with " << cell_id );
                            new_cells.find_insert ( cell_id );
                            new_cell_is_local[cell_id] = false;
                            q_ptr->set_refine ( false );
                            q_ptr->set_coarsen ( 0 );
                        }
                    }
                }

                // ***************************************************
                // refine forest 
                ForestRefineData refine_data ( forest->mpirank, this->pXest_db_ );
                refine_data.set_children_ids ( children_ids_per_parent );
                pXest_set_forest_user_ptr ( forest, &refine_data );
                p8est_refine_ext ( forest, 0, -1, pXest_refine_decide_fn, pXest_refine_init_fn, pXest_refine_replace_fn );

                // resort quadrant arrays in forest
                pXest_sort_quad_arrays_in_forest ( forest );

                LOG_DEBUG ( 1, "[" << forest->mpirank << "] Number of cells after refinement: " << new_cells.size ( ) );
                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Cell ids after refinement: " << string_from_range ( new_cells.begin ( ), new_cells.end ( ) ) );

                // ***************************************************
                // coarsen forest 
                ForestCoarsenData coarsen_data;
                pXest_set_forest_user_ptr ( forest, &coarsen_data );
                p8est_coarsen_ext ( forest, 0, 0, pXest_coarsen_decide_fn, pXest_coarsen_init_fn, pXest_coarsen_replace_fn );

                // resort quadrant arrays in forest
                pXest_sort_quad_arrays_in_forest ( forest );

                // remove coarsened cells from list of marked cells
                for ( int l = 0; l < coarsen_data.coarsened_cells.size ( ); ++l )
                {
                    int pos;
                    if ( new_cells.find ( coarsen_data.coarsened_cells[l], &pos ) )
                    {
                        new_cells.erase ( pos );
                        new_cell_is_local.erase ( coarsen_data.coarsened_cells[l] );
                    }
                    // remove entity_to_quad maps
                    // this->pXest_db_->remove_entity_to_quad_maps(tdim, coarsen_data.coarsened_cells[l], 1);
                }

                // add parent cells of coarsened chidlren
                for ( int l = 0; l < coarsen_data.parent_cells.size ( ); ++l )
                {
                    new_cells.find_insert ( coarsen_data.parent_cells[l] );
                    new_cell_is_local[coarsen_data.parent_cells[l]] = true;
                }

                LOG_DEBUG ( 1, "[" << forest->mpirank << "] Number of cells after coarsening: " << new_cells.size ( ) );
                LOG_DEBUG ( 2, "[" << forest->mpirank << "] Cell ids after coarsening: " << string_from_range ( new_cells.begin ( ), new_cells.end ( ) ) );
            }

            // ***************************************************
            // create intermediate mesh 
            // PERF_TODO fix possible memory leak occuring here
            std::vector<bool> locality;
            assert ( new_cells.size ( ) == new_cell_is_local.size ( ) );
            for ( int l = 0; l < new_cells.size ( ); ++l )
            {
                locality.push_back ( new_cell_is_local[new_cells.data ( ).at ( l )] );
            }

            MeshPtr refined_mesh = MeshPtr ( new MeshPXest ( this->tdim ( ),
                                                             this->gdim ( ),
                                                             this->pXest_db_,
                                                             new_cells.data ( ),
                                                             this->history_index_ + 1,
                                                             locality,
                                                             this->period_ ) );

            this->pXest_db_->set_mesh ( refined_mesh, this->history_index_ + 1 );
            LOG_DEBUG ( 1, " Built adapted intermediate mesh consisting of " << new_cells.size ( ) << " with history index " << this->history_index_ + 1 );

            boost::intrusive_ptr<MeshPXest> refined_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( refined_mesh );

            // ***************************************************
            // balance forest to ensure having at most one hanging node per facet -> add new cells add_cells
            std::vector<Id> add_cells;
            std::vector<Id> rm_cells;
            refined_mesh_pXest->balance ( add_cells, rm_cells, this->connection_mode_ );

            // remove balanced cells
            for ( int l = 0; l < rm_cells.size ( ); ++l )
            {
                int pos;
                if ( new_cells.find ( rm_cells[l], &pos ) )
                {
                    new_cells.erase ( pos );
                    new_cell_is_local.erase ( rm_cells[l] );
                }
            }

            // add new cells
            // CHECK_TODO: are locallly ownder clles balanced only? 
            for ( int l = 0; l < add_cells.size ( ); ++l )
            {
                new_cells.find_insert ( add_cells[l] );
                new_cell_is_local[add_cells[l]] = true;
            }
            LOG_DEBUG ( 1, "[" << rank << "] Number of cells after balancing: " << new_cells.size ( ) );
            LOG_DEBUG ( 2, "[" << rank << "] Cell ids after balancing: " << string_from_range ( new_cells.begin ( ), new_cells.end ( ) ) );

            // ***************************************************
            // create new mesh 
            // PERF_TODO fix possible memory leak occuring here
            locality.clear ( );
            assert ( new_cells.size ( ) == new_cell_is_local.size ( ) );
            for ( int l = 0; l < new_cells.size ( ); ++l )
            {
                locality.push_back ( new_cell_is_local[new_cells.data ( ).at ( l )] );
            }

            MeshPtr balanced_mesh = MeshPtr ( new MeshPXest ( this->tdim ( ),
                                                              this->gdim ( ),
                                                              this->pXest_db_,
                                                              new_cells.data ( ),
                                                              this->history_index_ + 2,
                                                              locality,
                                                              this->period_ ) );
            LOG_DEBUG ( 1, " Built adapted balanced mesh consisting of " << new_cells.size ( ) << " with history index " << this->history_index_ + 2 );

            // Put pointer to new mesh into database                                                  
            this->pXest_db_->set_mesh ( balanced_mesh, this->history_index_ + 2 );

            // ****************************************************
            // ensure that all cells are coarsenable
            if ( this->patch_mode_ )
            {
                boost::intrusive_ptr<MeshPXest> balanced_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( balanced_mesh );
                MeshPtr patched_mesh = balanced_mesh_pXest->make_mesh_coarsenable ( );

                boost::intrusive_ptr<MeshPXest> patched_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( patched_mesh );
                int patched_history_index = patched_mesh_pXest->history_index_;

                LOG_DEBUG ( 1, "[" << rank << "] Number of cells after patching: " << patched_mesh->num_entities ( tdim ) );
                LOG_DEBUG ( 1, " Built pacthed mesh with history index " << patched_history_index );

                return patched_mesh;
            }
            else
            {
                return balanced_mesh;
            }

            // PERF_TODO delete ForestXData?
#endif
        }

        Id MeshPXest::get_parent_cell_id ( EntityNumber cell_index ) const
        {
            assert ( cell_index >= 0 );
            const int tdim = this->tdim ( );
            Id cell_id = this->get_id ( tdim, cell_index );
            return this->get_ancestor_cell_id ( cell_id, 1 );
        }

        Entity MeshPXest::get_parent_cell ( EntityNumber cell_index ) const
        {
            assert ( cell_index >= 0 );
            const Id parent_id = this->get_parent_cell_id ( cell_index );
            const int parent_history_index = this->pXest_db_->get_last_mesh_history_index ( parent_id );

            if ( parent_history_index >= 0 )
            {
                Id cell_id = this->get_id ( this->tdim ( ), cell_index );
                const boost::intrusive_ptr<MeshPXest> parent_mesh = boost::static_pointer_cast<MeshPXest> ( this->pXest_db_->get_mesh ( parent_history_index ) );
                const EntityNumber parent_cell_index = parent_mesh->get_cell_index ( parent_id );

                return parent_mesh->get_entity ( this->tdim ( ), parent_cell_index );
            }
            else
            {
                assert ( false );
            }
            return Entity ( );
        }

        bool MeshPXest::cell_has_parent ( EntityNumber cell_index ) const
        {
            assert ( cell_index >= 0 );
            const Id parent_id = this->get_parent_cell_id ( cell_index );
            if ( parent_id >= 0 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        std::vector<Id> MeshPXest::get_children_cell_ids ( EntityNumber cell_index ) const
        {
            assert ( cell_index >= 0 );
            Id cell_id = this->get_id ( this->tdim ( ), cell_index );
            return this->get_descendant_cell_ids ( cell_id, 1 );
        }

        std::vector<EntityNumber> MeshPXest::get_sibling_cell_indices ( EntityNumber cell_index ) const
        {
            assert ( cell_index >= 0 );
            assert ( cell_index <= this->num_entities ( this->tdim ( ) ) );

            Id parent_id = this->get_parent_cell_id ( cell_index );
            if ( parent_id < 0 )
            {
                return std::vector<EntityNumber>( 0 );
            }

            std::vector<Id> sibling_ids = this->get_descendant_cell_ids ( parent_id, 1 );
            std::vector<EntityNumber> sibling_indices;
            for ( int l = 0; l < sibling_ids.size ( ); ++l )
            {
                assert ( sibling_ids[l] >= 0 );

                int sibling_index = this->get_cell_index ( sibling_ids[l] );
                if ( sibling_index >= 0 )
                {
                    sibling_indices.push_back ( sibling_index );
                }
            }
            return sibling_indices;
        }

        int MeshPXest::num_ghost_cells ( ) const
        {
#ifdef WITH_P4EST
            if ( this->tdim ( ) == 2 )
            {
                if ( this->pXest_db_->get_p4est_ghost ( ) != NULL )
                {
                    return this->pXest_db_->get_p4est_ghost ( )->ghosts.elem_count;
                }
            }
            else if ( this->tdim ( ) == 3 )
            {
                if ( this->pXest_db_->get_p8est_ghost ( ) != NULL )
                {
                    return this->pXest_db_->get_p8est_ghost ( )->ghosts.elem_count;
                }
            }
            return 0;
#else
            return 0;
#endif
        }

        bool MeshPXest::cell_is_local ( EntityNumber index ) const
        {
            assert ( index >= 0 );
            assert ( index < this->num_entities ( this->tdim ( ) ) );

            return this->cell_is_local_[index];
        }

        void MeshPXest::save ( std::string filename, const MPI_Comm& comm ) const
        {
#ifdef WITH_HDF5
            MeshDbView::save ( filename, comm );

            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );

            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "MeshPXest";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );

            write_map_parallel ( group_ptr, "cell_id_to_index", cell_id_to_index_, comm );

            write_array_parallel ( group_ptr, "history_index", &history_index_, 1, comm );

            std::vector<int> bool_to_int ( cell_is_local_.size ( ) );
            for ( int i = 0; i < bool_to_int.size ( ); ++i )
                bool_to_int[i] = ( int ) cell_is_local_[i];
            write_array_parallel<int>( group_ptr, "cell_is_local", bool_to_int, comm );

            write_array_parallel ( group_ptr, "connection_mode", &connection_mode_, 1, comm );

            int temp_patch_mode = (int) patch_mode_;
            write_array_parallel ( group_ptr, "patch_mode", &temp_patch_mode, 1, comm );
#endif
        }

        void MeshPXest::load ( std::string filename, const MPI_Comm& comm )
        {
#ifdef WITH_HDF5
            MeshDbView::load ( filename, comm );

            H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );

            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "MeshPXest";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );

            read_map_parallel ( group_ptr, "cell_id_to_index", cell_id_to_index_, comm );

            int* buffer;
            int dummy;
            read_array_parallel ( group_ptr, "history_index", buffer, dummy, comm );
            assert(dummy == 1);
            history_index_ = *buffer;
            delete buffer;

            std::vector<int> bool_to_int;

            read_array_parallel<int>( group_ptr, "cell_is_local", bool_to_int, comm );
            cell_is_local_.resize ( bool_to_int.size ( ) );
            for ( int i = 0; i < bool_to_int.size ( ); ++i )
                cell_is_local_[i] = ( bool ) bool_to_int[i];
            read_array_parallel ( group_ptr, "connection_mode", buffer, dummy, comm );
            assert(dummy == 1);
            connection_mode_ = *buffer;
            delete buffer;

            read_array_parallel ( group_ptr, "patch_mode", buffer, dummy, comm );
            assert(dummy == 1);
            patch_mode_ = ( bool ) *buffer;
            delete buffer;
#endif
        }

        void MeshPXest::copy_from ( const MeshPtr mesh )
        {
            // copy member variables defined in MeshDbView class
            MeshDbView::copy_from ( mesh );
            boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( mesh );

            const int tdim = mesh->tdim ( );

            for ( int dim = 0; dim < tdim + 1; ++dim )
            {
                assert ( this->num_entities ( dim ) == mesh->num_entities ( dim ) );
            }

            // copy some flags
            this->connection_mode_ = mesh_pXest->connection_mode_;
            this->patch_mode_ = mesh_pXest->patch_mode_;
            this->cell_is_local_ = mesh_pXest->cell_is_local_;
            this->history_index_ = mesh_pXest->history_index_;
            this->cell_id_to_index_ = mesh_pXest->cell_id_to_index_;

            // copy attribute
            Attribute* attr = mesh_pXest->sub_cell_number_attr_.get ( );
            if ( attr != NULL )
            {

                IntAttribute* int_attr = dynamic_cast < IntAttribute* > ( attr );

                int i = 0;
                std::vector<int> int_att ( int_attr->size ( ) );

                for ( std::vector<int>::iterator it_att = int_att.begin ( ); it_att < int_att.end ( ); ++it_att, ++i )
                {
                    *it_att = int_attr->get_int_value ( i );
                }

                AttributePtr new_attribute = AttributePtr ( new IntAttribute ( int_att ) );
                this->add_attribute ( std::string ( "__sub_cell_number__" ), tdim, new_attribute );

                this->sub_cell_number_attr_ = new_attribute;
            }

            // copy Database pointer
            this->pXest_db_ = mesh_pXest->pXest_db_;
        }

        void MeshPXest::deep_copy_from ( const MeshPtr mesh )
        {
            // standard copy
            this->copy_from ( mesh );

            // copy database
            MeshPXestDatabasePtr tmp_db ( new MeshPXestDatabase ( this->tdim ( ), this->gdim ( ) ) );
            boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( mesh );

            tmp_db->deep_copy_from ( mesh_pXest->get_db ( ) );
            this->set_db ( tmp_db );

            // replace pointer for self mesh object
            MeshPtr self_ptr ( this );
            this->get_db ( )->set_mesh ( self_ptr, this->history_index_ );
        }

        bool MeshPXest::is_uniformly_coarsenable ( ) const
        {
            const int tdim = this->tdim ( );
            int num_cells = this->num_entities ( tdim );

            // Find all cells that form a family which is coarsenable
            SortedArray<EntityNumber> coarsenable_cells;
            this->find_coarsenable_cells ( coarsenable_cells );

            // Find cells that are not coarsenable
            std::vector<EntityNumber> non_coarsenable_cells;
            for ( EntityNumber j = 0; j < num_cells; ++j )
            {
                int pos;
                if ( !coarsenable_cells.find ( j, &pos ) )
                {
                    if ( this->cell_is_local ( j ) )
                    {
                        non_coarsenable_cells.push_back ( j );
                    }
                }
            }

            LOG_DEBUG ( 1, "Number of non coarsenable cells: " << non_coarsenable_cells.size ( ) );
            if ( non_coarsenable_cells.size ( ) == 0 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /////////// public but not in mesh interface ////////////////////////////////

        void MeshPXest::set_connection_mode ( int mode )
        {
            this->connection_mode_ = mode;
            if ( this->patch_mode_ )
            {
                this->connection_mode_ = 2;
                if ( mode == 0 )
                {
                    LOG_DEBUG ( 0, "Caution: You set patch_mode = true. This choice requires connection_mode = 1 (2D) or 2 (3D)! " );
                }
            }
        }

        void MeshPXest::set_patch_mode ( bool flag )
        {
            this->patch_mode_ = flag;
            if ( flag )
            {
                this->connection_mode_ = 2;
            }
        }

        /////////// protected ///////////////////////////////////////////////////////

        void MeshPXest::create_interface_list ( InterfaceList* list ) const
        {
#ifdef WITH_P4EST
            LOG_DEBUG ( 1, "Create interface list by using p4est iterator" );
            list->clear ( );

            // iterate over forest and call function "p4est_interface_list_fn" on each facet
            // TODO facets zwischen ghost cells nötig?
            if ( this->tdim ( ) == 2 )
            {
                p4est_iterate ( this->pXest_db_->get_p4est_forest ( ),
                                this->pXest_db_->get_p4est_ghost ( ),
                                list,
                                NULL,
                                pXest_interface_list_fn,
                                NULL );
            }
            if ( this->tdim ( ) == 3 )
            {
                p8est_iterate ( this->pXest_db_->get_p8est_forest ( ),
                                this->pXest_db_->get_p8est_ghost ( ),
                                list,
                                NULL,
                                pXest_interface_list_fn,
                                NULL,
                                NULL );
            }
#endif
        }

        // ---------------------------------------------------------------------
        // BUG_TODO: check whether works fine in case of recursive refinement of balanced cell

        void MeshPXest::balance ( std::vector<Id>& add_cells, std::vector<Id>& rm_cells, int connection_mode ) const
        {
            LOG_DEBUG ( 1, "Balance Mesh" );
            add_cells.clear ( );
            rm_cells.clear ( );
#ifdef WITH_P4EST
            int tdim = this->tdim ( );

            // call p4est_balance routine.  Refined cell indices are stored in balance_data -> ancestor_cells 
            ForestBalanceData balance_data;
            p4est_t* forest4;
            p8est_t* forest8;

            if ( tdim == 2 )
            {
                forest4 = this->pXest_db_->get_p4est_forest ( );
                pXest_set_forest_user_ptr ( forest4, &balance_data );
                if ( connection_mode >= 1 )
                {
                    p4est_balance_ext ( forest4, P4EST_CONNECT_FULL, pXest_balance_init_fn, pXest_balance_replace_fn );
                }
                else
                {
                    p4est_balance_ext ( forest4, P4EST_CONNECT_FACE, pXest_balance_init_fn, pXest_balance_replace_fn );
                }

                // resort quadrant arrays in forest
                pXest_sort_quad_arrays_in_forest ( forest4 );

                assert ( pXest_quad_arrays_in_forest_sorted ( forest4 ) );
            }
            else if ( tdim == 3 )
            {
                forest8 = this->pXest_db_->get_p8est_forest ( );
                pXest_set_forest_user_ptr ( forest8, &balance_data );
                if ( connection_mode == 2 )
                {
                    p8est_balance_ext ( forest8, P8EST_CONNECT_FULL, pXest_balance_init_fn, pXest_balance_replace_fn );
                }
                else if ( connection_mode == 1 )
                {
                    p8est_balance_ext ( forest8, P8EST_CONNECT_EDGE, pXest_balance_init_fn, pXest_balance_replace_fn );
                }
                else
                {
                    p8est_balance_ext ( forest8, P8EST_CONNECT_FACE, pXest_balance_init_fn, pXest_balance_replace_fn );
                }

                // resort quadrant arrays in forest
                pXest_sort_quad_arrays_in_forest ( forest8 );

                assert ( pXest_quad_arrays_in_forest_sorted ( forest8 ) );
            }
            LOG_DEBUG ( 1, " P4est balancing done " );

            // get cells that were replaced
            int num_replaced_cells = balance_data.replaced_cell_indices.size ( );

            LOG_DEBUG ( 1, " Number of replaced cells " << balance_data.replaced_cell_indices.size ( ) );

            // Loop over all replaced ancestor cells and refine
            for ( int l = 0; l < num_replaced_cells; ++l )
            {
                // get ancestor cell, index and id *************************
                EntityNumber ancestor_index = balance_data.replaced_cell_indices[l];
                Id ancestor_id = this->get_id ( tdim, ancestor_index );
                const Entity ancestor_cell = this->get_entity ( tdim, ancestor_index );

                LOG_DEBUG ( 2, " Balance (i.e. refine) ancestor index " << ancestor_index << " with cell id " << ancestor_id );
                rm_cells.push_back ( ancestor_id );

                // get quadrant of replaced ancestor cell ******************
                QuadCoord ancestor_coord = this->pXest_db_->get_entity_to_quad_coord ( tdim, ancestor_id, 1 );

                if ( DEBUG_LEVEL >= 2 )
                {
                    ancestor_coord.print ( );
                }
                assert ( ancestor_coord.tree >= 0 );

                // setup refinementree and match its structure with newly created quadrants
                RefinementTree ref_tree ( tdim, ancestor_cell.cell_type ( ) );

                p4est_quadrant_t ancestor_quad4;
                p8est_quadrant_t ancestor_quad8;
                std::vector<p4est_quadrant_t> tree_quads4 ( 1 );
                std::vector<p8est_quadrant_t> tree_quads8 ( 1 );
                int num_desc;

                if ( tdim == 2 )
                {
                    ancestor_quad4.level = ancestor_coord.level;
                    ancestor_quad4.x = ancestor_coord.x;
                    ancestor_quad4.y = ancestor_coord.y;
                    ancestor_quad4.p.which_tree = ancestor_coord.tree;

                    int ancestor_level = ancestor_coord.level;
                    LOG_DEBUG ( 2, " Ancestor is on level " << ( int ) ancestor_quad4.level << " and on tree " << ancestor_quad4.p.which_tree );

                    // get descendant quadrants of current ancestor ************
                    num_desc = balance_data.new_quads4[ancestor_index].size ( );

                    tree_quads4[0] = ancestor_quad4;
                    tree_quads4.insert ( tree_quads4.end ( ), balance_data.new_quads4[ancestor_index].begin ( ), balance_data.new_quads4[ancestor_index].end ( ) );

                    assert ( tree_quads4.size ( ) == num_desc + 1 );

                    // create refined cells ************************************
                    // Note: not all of these cells are put into the mesh. Only those
                    // corresponding to a leaf quadrant stored in ForestBalanceData
                    pXest_build_ref_tree ( tree_quads4, &ref_tree );
                }
                else if ( tdim == 3 )
                {
                    ancestor_quad8.level = ancestor_coord.level;
                    ancestor_quad8.x = ancestor_coord.x;
                    ancestor_quad8.y = ancestor_coord.y;
                    ancestor_quad8.z = ancestor_coord.z;
                    ancestor_quad8.p.which_tree = ancestor_coord.tree;

                    int ancestor_level = ancestor_coord.level;
                    LOG_DEBUG ( 2, " Ancestor is on level " << ( int ) ancestor_quad8.level << " and on tree " << ancestor_quad8.p.which_tree );

                    // get descendant quadrants of current ancestor ************
                    num_desc = balance_data.new_quads8[ancestor_index].size ( );
                    tree_quads8[0] = ancestor_quad8;
                    tree_quads8.insert ( tree_quads8.end ( ), balance_data.new_quads8[ancestor_index].begin ( ), balance_data.new_quads8[ancestor_index].end ( ) );

                    assert ( tree_quads8.size ( ) == num_desc + 1 );

                    // create refined cells ************************************
                    // Note: not all of these cells are put into the mesh. Only those
                    // corresponding to a leaf quadrant stored in ForestBalanceData
                    pXest_build_ref_tree ( tree_quads8, &ref_tree );
                }

                int num_cells = ref_tree.num_cells ( );
                assert ( num_cells == num_desc + 1 );

                LOG_DEBUG ( 2, "RefTree has " << ref_tree.num_cells ( ) << " cells" );
                if ( DEBUG_LEVEL >= 2 )
                {
                    for ( int jn = 0; jn < num_cells; ++jn )
                    {
                        const RefinedCell& cur_cell = ref_tree.cell ( jn );

                        LOG_DEBUG ( 2, "Tree node " << jn << " has parent cell " << cur_cell.parent ( ) << " and "
                                    << cur_cell.num_children ( ) << " children and sub_cell_number " << cur_cell.sub_cell_number ( ) );

                        if ( cur_cell.parent ( ) >= 0 )
                        {
                            const RefinedCell& parent_cell = ref_tree.cell ( cur_cell.parent ( ) );
                            if ( tdim == 2 )
                            {
                                p4est_quadrant_t cur_quad = tree_quads4[cur_cell.sub_cell_number ( )];
                                p4est_quadrant_t parent_quad = tree_quads4[parent_cell.sub_cell_number ( )];
                                LOG_DEBUG ( 2, " This should be true: " << p4est_quadrant_is_parent ( &parent_quad, &cur_quad ) );
                            }
                            else if ( tdim == 3 )
                            {
                                p8est_quadrant_t cur_quad = tree_quads8[cur_cell.sub_cell_number ( )];
                                p8est_quadrant_t parent_quad = tree_quads8[parent_cell.sub_cell_number ( )];
                                LOG_DEBUG ( 2, " This should be true: " << p8est_quadrant_is_parent ( &parent_quad, &cur_quad ) );
                            }
                        }
                    }
                }

                // refine recursively
                std::vector<Id> desc_ids;
                std::vector<int> local_sub_cell_numbers;
                std::vector<int> tree_node_numbers;
                desc_ids = this->compute_refined_cells_recursive ( ancestor_cell, &ref_tree, local_sub_cell_numbers, tree_node_numbers, NULL, false );
                assert ( local_sub_cell_numbers.size ( ) == num_desc );

                LOG_DEBUG ( 2, "Refining cell " << ancestor_index << " yields " << desc_ids.size ( ) << " new cells with ids " );
                if ( DEBUG_LEVEL >= 2 )
                {
                    for ( int c = 0; c < desc_ids.size ( ); ++c )
                    {
                        LOG_DEBUG ( 2, desc_ids[c] << " " );
                    }
                }

                // put the data into the leaf quadrants ****************************
                // and select those newly created cells for the new mesh that correspond to a leaf quad

                // Loop over all cells in refinementTree except of root
                for ( int j = 0; j < tree_node_numbers.size ( ); ++j )
                {
                    // index of current cell in refined_cells array of tree
                    int jn = tree_node_numbers[j];
                    assert ( jn > 0 );

                    // get current cell
                    const RefinedCell& cur_cell = ref_tree.cell ( jn );
                    int jq = cur_cell.quadrant_number ( );
                    Id cur_id = desc_ids[j];

                    QuadCoord coord;
                    int cur_level;

                    p4est_quadrant_t cur_quad4;
                    p8est_quadrant_t cur_quad8;

                    if ( tdim == 2 )
                    {
                        assert ( jq >= 0 && jq < tree_quads4.size ( ) );
                        cur_quad4 = tree_quads4[jq];
                        cur_level = cur_quad4.level;
                        coord = pXest_get_quad_coord ( &cur_quad4, ancestor_coord.tree );
                        this->pXest_db_->add_entity_to_quad_coord ( 2, cur_id, coord, 1 );
                    }
                    else if ( tdim == 3 )
                    {
                        assert ( jq >= 0 && jq < tree_quads8.size ( ) );
                        cur_quad8 = tree_quads8[jq];
                        cur_level = cur_quad8.level;
                        coord = pXest_get_quad_coord ( &cur_quad8, ancestor_coord.tree );
                        this->pXest_db_->add_entity_to_quad_coord ( 3, cur_id, coord, 1 );
                    }

                    LOG_DEBUG ( 2, "Current quad with quadrant nr " << jq << " has level " << cur_level << ", tree_node_nr: " << jn << ", cell id: " << cur_id
                                << ", is leaf: " << !cur_cell.is_refined ( ) );

                    if ( !cur_cell.is_refined ( ) )
                    { // cur_quad is leaf

                        // put corresponding cell_id into add_cells 
                        add_cells.push_back ( cur_id );

                        // get QuadData pointer of current desc quadrant 
                        QuadData* q_ptr;

                        if ( tdim == 2 )
                        {
                            treeId tree_id = ancestor_quad4.p.which_tree;
                            p4est_tree_t* tree = pXest_get_tree_in_forest ( forest4, tree_id );
                            LOG_DEBUG ( 2, "Current leaf: Tree: " << tree_id << " level " << ( int ) cur_quad4.level << " x " << cur_quad4.x << " y " << cur_quad4.y );

                            q_ptr = pXest_get_quad_data_ptr ( tree, &cur_quad4 );
                        }
                        else if ( tdim == 3 )
                        {
                            treeId tree_id = ancestor_quad8.p.which_tree;
                            p8est_tree_t* tree = pXest_get_tree_in_forest ( forest8, tree_id );
                            q_ptr = pXest_get_quad_data_ptr ( tree, &cur_quad8 );
                        }
                        assert ( q_ptr != NULL );

                        // get associated cell ids of all quadrants between ancestor_quad and cur_quad 
                        // and store them in QuadData of leaf quadrant
                        RefinedCell const * it_cell = &cur_cell;
                        int it_level = cur_level;
                        int it_nr = jn;
                        int it_index = j;
                        bool reached_root = false;

                        // traverse tree from cur_cell to root
                        while ( !reached_root )
                        {
                            // get cell_id 
                            Id it_id;
                            if ( it_index == -1 )
                            {
                                it_id = ancestor_id;
                            }
                            else
                            {
                                it_id = desc_ids[it_index];
                            }

                            q_ptr->set_cell_id ( it_id, it_level );
                            q_ptr->set_remote_cell_id ( it_id, it_level );
                            LOG_DEBUG ( 2, " Traverse: it_cell_id: " << it_id << " level " << it_level );

                            // check termination criterion
                            if ( it_cell->is_root ( ) )
                            {
                                reached_root = true;
                            }
                            else
                            {
                                // get parent cell
                                it_nr = it_cell->parent ( );

                                // TODO effizienter
                                if ( it_nr == 0 )
                                {
                                    it_index = -1;
                                }
                                else
                                {
                                    for ( int l = 0; l < tree_node_numbers.size ( ); ++l )
                                    {
                                        if ( tree_node_numbers[l] == it_nr )
                                        {
                                            it_index = l;
                                            break;
                                        }
                                    }
                                }
                                LOG_DEBUG ( 2, " Traverse: parent nr: " << it_nr << " index " << it_index );
                                it_cell = &( ref_tree.cell ( it_nr ) );
                                --it_level;
                            }
                        }
                        if ( DEBUG_LEVEL >= 3 )
                        {
                            q_ptr->print ( );
                        }
                    }
                }
                // PERF_TODO delete ForestBalanceData struct?
            }
#endif
        }

        // BUG_TODO fails for num_ref_steps > 1

        MeshPtr MeshPXest::refine_uniform_seq ( int num_ref_steps ) const
        {
            assert ( num_ref_steps > 0 );
            MeshPtr refined_master_mesh = this->refine_uniform_seq ( );

            for ( int l = 1; l < num_ref_steps; ++l )
            {
                refined_master_mesh = refined_master_mesh->refine_uniform_seq ( );
            }
            return refined_master_mesh;
        }

        MeshPtr MeshPXest::refine_uniform_seq ( ) const
        {
            int num_ref_steps = 1;
            assert ( num_ref_steps == 1 );

            LOG_DEBUG ( 1, "Refine mesh uniformly " << num_ref_steps << " times" );

            SortedArray<Id> new_cells;
            const int tdim = this->tdim ( );
            const int gdim = this->gdim ( );
            int num_cells = this->num_entities ( tdim );

            this->get_db ( )->clear_entity_to_quad_map ( tdim, 0 );

            // Loop over all cells in current mesh
            for ( int index = 0; index < num_cells; ++index )
            {
                // refine cell entity recursively
                const Entity current_cell = this->get_entity ( tdim, index );
                RefinementTree ref_tree ( tdim, current_cell.cell_type ( ), num_ref_steps, static_cast < int > ( std::pow ( static_cast < double > ( 2 ), tdim ) ) );

                std::vector<int> local_sub_cell_numbers;
                std::vector<int> tree_node_numbers;
                std::vector<Id> desc_ids = this->compute_refined_cells_recursive
                        ( current_cell, &ref_tree, local_sub_cell_numbers, tree_node_numbers, NULL, true );

                // select leaf cells for new mesh
                // Loop over all cells in refinementTree except of root
                for ( int j = 0; j < tree_node_numbers.size ( ); ++j )
                {
                    // index of current cell in refined_cells array of tree
                    int jn = tree_node_numbers[j];
                    assert ( jn > 0 );

                    // get current cell
                    const RefinedCell& cur_cell = ref_tree.cell ( jn );

                    if ( !cur_cell.is_refined ( ) )
                    { // cell is leaf
                        new_cells.find_insert ( desc_ids[j] );
                    }
                    else
                    {
                        assert ( false );
                    }
                }
            }

            // create new mesh
            std::vector<bool> locality;
            MeshPtr refined_mesh = MeshPtr ( new MeshPXest ( tdim, gdim, this->pXest_db_, new_cells.data ( ), 0, locality, this->period_ ) );
            boost::intrusive_ptr<MeshPXest> refined_mesh_p4est = boost::static_pointer_cast<MeshPXest> ( refined_mesh );

            LOG_DEBUG ( 1, "Refined mesh consists of " << new_cells.size ( ) << " cells" );

            // adapt connectivity according to refined mesh
            if ( tdim == 2 )
            {
                assert ( num_cells * static_cast < int > ( std::pow ( static_cast < double > ( 4 ), num_ref_steps ) ) == refined_mesh_p4est->num_entities ( tdim ) );
            }
            if ( tdim == 3 )
            {
                assert ( num_cells * static_cast < int > ( std::pow ( static_cast < double > ( 8 ), num_ref_steps ) ) == refined_mesh_p4est->num_entities ( tdim ) );
            }

            num_cells = refined_mesh_p4est->num_entities ( tdim );
            int num_vertices = refined_mesh_p4est->num_entities ( 0 );
            assert ( num_vertices == this->pXest_db_->num_vertices ( ) );

            int num_vert_per_cell = static_cast < int > ( std::pow ( static_cast < double > ( 2 ), tdim ) );
            SortedArray<Id> refined_vertices;

            std::vector<double> p4est_vertices ( num_vertices * 3, 0. );
            std::vector<int> p4est_tree_to_vertices;

            // loop over all vertices in refined_mesh -> get vertex coordinates
            // loop over all cells in refined_mesh
            for ( EntityIterator it = refined_mesh_p4est->begin ( tdim ), end_it = refined_mesh_p4est->end ( tdim ); it != end_it; ++it )
            {
                // get vertex ids of current cell
                for ( VertexIdIterator it_v = it->begin_vertex_ids ( ); it_v != it->end_vertex_ids ( ); ++it_v )
                {
                    refined_vertices.find_insert ( *it_v );
                }
            }
            assert ( num_vertices = refined_vertices.size ( ) );

            for ( int jv = 0; jv < num_vertices; ++jv )
            {
                //Id v_id = refined_vertices.data().at(jv);
                //std::vector<Coordinate> coords = refined_mesh_p4est->get_db()->get_coordinates(v_id);
                std::vector<Coordinate> coords = refined_mesh_p4est->get_db ( )->get_coordinates ( jv );
                // 3rd coordinate stays 0 in case of 2d
                for ( int c = 0; c < gdim; ++c )
                {
                    p4est_vertices[3 * jv + c] = coords[c];
                }
            }

            // loop over all cells in refined_mesh
            int counter = 0;
            for ( EntityIterator it = refined_mesh_p4est->begin ( tdim ), end_it = refined_mesh_p4est->end ( tdim ); it != end_it; ++it )
            {
                // get vertex ids of current cell
                std::vector<Id> v_ids;
                for ( VertexIdIterator it_v = it->begin_vertex_ids ( ); it_v != it->end_vertex_ids ( ); ++it_v )
                {
                    v_ids.push_back ( *it_v );
                }
                assert ( v_ids.size ( ) == num_vert_per_cell );

                // put vertex ids in connectivity data
                std::vector<Id> v ( num_vert_per_cell );
                v[0] = v_ids[0];
                v[1] = v_ids[1];
                v[2] = v_ids[3];
                v[3] = v_ids[2];
                if ( tdim == 3 )
                {
                    v[4] = v_ids[4];
                    v[5] = v_ids[5];
                    v[6] = v_ids[7];
                    v[7] = v_ids[6];
                }

                p4est_tree_to_vertices.insert ( p4est_tree_to_vertices.end ( ), v.begin ( ), v.end ( ) );

                // store entitiy(tdim)-to-quad map
                // [tree_id, ref_level, x, y, z, (as topological coordinates inside of cell), localId, topological dim, id in coarse mesh]

                QuadCoord coord ( it->index ( ), 0, 0, 0, 0, 0, tdim, it->id ( ) );
                refined_mesh_p4est->get_db ( )->add_entity_to_quad_coord ( tdim, it->id ( ), coord, 0 );
                LOG_DEBUG ( 2, "added coarse entity_to_quad map for cell id " << it->id ( ) << " with corresponding tree id " << it->index ( ) );

                assert ( it->index ( ) == counter );
                counter++;
            }
            refined_mesh_p4est->get_db ( )->set_conn_data ( num_vertices, num_cells, p4est_vertices, p4est_tree_to_vertices );

            return refined_mesh;
        }

        // ---------------------------------------------------------------------

        void MeshPXest::set_cell_flags ( std::vector<bool>& flags )
        {
            this->cell_is_local_ = flags;
        }

        MeshPtr MeshPXest::add_ghost_cells ( MPI_Comm comm, int layer_width ) const
        {
            assert ( layer_width > 0 );

            LOG_DEBUG ( 1, "Compute ghost cells with lyer width " << layer_width );
            int rank = -1;
            MPI_Comm_rank ( comm, &rank );
            const TDim tdim = this->tdim ( );

            // get Ids of old ghost cells
            SortedArray<Id> old_ghosts = this->get_ghost_ids ( );

            // Initialize p4est ghost layer
            this->init_ghost_layer ( layer_width );

            // Pack data correspodning to mirror quads 
            this->init_ghost_data ( );

            // exchange data in ghost layers and put it into database, sub_domains and remote_indices
            std::vector< SortedArray<Id> > new_cells;
            this->exchange_ghost_data ( new_cells );

            // Add new ghost cells to mesh
            SortedArray<Id> fine_cells = this->get_entities ( this->tdim ( ) );
            for ( int l = 0; l < new_cells[this->history_index_].size ( ); ++l )
            {
                fine_cells.find_insert ( new_cells[this->history_index_][l] );
            }

            // TODOS :: zug. cell_to_mesh_map entfernen oder überschreiben
            // Remove outdated ghost cells from current mesh
            for ( int l = 0; l < old_ghosts.size ( ); ++l )
            {
                Id ghost_id = old_ghosts[l];
                int pos;
                if ( !new_cells[this->history_index_].find ( ghost_id, &pos ) )
                {
                    if ( fine_cells.find ( ghost_id, &pos ) )
                    {
                        fine_cells.erase ( pos );
                    }
                }
            }

            std::vector<bool> cell_is_local ( fine_cells.size ( ), true );
            for ( int l = 0; l < fine_cells.size ( ); ++l )
            {
                Id cur_cell = fine_cells.data ( ).at ( l );
                int pos;
                if ( new_cells[this->history_index_].find ( cur_cell, &pos ) )
                {
                    cell_is_local[l] = false;
                }
            }

            LOG_DEBUG ( 3, "[" << rank << "] Mesh " << this->history_index_
                        << " : cell ids after adding ghost cells: " << string_from_range ( fine_cells.begin ( ), fine_cells.end ( ) ) );

            MeshPtr mesh_with_ghost = MeshPtr ( new MeshPXest ( this->tdim ( ),
                                                                this->gdim ( ),
                                                                this->pXest_db_,
                                                                fine_cells,
                                                                this->history_index_,
                                                                cell_is_local,
                                                                this->period_ ) );

            LOG_DEBUG ( 1, "[" << rank << "] Number of cells after adding ghost layer: " << mesh_with_ghost->num_entities ( tdim ) );

            //boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( mesh_with_ghost );                                                    
            //mesh_pXest->set_cell_flags(cell_is_local);

            // Update pointer to this mesh in database                                                 
            this->pXest_db_->set_mesh ( mesh_with_ghost, this->history_index_ );

            assert ( this->pXest_db_->get_mesh ( this->history_index_ ) == mesh_with_ghost );
            assert ( this->pXest_db_->get_mesh ( this->history_index_ )->num_entities ( this->tdim ( ) ) == fine_cells.size ( ) );

            // update old meshes if new ghost cells were included        
            assert ( new_cells.size ( ) == this->history_index_ + 1 );
            for ( int hi = 0; hi<this->history_index_; ++hi )
            {
                if ( new_cells[hi].size ( ) > 0 )
                {
                    boost::intrusive_ptr<MeshPXest> cur_mesh = boost::static_pointer_cast<MeshPXest> ( this->pXest_db_->get_mesh ( hi ) );
                    LOG_DEBUG ( 2, " Mesh " << hi << " gets possibly " << new_cells[hi].size ( ) << " new ghost cells " );

                    cur_mesh->add_cells_and_rebuild ( new_cells[hi] );
                    for ( int l = 0; l < new_cells[hi].size ( ); ++l )
                    {
                        LOG_DEBUG ( 3, " Mesh: " << hi << ", new cell: " << new_cells[hi][l] );
                    }
                }
            }

            // update remote cell indices    
            std::vector<SubDomainId> sub_domains;
            std::vector<EntityNumber> remote_indices;

            boost::static_pointer_cast<MeshPXest>( mesh_with_ghost )->init_ghost_indices ( );
            boost::static_pointer_cast<MeshPXest>( mesh_with_ghost )->exchange_ghost_indices ( );
            boost::static_pointer_cast<MeshPXest>( mesh_with_ghost )->compute_subdomains_and_remote_indices ( sub_domains, remote_indices );

            // add attributes to mesh 
            AttributePtr sub_domain_attr = AttributePtr ( new IntAttribute ( sub_domains ) );
            mesh_with_ghost->add_attribute ( "_sub_domain_", tdim, sub_domain_attr );
            AttributePtr remote_index_attr = AttributePtr ( new IntAttribute ( remote_indices ) );
            mesh_with_ghost->add_attribute ( "_remote_index_", tdim, remote_index_attr );

            return mesh_with_ghost;
        }

        SortedArray<Id> MeshPXest::get_ghost_ids ( ) const
        {
            SortedArray<Id> ghost_ids;
            const int tdim = this->tdim ( );

#ifdef WITH_P4EST
            p4est_ghost_t* ghost4 = this->pXest_db_->get_p4est_ghost ( );
            p4est_t* forest4 = this->pXest_db_->get_p4est_forest ( );
            p8est_ghost_t* ghost8 = this->pXest_db_->get_p8est_ghost ( );
            p8est_t* forest8 = this->pXest_db_->get_p8est_forest ( );
            int num_ghosts = 0;

            if ( tdim == 2 )
            {
                if ( ghost4 == NULL )
                {
                    return ghost_ids;
                }
                num_ghosts = ghost4->ghosts.elem_count;
            }
            else if ( tdim == 3 )
            {
                if ( ghost8 == NULL )
                {
                    return ghost_ids;
                }
                num_ghosts = ghost8->ghosts.elem_count;
            }

            // Loop over all ghost quadrants
            for ( int jq = 0; jq < num_ghosts; ++jq )
            {
                // get ghost quad data pointer
                QuadData* q_ptr;
                int level = -1;
                if ( tdim == 2 )
                {
                    p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost4->ghosts, jq );
                    q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                    level = ghost_quad->level;
                }
                else if ( tdim == 3 )
                {
                    p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost8->ghosts, jq );
                    q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                    level = ghost_quad->level;
                }

                const Id cell_id = q_ptr->get_cell_id ( level );

                assert ( cell_id >= 0 );
                ghost_ids.insert ( cell_id );
            }
#endif
            return ghost_ids;
        }

        void MeshPXest::init_ghost_layer ( int layer_width ) const
        {
            assert ( layer_width > 0 );
            this->ghost_layer_width_ = layer_width;

#ifdef WITH_P4EST
            LOG_DEBUG ( 1, "Init p4est ghost layer of width " << layer_width );
            if ( this->tdim ( ) == 2 )
            {
                // create new ghost layer
                p4est_ghost_t* ghost = p4est_ghost_new ( this->pXest_db_->get_p4est_forest ( ), P4EST_CONNECT_FULL );

                // expand ghost layer
                for ( int l = 1; l < layer_width; ++l )
                {
                    p4est_ghost_expand ( this->pXest_db_->get_p4est_forest ( ), ghost );
                }
                this->get_db ( )->set_p4est_ghost ( ghost );
            }
            if ( this->tdim ( ) == 3 )
            {
                // create new ghost layer
                p8est_ghost_t* ghost = p8est_ghost_new ( this->pXest_db_->get_p8est_forest ( ), P8EST_CONNECT_FULL );

                // expand ghost layer
                for ( int l = 1; l < layer_width; ++l )
                {
                    p8est_ghost_expand ( this->pXest_db_->get_p8est_forest ( ), ghost );
                }
                this->get_db ( )->set_p8est_ghost ( ghost );
            }

            this->get_db ( )->set_layer_width ( layer_width );
#endif        
        }

        void MeshPXest::init_ghost_data ( ) const
        {
#ifdef WITH_P4EST
            int tdim = this->tdim ( );
            int num_mirrors = 0;
            std::vector<EntityNumber> mirror_indices;

            // get mirror cells -> local cells that are in ghost layer of any other process
            if ( tdim == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                p4est_ghost_t* ghost = this->pXest_db_->get_p4est_ghost ( );
                num_mirrors = ghost->mirrors.elem_count;
                this->mirror_data_22_.resize ( num_mirrors );

                // Loop over all mirror quadrants
                for ( int jq = 0; jq < num_mirrors; ++jq )
                {
                    // get mirror quad
                    p4est_quadrant_t* mirror_quad = p4est_quadrant_array_index ( &ghost->mirrors, jq );
                    treeId jt = mirror_quad->p.piggy3.which_tree;

                    // extract quad from corresponding local tree
                    p4est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                    p4est_locidx_t which_quad = mirror_quad->p.piggy3.local_num - tree->quadrants_offset;
                    p4est_quadrant_t* mirror_quad_with_data = p4est_quadrant_array_index ( &tree->quadrants, which_quad );
                    QuadData* q_ptr = pXest_get_quad_data_ptr ( mirror_quad_with_data );

                    EntityNumber index = q_ptr->get_cell_index ( );
                    assert ( index >= 0 );

                    // Check if quadrants corresponds to active cell 
                    if ( index >= 0 )
                    {
                        // transfer QuadData and entity packages
                        mirror_indices.push_back ( index );

                        // pack cell and facet entities of mirror cell and all its ancestors                        
                        int level = mirror_quad_with_data->level;
                        std::vector<EntityPackage> cell_packs ( level + 1 );
                        std::vector<EntityPackage> facet_packs ( ( level + 1 )*4 );
                        int facet_pack_offset = 0;

                        GhostCommData22* mirror_pack = new GhostCommData22;

                        // Loop over ancestors
                        for ( int l = level; l >= 0; --l )
                        {
                            Id cell_id = q_ptr->get_cell_id ( l );

                            // Put most recent mesh history index of cell into data package
                            int history_index = this->pXest_db_->get_last_mesh_history_index ( cell_id );
                            mirror_pack->mesh_history[l] = history_index;

                            LOG_DEBUG ( 3, "Pack cell id " << cell_id << " on level " << l << " with history index " << history_index );
                            // Pack cell entity
                            SortedArray<Id> ids ( 1 );
                            ids[0] = cell_id;
                            this->pXest_db_->create_entity_package ( tdim, ids, &cell_packs[l] );

                            // pack incident facet entities
                            SortedArray<Id> incident_facets;
                            Connectivity cell_d_connections;
                            this->pXest_db_->build ( tdim - 1, ids, incident_facets, cell_d_connections );

                            for ( int jf = 0; jf < incident_facets.size ( ); ++jf )
                            {
                                std::vector<Id> facet_ids ( 1 );
                                facet_ids[0] = incident_facets[jf];
                                this->pXest_db_->create_entity_package ( tdim - 1, facet_ids, &facet_packs[facet_pack_offset] );
                                facet_pack_offset++;
                            }
                        }

                        // put packs into mirror data
                        pack_ghost_data ( mirror_pack, *q_ptr, cell_packs, facet_packs );

                        this->mirror_data_22_[jq] = mirror_pack;
                    }
                        // BUG_TODO else nötig?
                    else
                    {
                        // transfer only QuadData
                        GhostCommData22* mirror_pack = new GhostCommData22;
                        pack_ghost_data ( mirror_pack, *q_ptr );
                        this->mirror_data_22_[jq] = mirror_pack;
                    }
                }
            }
            else if ( tdim == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                p8est_ghost_t* ghost = this->pXest_db_->get_p8est_ghost ( );
                num_mirrors = ghost->mirrors.elem_count;
                this->mirror_data_33_.resize ( num_mirrors );

                // Loop over all mirror quadrants
                for ( int jq = 0; jq < num_mirrors; ++jq )
                {
                    // get mirror quad
                    p8est_quadrant_t* mirror_quad = p8est_quadrant_array_index ( &ghost->mirrors, jq );
                    treeId jt = mirror_quad->p.piggy3.which_tree;

                    // extract quad from corresponding local tree
                    p8est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                    p4est_locidx_t which_quad = mirror_quad->p.piggy3.local_num - tree->quadrants_offset;
                    p8est_quadrant_t* mirror_quad_with_data = p8est_quadrant_array_index ( &tree->quadrants, which_quad );
                    QuadData* q_ptr = pXest_get_quad_data_ptr ( mirror_quad_with_data );

                    EntityNumber index = q_ptr->get_cell_index ( );
                    assert ( index >= 0 );

                    // Check if quadrants corresponds to active cell 
                    if ( index >= 0 )
                    {
                        // transfer QuadData and entity packages
                        mirror_indices.push_back ( index );

                        // pack cell and facet entities of mirror cell and all its ancestors                        
                        int level = mirror_quad_with_data->level;
                        std::vector<EntityPackage> cell_packs ( level + 1 );
                        std::vector<EntityPackage> facet_packs ( ( level + 1 )*6 );
                        int facet_pack_offset = 0;

                        GhostCommData33* mirror_pack = new GhostCommData33;

                        // Loop over ancestors
                        for ( int l = level; l >= 0; --l )
                        {
                            Id cell_id = q_ptr->get_cell_id ( l );

                            // Put most recent mesh history index of cell into data package
                            int history_index = this->pXest_db_->get_last_mesh_history_index ( cell_id );
                            mirror_pack->mesh_history[l] = history_index;

                            LOG_DEBUG ( 3, "Pack cell id " << cell_id << " on level " << l << " with history index " << history_index );
                            // Pack cell entity
                            SortedArray<Id> ids ( 1 );
                            ids[0] = cell_id;
                            this->pXest_db_->create_entity_package ( tdim, ids, &cell_packs[l] );

                            // pack incident facet entities
                            SortedArray<Id> incident_facets;
                            Connectivity cell_d_connections;
                            this->pXest_db_->build ( tdim - 1, ids, incident_facets, cell_d_connections );

                            for ( int jf = 0; jf < incident_facets.size ( ); ++jf )
                            {
                                std::vector<Id> facet_ids ( 1 );
                                facet_ids[0] = incident_facets[jf];
                                this->pXest_db_->create_entity_package ( tdim - 1, facet_ids, &facet_packs[facet_pack_offset] );
                                facet_pack_offset++;
                            }
                        }

                        // put packs into mirror data
                        pack_ghost_data ( mirror_pack, *q_ptr, cell_packs, facet_packs );

                        this->mirror_data_33_[jq] = mirror_pack;
                    }
                        // BUG_TODO else nötig?
                    else
                    {
                        // transfer only QuadData
                        GhostCommData33* mirror_pack = new GhostCommData33;
                        pack_ghost_data ( mirror_pack, *q_ptr );
                        this->mirror_data_33_[jq] = mirror_pack;
                    }
                }
            }
#endif
        }

        void MeshPXest::exchange_ghost_data ( std::vector< SortedArray<Id> >& new_cells ) const
        {
            LOG_DEBUG ( 1, "Exchange ghost cells " );
            MeshPXestBuilder builder ( this->tdim ( ), this->gdim ( ), this->get_db ( ), this->get_entities ( this->tdim ( ) ) );
#ifdef WITH_P4EST
            int tdim = this->tdim ( );
            int num_ghosts;

            // exchange ghost data with p4est and retrieve entity_packages from ghost quadrants
            std::vector<EntityPackage> cell_packs;
            new_cells.resize ( this->history_index_ + 1 );

            if ( tdim == 2 )
            {
                p4est_ghost_t* ghost = this->pXest_db_->get_p4est_ghost ( );
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );

                // init some data structures
                num_ghosts = ghost->ghosts.elem_count;
                GhostCommData22* ghost_data_22;
                ghost_data_22 = ( GhostCommData22* ) malloc ( num_ghosts * sizeof (GhostCommData22 ) );

                p4est_ghost_exchange_custom ( forest, ghost, sizeof (GhostCommData22 ), ( void** ) & this->mirror_data_22_[0], ( void* ) ghost_data_22 );
                num_ghosts = ghost->ghosts.elem_count;

                // Loop over all ghost quadrants
                for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                {
                    // get ghost quad
                    p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost->ghosts, jq );
                    treeId jt = ghost_quad->p.piggy3.which_tree;

                    // get pointer to data of current ghost quad
                    GhostCommData22* g_ptr = &ghost_data_22[jq];
                    EntityNumber remote_index = g_ptr->cell_index_remote;

                    // check if quadrant corresponds to cell in current mesh
                    if ( remote_index >= 0 )
                    {
                        // extract data from ghost_data
                        QuadData* q_ptr = new QuadData ( );
                        std::vector<EntityPackage> cell_packs;
                        std::vector<EntityPackage> facet_packs;

                        unpack_ghost_data ( *g_ptr, *q_ptr, cell_packs, facet_packs );

                        // put QuadData into ghost quadrant
                        pXest_set_quad_data_ptr ( ghost_quad, q_ptr );
                        int num_cells = cell_packs.size ( );
                        int level = ghost_quad->level;
                        std::vector<Id> ghost_ids;

                        assert ( level + 1 == num_cells );

                        // put cells into data base    
                        for ( int l = 0; l <= level; ++l )
                        {
                            std::vector<Id> new_cell_id = unpack_entities ( cell_packs[l], builder );
                            q_ptr->set_cell_id ( new_cell_id[0], l );

                            // Put history indices into database
                            std::vector<int> history_index ( 1, g_ptr->mesh_history[l] );
                            assert ( history_index[0] <= this->history_index_ );

                            new_cells[history_index[0]].find_insert ( new_cell_id[0] );
                            ghost_ids.push_back ( new_cell_id[0] );
                            this->pXest_db_->add_cell_to_mesh_maps ( new_cell_id, history_index );

                            LOG_DEBUG ( 3, "UnPack cell id " << new_cell_id[0] << " on level " << l << " with history index " << history_index[0] );
                        }

                        // put facets into database
                        for ( int l = 0; l < facet_packs.size ( ); ++l )
                        {
                            unpack_entities ( facet_packs[l], builder );
                        }

                        // add new entity_to_quad maps
                        p4est_quadrant_t ancestor_quad;
                        ancestor_quad.x = ghost_quad->x;
                        ancestor_quad.y = ghost_quad->y;
                        ancestor_quad.level = ghost_quad->level;

                        for ( int l = level; l >= 0; --l )
                        {
                            QuadCoord coord = pXest_get_quad_coord ( &ancestor_quad, jt );
                            this->pXest_db_->add_entity_to_quad_coord ( tdim, ghost_ids[l], coord, 1 );

                            if ( l >= 1 )
                            {
                                p4est_quadrant_t tmp_quad;
                                p4est_quadrant_parent ( &ancestor_quad, &tmp_quad );
                                ancestor_quad.x = tmp_quad.x;
                                ancestor_quad.y = tmp_quad.y;
                                ancestor_quad.level = tmp_quad.level;
                            }
                        }
                    }
                }

                // free ghost_data
                free ( ghost_data_22 );
            }
            else if ( tdim == 3 )
            {
                p8est_ghost_t* ghost = this->pXest_db_->get_p8est_ghost ( );
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );

                // init some data structures
                num_ghosts = ghost->ghosts.elem_count;
                GhostCommData33* ghost_data_33;
                ghost_data_33 = ( GhostCommData33* ) malloc ( num_ghosts * sizeof (GhostCommData33 ) );

                p8est_ghost_exchange_custom ( forest, ghost, sizeof (GhostCommData33 ), ( void** ) & this->mirror_data_33_[0], ( void* ) ghost_data_33 );
                num_ghosts = ghost->ghosts.elem_count;

                // Loop over all ghost quadrants
                for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                {
                    // get ghost quad
                    p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost->ghosts, jq );
                    treeId jt = ghost_quad->p.piggy3.which_tree;

                    // get pointer to data of current ghost quad
                    GhostCommData33* g_ptr = &ghost_data_33[jq];
                    EntityNumber remote_index = g_ptr->cell_index_remote;

                    // check if quadrant corresponds to cell in current mesh
                    if ( remote_index >= 0 )
                    {
                        // extract data from ghost_data
                        QuadData* q_ptr = new QuadData ( );
                        std::vector<EntityPackage> cell_packs;
                        std::vector<EntityPackage> facet_packs;

                        unpack_ghost_data ( *g_ptr, *q_ptr, cell_packs, facet_packs );

                        // put QuadData into ghost quadrant
                        pXest_set_quad_data_ptr ( ghost_quad, q_ptr );
                        int num_cells = cell_packs.size ( );
                        int level = ghost_quad->level;
                        std::vector<Id> ghost_ids;

                        assert ( level + 1 == num_cells );

                        // put cells into data base    
                        for ( int l = 0; l <= level; ++l )
                        {
                            std::vector<Id> new_cell_id = unpack_entities ( cell_packs[l], builder );
                            q_ptr->set_cell_id ( new_cell_id[0], l );

                            // Put history indices into database
                            std::vector<int> history_index ( 1, g_ptr->mesh_history[l] );
                            assert ( history_index[0] <= this->history_index_ );

                            new_cells[history_index[0]].find_insert ( new_cell_id[0] );
                            ghost_ids.push_back ( new_cell_id[0] );
                            this->pXest_db_->add_cell_to_mesh_maps ( new_cell_id, history_index );

                            LOG_DEBUG ( 3, "UnPack cell id " << new_cell_id[0] << " on level " << l << " with history index " << history_index[0] );
                        }

                        // put facets into database
                        for ( int l = 0; l < facet_packs.size ( ); ++l )
                        {
                            unpack_entities ( facet_packs[l], builder );
                        }

                        // add new entity_to_quad maps
                        p8est_quadrant_t ancestor_quad;
                        ancestor_quad.x = ghost_quad->x;
                        ancestor_quad.y = ghost_quad->y;
                        ancestor_quad.z = ghost_quad->z;
                        ancestor_quad.level = ghost_quad->level;

                        for ( int l = level; l >= 0; --l )
                        {
                            QuadCoord coord = pXest_get_quad_coord ( &ancestor_quad, jt );
                            this->pXest_db_->add_entity_to_quad_coord ( tdim, ghost_ids[l], coord, 1 );

                            if ( l >= 1 )
                            {
                                p8est_quadrant_t tmp_quad;
                                p8est_quadrant_parent ( &ancestor_quad, &tmp_quad );
                                ancestor_quad.x = tmp_quad.x;
                                ancestor_quad.y = tmp_quad.y;
                                ancestor_quad.z = tmp_quad.z;
                                ancestor_quad.level = tmp_quad.level;
                            }
                        }
                    }
                }

                // free ghost_data
                free ( ghost_data_33 );
            }
#endif
        }

        void MeshPXest::init_ghost_indices ( ) const
        {
#ifdef WITH_P4EST
            int tdim = this->tdim ( );
            int num_mirrors = 0;
            std::vector<EntityNumber> mirror_indices;

            p4est_ghost_t* ghost4 = this->pXest_db_->get_p4est_ghost ( );
            p4est_t* forest4 = this->pXest_db_->get_p4est_forest ( );
            p8est_ghost_t* ghost8 = this->pXest_db_->get_p8est_ghost ( );
            p8est_t* forest8 = this->pXest_db_->get_p8est_forest ( );

            // get mirror cells -> local cells that are in ghost layer of any other process
            if ( tdim == 2 )
            {
                assert ( ghost4 != NULL );
                num_mirrors = ghost4->mirrors.elem_count;
            }
            else if ( tdim == 3 )
            {
                assert ( ghost8 != NULL );
                num_mirrors = ghost8->mirrors.elem_count;
            }

            this->mirror_data_.resize ( num_mirrors );

            // Loop over all mirror quadrants
            for ( int jq = 0; jq < num_mirrors; ++jq )
            {
                // BUG_TODO take care of level? 

                QuadData* q_ptr;

                // get mirror quad data pointer
                if ( tdim == 2 )
                {
                    p4est_quadrant_t* mirror_quad = p4est_quadrant_array_index ( &ghost4->mirrors, jq );
                    treeId jt = mirror_quad->p.piggy3.which_tree;

                    // extract quad from corresponding local tree
                    p4est_tree_t* tree = pXest_get_tree_in_forest ( forest4, jt );
                    p4est_locidx_t which_quad = mirror_quad->p.piggy3.local_num - tree->quadrants_offset;
                    p4est_quadrant_t* mirror_quad_with_data = p4est_quadrant_array_index ( &tree->quadrants, which_quad );
                    q_ptr = pXest_get_quad_data_ptr ( mirror_quad_with_data );
                }
                else if ( tdim == 3 )
                {
                    p8est_quadrant_t* mirror_quad = p8est_quadrant_array_index ( &ghost8->mirrors, jq );
                    treeId jt = mirror_quad->p.piggy3.which_tree;

                    // extract quad from corresponding local tree
                    p8est_tree_t* tree = pXest_get_tree_in_forest ( forest8, jt );
                    p4est_locidx_t which_quad = mirror_quad->p.piggy3.local_num - tree->quadrants_offset;
                    p8est_quadrant_t* mirror_quad_with_data = p8est_quadrant_array_index ( &tree->quadrants, which_quad );
                    q_ptr = pXest_get_quad_data_ptr ( mirror_quad_with_data );
                }

                EntityNumber index = q_ptr->get_cell_index ( );

                GhostCommData* mirror_pack = new GhostCommData;
                pack_ghost_data ( mirror_pack, *q_ptr );
                this->mirror_data_[jq] = mirror_pack;
            }
#endif
        }

        void MeshPXest::exchange_ghost_indices ( ) const
        {
#ifdef WITH_P4EST
            int tdim = this->tdim ( );

            p4est_ghost_t* ghost4 = this->pXest_db_->get_p4est_ghost ( );
            p4est_t* forest4 = this->pXest_db_->get_p4est_forest ( );
            p8est_ghost_t* ghost8 = this->pXest_db_->get_p8est_ghost ( );
            p8est_t* forest8 = this->pXest_db_->get_p8est_forest ( );
            int num_ghosts = 0;
            GhostCommData* ghost_data;

            // init some data structures and exchange data between ghosts
            if ( tdim == 2 )
            {
                assert ( ghost4 != NULL );

                num_ghosts = ghost4->ghosts.elem_count;
                ghost_data = ( GhostCommData* ) malloc ( num_ghosts * sizeof (GhostCommData ) );
                p4est_ghost_exchange_custom ( forest4, ghost4, sizeof (GhostCommData ), ( void** ) & this->mirror_data_[0], ( void* ) ghost_data );
            }
            if ( tdim == 3 )
            {
                assert ( ghost8 != NULL );

                num_ghosts = ghost8->ghosts.elem_count;
                ghost_data = ( GhostCommData* ) malloc ( num_ghosts * sizeof (GhostCommData ) );
                p8est_ghost_exchange_custom ( forest8, ghost8, sizeof (GhostCommData ), ( void** ) & this->mirror_data_[0], ( void* ) ghost_data );
            }

            // Loop over all ghost quadrants
            for ( int jq = 0; jq < num_ghosts; ++jq )
            {
                // get ghost quad data pointer
                QuadData* q_ptr;

                if ( tdim == 2 )
                {
                    p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost4->ghosts, jq );
                    q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                }
                else if ( tdim == 3 )
                {
                    p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost8->ghosts, jq );
                    q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                }

                GhostCommData* g_ptr = &ghost_data[jq];
                EntityNumber index = g_ptr->cell_index;

                // put cell index into ghost quadrant
                q_ptr->set_remote_cell_index ( index );
            }

            // free ghost_data
            free ( ghost_data );
#endif
        }

        int MeshPXest::get_rank ( ) const
        {
#ifdef WITH_P4EST
            if ( this->tdim ( ) == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                return forest->mpirank;
            }
            if ( this->tdim ( ) == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                return forest->mpirank;
            }
#endif
        }

        void MeshPXest::compute_subdomains_and_remote_indices ( std::vector<SubDomainId>& sub_domains, std::vector<EntityNumber>& remote_indices ) const
        {
#ifdef WITH_P4EST
            int tdim = this->tdim ( );
            int my_rank = this->get_rank ( );

            int num_cells = this->num_entities ( tdim );
            sub_domains .resize ( num_cells, my_rank );
            remote_indices.resize ( num_cells, -1 );
            std::vector<QuadData*> q_ptrs;

            if ( tdim == 2 )
            {
                p4est_ghost_t* ghost = this->pXest_db_->get_p4est_ghost ( );

                if ( ghost != NULL )
                {
                    // Loop over all ghost quads
                    for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                    {
                        // get data from ghost quad
                        p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost->ghosts, jq );
                        q_ptrs.push_back ( pXest_get_quad_data_ptr ( ghost_quad ) );
                    }
                }
            }
            if ( tdim == 3 )
            {
                p8est_ghost_t* ghost = this->pXest_db_->get_p8est_ghost ( );

                if ( ghost != NULL )
                {
                    // Loop over all ghost quads
                    for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                    {
                        // get data from ghost quad
                        p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost->ghosts, jq );
                        q_ptrs.push_back ( pXest_get_quad_data_ptr ( ghost_quad ) );
                    }
                }
            }

            // extract data from quads
            for ( int j = 0; j < q_ptrs.size ( ); ++j )
            {
                EntityNumber index = q_ptrs[j]->get_cell_index ( );
                EntityNumber remote_index = q_ptrs[j]->get_remote_cell_index ( );
                int owner_rank = q_ptrs[j]->get_owner_rank ( );

                assert ( index >= 0 && index < num_cells );
                assert ( remote_index >= 0 );

                sub_domains[index] = owner_rank;
                remote_indices[index] = remote_index;

                LOG_DEBUG ( 3, "Index: " << index << " remote index: " << remote_index << " subdomain: " << owner_rank << " my rank " << my_rank );

                assert ( remote_index == -1 || ( owner_rank != my_rank ) );
                assert ( ( owner_rank == my_rank ) || remote_index > -1 );
            }
#endif
        }

        void MeshPXest::find_coarsenable_cells ( SortedArray<EntityNumber>& coarsenable_cells ) const
        {
            coarsenable_cells.clear ( );
            ForestPatchData patch_data;
#ifdef WITH_P4EST
            if ( this->tdim ( ) == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                pXest_set_forest_user_ptr ( forest, &patch_data );
                p4est_coarsen_ext ( forest, 0, 0, pXest_coarsen_find_fn, pXest_coarsen_init_fn, pXest_coarsen_replace_fn );
            }
            else if ( this->tdim ( ) == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                pXest_set_forest_user_ptr ( forest, &patch_data );
                p8est_coarsen_ext ( forest, 0, 0, pXest_coarsen_find_fn, pXest_coarsen_init_fn, pXest_coarsen_replace_fn );
            }

            // PERF_TODO: use SortedArray copy operator
            for ( int l = 0; l < patch_data.coarsenable_cells.size ( ); ++l )
            {
                coarsenable_cells.find_insert ( patch_data.coarsenable_cells.data ( ).at ( l ) );
            }
#endif
        }

        MeshPtr MeshPXest::make_mesh_coarsenable ( )
        {
            LOG_DEBUG ( 1, "Make mesh coarsenable " );
            const int tdim = this->tdim ( );
            int num_cells = this->num_entities ( tdim );
            bool finished = false;
            MeshPtr tmp_mesh = this;
            boost::intrusive_ptr<MeshPXest> tmp_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( tmp_mesh );

            tmp_mesh_pXest->set_connection_mode ( 2 );
            tmp_mesh_pXest->set_patch_mode ( false );

            int iter = 1;

            MPI_Comm* mpi_comm;
#ifdef WITH_P4EST
            if ( this->tdim ( ) == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                mpi_comm = &forest->mpicomm;
            }
            else if ( this->tdim ( ) == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                mpi_comm = &forest->mpicomm;
            }
#endif
            while ( !finished )
            {
                // Find all cells that form a family which is coarsenable
                SortedArray<EntityNumber> coarsenable_cells;
                tmp_mesh_pXest->find_coarsenable_cells ( coarsenable_cells );

                // Find cells that are not coarsenable
                std::vector<EntityNumber> non_coarsenable_cells;
                for ( EntityNumber j = 0; j < num_cells; ++j )
                {
                    int pos;
                    if ( !coarsenable_cells.find ( j, &pos ) )
                    {
                        if ( tmp_mesh_pXest->cell_is_local ( j ) )
                        {
                            non_coarsenable_cells.push_back ( j );
                        }
                    }
                }

                LOG_DEBUG ( 1, "Iteration: " << iter << " total number of cells: " << num_cells << " coarsenable cells: " << coarsenable_cells.size ( )
                            << " non coarsenable local cells: " << non_coarsenable_cells.size ( ) );

                int local_num_non_coarsenable_cells = non_coarsenable_cells.size ( );
                int global_num_non_coarsenable_cells = 0;
                MPI_Allreduce ( &local_num_non_coarsenable_cells, &global_num_non_coarsenable_cells,
                                1, MPI_INT, MPI_SUM, *mpi_comm );
                if ( global_num_non_coarsenable_cells == 0 )
                {
                    finished = true;
                    break;
                }

                // Mark non_coarsenable_cells for refinement
                std::vector<int> refinement ( num_cells, 0 );
                for ( int l = 0; l < non_coarsenable_cells.size ( ); ++l )
                {
                    refinement[non_coarsenable_cells[l]] = 1;
                }

                // Refine and update variables
                tmp_mesh = tmp_mesh->refine ( refinement );
                tmp_mesh_pXest = boost::static_pointer_cast<MeshPXest> ( tmp_mesh );
                tmp_mesh_pXest->set_connection_mode ( 2 );
                tmp_mesh_pXest->set_patch_mode ( false );
                num_cells = tmp_mesh->num_entities ( tdim );
                iter++;
            }
            return tmp_mesh;
        }

        // ---------------------------------------------------------------------

        Id MeshPXest::get_ancestor_cell_id ( Id cell_id, int level_diff ) const
        {
#ifdef WITH_P4EST
            assert ( level_diff > 0 );
            assert ( cell_id >= 0 );
            int pos;
            assert ( this->get_entities ( this->tdim ( ) ).find ( cell_id, &pos ) );

            const int tdim = this->tdim ( );

            // *********************************************
            // get quad in forst based on entity_to_quad mapping
            QuadCoord coord = pXest_db_-> get_entity_to_quad_coord ( tdim, cell_id, 1 );
            assert ( coord.tree >= 0 );

            int level = coord.level;
            int32_t x = coord.x;
            int32_t y = coord.y;
            int32_t z = coord.z;
            treeId tree_id = coord.tree;

            LOG_DEBUG ( 3, "id " << cell_id << " level " << level << " x " << x << " y " << y << " z " << z << " tree " << tree_id );
            if ( level < 1 )
            {
                return -1;
            }
            if ( tdim == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                bool quad_is_local = false;
                QuadData* q_ptr;

                assert ( pXest_quad_arrays_in_forest_sorted ( forest ) );

                // check if tree is local
                treeId first_tree = pXest_get_first_local_treeId ( forest );
                treeId last_tree = pXest_get_last_local_treeId ( forest );

                // create quadrant object for corresponding cell
                p4est_quadrant_t tmp_quad;
                tmp_quad.level = level;
                tmp_quad.x = x;
                tmp_quad.y = y;

                if ( tree_id >= first_tree && tree_id <= last_tree )
                {// tree is local
                    LOG_DEBUG ( 3, "Search local trees " );

                    // get tree
                    p4est_tree_t* tree = pXest_get_tree_in_forest ( pXest_db_->get_p4est_forest ( ), tree_id );

                    q_ptr = pXest_get_quad_data_ptr ( tree, &tmp_quad );
                    if ( q_ptr != NULL )
                    {// tmp_quad contains required data
                        if ( q_ptr->get_cell_id ( level ) == cell_id )
                        {
                            LOG_DEBUG ( 3, "Requested cells corresponds to local leaf " );
                            quad_is_local = true;
                        }
                    }
                    else
                    {// either tmp_quad is in the ghost layer, or tmp_quad is not a leaf, i.e. it contains no data
                        // try to find descendant of tmp_quad which is a local leaf 

                        // Loop over all quads in current tree
                        int num_quads = tree->quadrants.elem_count;
                        for ( int jq = 0; jq < num_quads; ++jq )
                        {
                            p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, jq );
                            // check if cur_quad is descendant of tmp_quad
                            if ( pXest_quad_is_ancestor ( &tmp_quad, cur_quad ) )
                            {
                                LOG_DEBUG ( 3, "Requested cells corresponds to ancestor of local leaf " );
                                q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                                quad_is_local = true;
                                break;
                            }
                        }
                    }
                }
                if ( !quad_is_local )
                { // search ghost layer
                    LOG_DEBUG ( 3, "Search ghost layer " );

                    p4est_ghost_t* ghost = pXest_db_->get_p4est_ghost ( );

                    q_ptr = pXest_get_quad_data_ptr ( ghost, &tmp_quad, tree_id );

                    if ( q_ptr == NULL )
                    {
                        for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                        {
                            // get ghost quad and data
                            p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost->ghosts, jq );
                            q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                            assert ( q_ptr != NULL );

                            // test if ghost quadrant is correct quad
                            if ( q_ptr->get_cell_id ( level ) == cell_id )
                            {
                                LOG_DEBUG ( 3, "Requested cell found in ghost layer " );
                                /*
                                if (DEBUG_LEVEL > 0)
                                {
                                    if (!p4est_quadrant_is_equal(ghost_quad, &tmp_quad))
                                    {
                                        LOG_DEBUG(1, "x : " << ghost_quad->x << " <-> " << tmp_quad.x );
                                        LOG_DEBUG(1, "y : " << ghost_quad->y << " <-> " << tmp_quad.y );
                                        LOG_DEBUG(1, "level : " << (int)ghost_quad->level << " <-> " << (int)tmp_quad.level );
                                        q_ptr->print();
                                        coord.print();
                                    }
                                }
                                 * */
                                break;
                            }
                        }
                    }
                }

                assert ( q_ptr != NULL );
                /*
                if ( DEBUG_LEVEL >= 3 )
                {
                    LOG_DEBUG ( 3, "Data used for determination of ancestor id " );
                    q_ptr->print ( );
                }
                 */
                assert ( q_ptr->get_cell_id ( level ) == cell_id );

                // access user data of parent quad
                return q_ptr->get_cell_id ( level - level_diff );
            }
            else if ( tdim == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                bool quad_is_local = false;
                QuadData* q_ptr;

                assert ( pXest_quad_arrays_in_forest_sorted ( forest ) );

                // get mortonId of associated quadrant
                p8est_quadrant_t tmp_quad;
                tmp_quad.level = level;
                tmp_quad.x = x;
                tmp_quad.y = y;
                tmp_quad.z = z;

                // check if tree is local
                treeId first_tree = pXest_get_first_local_treeId ( forest );
                treeId last_tree = pXest_get_last_local_treeId ( forest );

                if ( tree_id >= first_tree && tree_id <= last_tree )
                {// tree is local
                    LOG_DEBUG ( 3, "Search local trees " );

                    // get tree
                    p8est_tree_t* tree = pXest_get_tree_in_forest ( pXest_db_->get_p8est_forest ( ), tree_id );

                    q_ptr = pXest_get_quad_data_ptr ( tree, &tmp_quad );
                    if ( q_ptr != NULL )
                    {// tmp_quad contains required data
                        if ( q_ptr->get_cell_id ( level ) == cell_id )
                        {
                            LOG_DEBUG ( 3, "Requested cells corresponds to local leaf " );
                            quad_is_local = true;
                        }
                    }
                    else
                    {// either tmp_quad is in the ghost layer, or tmp_quad is not a leaf, i.e. it contains no data
                        // try to find descendant of tmp_quad which is a local leaf 

                        // Loop over all quads in current tree
                        int num_quads = tree->quadrants.elem_count;
                        for ( int jq = 0; jq < num_quads; ++jq )
                        {
                            p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, jq );
                            // check if cur_quad is descendant of tmp_quad
                            if ( pXest_quad_is_ancestor ( &tmp_quad, cur_quad ) )
                            {
                                LOG_DEBUG ( 3, "Requested cells corresponds to ancestor of local leaf " );
                                q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                                quad_is_local = true;
                                break;
                            }
                        }
                    }
                }
                if ( !quad_is_local )
                { // search ghost layer
                    LOG_DEBUG ( 3, "Search ghost layer " );

                    p8est_ghost_t* ghost = pXest_db_->get_p8est_ghost ( );

                    q_ptr = pXest_get_quad_data_ptr ( ghost, &tmp_quad, tree_id );

                    if ( q_ptr == NULL )
                    {
                        for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                        {
                            // get ghost quad and data
                            p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost->ghosts, jq );
                            q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                            assert ( q_ptr != NULL );

                            // test if ghost quadrant is correct quad
                            if ( q_ptr->get_cell_id ( level ) == cell_id )
                            {
                                LOG_DEBUG ( 3, "Requested cell is ancestor of ghost cell " );
                                /*
                                if (DEBUG_LEVEL > 0)
                                {
                                    if (!p8est_quadrant_is_equal(ghost_quad, &tmp_quad))
                                    {
                                        LOG_DEBUG(1, "x : " << ghost_quad->x << " <-> " << tmp_quad.x );
                                        LOG_DEBUG(1, "y : " << ghost_quad->y << " <-> " << tmp_quad.y );
                                        LOG_DEBUG(1, "z : " << ghost_quad->z << " <-> " << tmp_quad.z );
                                        LOG_DEBUG(1, "level : " << (int)ghost_quad->level << " <-> " << (int)tmp_quad.level );
                                        q_ptr->print();
                                        coord.print();
                                    }
                                }
                                 * */
                                break;
                            }
                        }
                    }
                }

                assert ( q_ptr != NULL );
                /*
                if ( DEBUG_LEVEL >= 3 )
                {
                    LOG_DEBUG ( 3, "Data used for determination of ancestor id " );
                    q_ptr->print ( );
                }
                 */
                assert ( q_ptr->get_cell_id ( level ) == cell_id );

                // access user data of parent quad
                return q_ptr->get_cell_id ( level - level_diff );
            }
#else
            return -1;
#endif
        }

        std::vector<Id> MeshPXest::get_descendant_cell_ids ( Id cell_id, int level_diff ) const
        {
            assert ( cell_id >= 0 );
            assert ( level_diff > 0 );

            const int tdim = this->tdim ( );
            SortedArray<Id> children_ids;
#ifdef WITH_P4EST
            // *********************************************
            // get quad in forest based on entity_to_quad mapping
            QuadCoord coord = pXest_db_-> get_entity_to_quad_coord ( tdim, cell_id, 1 );
            assert ( coord.tree >= 0 );

            int level = coord.level;
            int32_t x = coord.x;
            int32_t y = coord.y;
            int32_t z = coord.z;
            treeId tree_id = coord.tree;

            LOG_DEBUG ( 3, "id " << cell_id << " level " << level << " x " << x << " y " << y << " z " << z << " tree " << tree_id );
            if ( tdim == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                p4est_ghost_t* ghost = this->pXest_db_->get_p4est_ghost ( );
                bool quad_is_local = false;
                QuadData* q_ptr;

                assert ( pXest_quad_arrays_in_forest_sorted ( forest ) );

                // create associated quadrant
                p4est_quadrant_t tmp_quad;
                tmp_quad.level = level;
                tmp_quad.x = x;
                tmp_quad.y = y;

                // check if tree is local
                treeId first_tree = pXest_get_first_local_treeId ( forest );
                treeId last_tree = pXest_get_last_local_treeId ( forest );

                if ( tree_id >= first_tree && tree_id <= last_tree )
                {// tree is local

                    // get tree
                    p4est_tree_t* tree = pXest_get_tree_in_forest ( pXest_db_->get_p4est_forest ( ), tree_id );

                    // check if tmp_quad is leaf in tree
                    int64_t pos = pXest_get_quad_pos_in_tree_array ( tree, &tmp_quad );
                    if ( pos >= 0 )
                    {// tmp_quad contains data -> is a leaf, i.e. has no children
                        LOG_DEBUG ( 3, " Cell is a leaf -> has no descendants, position " << pos );
                        return children_ids;
                    }

                    // search descendants of tmp_quad which are leafs
                    std::vector<int64_t> leaf_pos_in_tree = pXest_get_leaf_desc_in_tree_array ( tree, &tmp_quad );
                    LOG_DEBUG ( 3, " Search " << leaf_pos_in_tree.size ( ) << " leaf descendants" );
                    for ( int l = 0; l < leaf_pos_in_tree.size ( ); ++l )
                    {
                        p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, leaf_pos_in_tree[l] );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                        assert ( q_ptr != NULL );
                        Id desc_id = q_ptr->get_cell_id ( level + level_diff );
                        if ( desc_id >= 0 )
                        {
                            LOG_DEBUG ( 3, " Found descendant with id " << desc_id );
                            children_ids.find_insert ( desc_id );
                        }
                    }
                }
                // search ghost layer
                if ( ghost != NULL )
                {
                    // check if tmp_quad is leaf in ghost 
                    int64_t pos = pXest_get_quad_pos_in_ghost_array ( ghost, &tmp_quad, tree_id );

                    std::vector<int64_t> leaf_pos_in_ghost = pXest_get_leaf_desc_in_ghost_array ( ghost, &tmp_quad, tree_id );
                    LOG_DEBUG ( 3, " Search " << leaf_pos_in_ghost.size ( ) << " ghost descendants" );
                    for ( int l = 0; l < leaf_pos_in_ghost.size ( ); ++l )
                    {
                        p4est_quadrant_t* cur_quad = pXest_get_local_quad_in_ghost ( ghost, leaf_pos_in_ghost[l] );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                        assert ( q_ptr != NULL );
                        Id desc_id = q_ptr->get_cell_id ( level + level_diff );
                        if ( desc_id >= 0 )
                        {
                            LOG_DEBUG ( 3, " Found descendant with id " << desc_id );
                            children_ids.find_insert ( desc_id );
                        }
                    }
                }
            }
            else if ( tdim == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                p8est_ghost_t* ghost = this->pXest_db_->get_p8est_ghost ( );
                bool quad_is_local = false;
                QuadData* q_ptr;

                assert ( pXest_quad_arrays_in_forest_sorted ( forest ) );

                // create associated quadrant
                p8est_quadrant_t tmp_quad;
                tmp_quad.level = level;
                tmp_quad.x = x;
                tmp_quad.y = y;
                tmp_quad.z = z;

                // check if tree is local
                treeId first_tree = pXest_get_first_local_treeId ( forest );
                treeId last_tree = pXest_get_last_local_treeId ( forest );

                if ( tree_id >= first_tree && tree_id <= last_tree )
                {// tree is local

                    // get tree
                    p8est_tree_t* tree = pXest_get_tree_in_forest ( pXest_db_->get_p8est_forest ( ), tree_id );

                    // check if tmp_quad is leaf in tree
                    int64_t pos = pXest_get_quad_pos_in_tree_array ( tree, &tmp_quad );
                    if ( pos >= 0 )
                    {// tmp_quad contains data -> is a leaf, i.e. has no children
                        LOG_DEBUG ( 3, " Cell is a leaf -> has no descendants, position " << pos );
                        return children_ids;
                    }

                    // search descendants of tmp_quad which are leafs
                    std::vector<int64_t> leaf_pos_in_tree = pXest_get_leaf_desc_in_tree_array ( tree, &tmp_quad );
                    LOG_DEBUG ( 3, " Search " << leaf_pos_in_tree.size ( ) << " leaf descendants" );
                    for ( int l = 0; l < leaf_pos_in_tree.size ( ); ++l )
                    {
                        p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_tree ( tree, leaf_pos_in_tree[l] );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                        assert ( q_ptr != NULL );
                        Id desc_id = q_ptr->get_cell_id ( level + level_diff );
                        if ( desc_id >= 0 )
                        {
                            LOG_DEBUG ( 3, " Found descendant with id " << desc_id );
                            children_ids.find_insert ( desc_id );
                        }
                    }
                }
                // search ghost layer
                if ( ghost != NULL )
                {
                    // check if tmp_quad is leaf in ghost 
                    int64_t pos = pXest_get_quad_pos_in_ghost_array ( ghost, &tmp_quad, tree_id );

                    std::vector<int64_t> leaf_pos_in_ghost = pXest_get_leaf_desc_in_ghost_array ( ghost, &tmp_quad, tree_id );
                    LOG_DEBUG ( 3, " Search " << leaf_pos_in_ghost.size ( ) << " ghost descendants" );
                    for ( int l = 0; l < leaf_pos_in_ghost.size ( ); ++l )
                    {
                        p8est_quadrant_t* cur_quad = pXest_get_local_quad_in_ghost ( ghost, leaf_pos_in_ghost[l] );
                        QuadData* q_ptr = pXest_get_quad_data_ptr ( cur_quad );
                        assert ( q_ptr != NULL );
                        Id desc_id = q_ptr->get_cell_id ( level + level_diff );
                        if ( desc_id >= 0 )
                        {
                            LOG_DEBUG ( 3, " Found descendant with id " << desc_id );
                            children_ids.find_insert ( desc_id );
                        }
                    }
                }
            }
#endif
            return children_ids.data ( );
        }

        // ---------------------------------------------------------------------

        void MeshPXest::update_cell_index_in_quads ( bool update_ghost )
        {
#ifdef WITH_P4EST
            treeId first_treeId = -1;
            treeId last_treeId = -2;

            // update mapping id -> index for cells
            this->update_cell_id_to_index_map ( );

            if ( this->tdim ( ) == 2 )
            {
                p4est_t* forest = this->pXest_db_->get_p4est_forest ( );
                p4est_ghost_t* ghost = this->pXest_db_->get_p4est_ghost ( );

                if ( forest == NULL )
                {
                    return;
                }
                first_treeId = pXest_get_first_local_treeId ( forest );
                last_treeId = pXest_get_last_local_treeId ( forest );

                // Update local quadrants
                // Loop over all local trees
                for ( treeId jt = first_treeId; jt <= last_treeId; ++jt )
                {
                    p4est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                    int num_quads = tree->quadrants.elem_count;

                    // Loop over all local quads in tree
                    for ( int jq = 0; jq < num_quads; ++jq )
                    {
                        p4est_quadrant_t* quad = p4est_quadrant_array_index ( &tree->quadrants, jq );

                        QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
                        q_ptr->set_cell_index ( this->get_cell_index ( q_ptr->get_cell_id ( quad->level ) ) );
                        q_ptr->set_remote_cell_index ( this->get_cell_index ( q_ptr->get_cell_id ( quad->level ) ) );
                    }
                }
                if ( update_ghost && ghost != NULL )
                {
                    // Update ghost quadrants
                    // Loop over all ghost quadrants
                    for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                    {
                        // get ghost quad
                        p4est_quadrant_t* ghost_quad = p4est_quadrant_array_index ( &ghost->ghosts, jq );

                        QuadData* q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                        q_ptr->set_cell_index ( this->get_cell_index ( q_ptr->get_cell_id ( ghost_quad->level ) ) );
                    }
                }
            }
            else if ( this->tdim ( ) == 3 )
            {
                p8est_t* forest = this->pXest_db_->get_p8est_forest ( );
                p8est_ghost_t* ghost = this->pXest_db_->get_p8est_ghost ( );

                if ( forest == NULL )
                {
                    return;
                }
                first_treeId = pXest_get_first_local_treeId ( forest );
                last_treeId = pXest_get_last_local_treeId ( forest );

                // Update local quadrants
                // Loop over all local trees
                for ( treeId jt = first_treeId; jt <= last_treeId; ++jt )
                {
                    p8est_tree_t* tree = pXest_get_tree_in_forest ( forest, jt );
                    int num_quads = tree->quadrants.elem_count;

                    // Loop over all local quads in tree
                    for ( int jq = 0; jq < num_quads; ++jq )
                    {
                        p8est_quadrant_t* quad = p8est_quadrant_array_index ( &tree->quadrants, jq );

                        QuadData* q_ptr = pXest_get_quad_data_ptr ( quad );
                        q_ptr->set_cell_index ( this->get_cell_index ( q_ptr->get_cell_id ( quad->level ) ) );
                        q_ptr->set_remote_cell_index ( this->get_cell_index ( q_ptr->get_cell_id ( quad->level ) ) );
                    }
                }
                if ( update_ghost && ghost != NULL )
                {
                    // Update ghost quadrants
                    // Loop over all ghost quadrants
                    for ( int jq = 0; jq < ghost->ghosts.elem_count; ++jq )
                    {
                        // get ghost quad
                        p8est_quadrant_t* ghost_quad = p8est_quadrant_array_index ( &ghost->ghosts, jq );

                        QuadData* q_ptr = pXest_get_quad_data_ptr ( ghost_quad );
                        q_ptr->set_cell_index ( this->get_cell_index ( q_ptr->get_cell_id ( ghost_quad->level ) ) );
                    }
                }
            }
#endif
        }

        void MeshPXest::update_cell_id_to_index_map ( )
        {
            this->cell_id_to_index_.clear ( );
            for ( EntityNumber index = 0; index != this->num_entities ( this->tdim ( ) ); ++index )
            {
                this->cell_id_to_index_.insert ( std::pair<Id, EntityNumber> ( this->get_id ( this->tdim ( ), index ), index ) );
            }
        }

        int MeshPXest::get_sub_cell_number ( Id cell_id ) const
        {
#ifdef WITH_P4EST
            assert ( cell_id >= 0 );

            int tdim = this->tdim ( );

            // get quad in forst based on entity_to_quad mapping
            QuadCoord coord = pXest_db_-> get_entity_to_quad_coord ( tdim, cell_id, 1 );
            assert ( coord.tree >= 0 );

            int level = coord.level;
            int32_t x = coord.x;
            int32_t y = coord.y;
            int32_t z = coord.z;
            treeId tree_id = coord.tree;

            // check if cell is not refined
            if ( level == 0 )
            {
                return 0;
            }
            int quad_sibling_number = -1;
            if ( tdim == 2 )
            {
                // create associated quadrant 
                p4est_quadrant_t quad;
                quad.level = level;
                quad.x = x;
                quad.y = y;
                quad_sibling_number = pXest_get_sibling_nr ( &quad );
            }
            if ( tdim == 3 )
            {
                // create associated quadrant 
                p8est_quadrant_t quad;
                quad.level = level;
                quad.x = x;
                quad.y = y;
                quad.z = z;
                quad_sibling_number = pXest_get_sibling_nr ( &quad );
            }

            int sub_cell_number = pXest_z_order_to_ccw_order ( quad_sibling_number ) + 1;
            LOG_DEBUG ( 3, "id " << cell_id << " level " << level << " x " << x << " y " << y << " z " << z
                        << " tree " << tree_id << "Quad sibling number: " << quad_sibling_number << " sub cell number: " << sub_cell_number );
            return sub_cell_number;
#else
            return 0;
#endif
        }

        void MeshPXest::update_sub_cell_numbers ( )
        {
            //bool has_attr = this->has_attribute("__sub_cell_number__", this->tdim()); 
            int num_cells = this->num_entities ( this->tdim ( ) );
            std::vector<int> sub_cell_numbers ( num_cells );

            // Loop over all cell entities in mesh
            for ( EntityNumber index = 0; index < num_cells; ++index )
            {
                Id id = this->get_id ( this->tdim ( ), index );
                sub_cell_numbers[index] = this->get_sub_cell_number ( id );
            }
            this->sub_cell_number_attr_ = AttributePtr ( new IntAttribute ( sub_cell_numbers ) );
            this->add_attribute ( std::string ( "__sub_cell_number__" ), this->tdim ( ), sub_cell_number_attr_ );
        }

        inline EntityNumber MeshPXest::get_cell_index ( Id cell_id ) const
        {
            assert ( cell_id >= 0 );

            std::map<Id, EntityNumber>::const_iterator it;
            it = this->cell_id_to_index_.find ( cell_id );
            if ( it != this->cell_id_to_index_.end ( ) )
            {
                return it->second;
            }
            else
            {
                return -1;
            }
        }

        void MeshPXest::add_cells_and_rebuild ( SortedArray<Id>& new_cells )
        {
            bool new_cells_inserted = false;
            SortedArray<Id> cells = this->get_entities ( this->tdim ( ) );

            LOG_DEBUG ( 3, "Mesh " << this->history_index_ << ": old cells " << string_from_range ( cells.begin ( ), cells.end ( ) )
                        << ", added cells: " << string_from_range ( new_cells.begin ( ), new_cells.end ( ) ) );

            // Loop over cells to be added
            for ( int l = 0; l < new_cells.size ( ); ++l )
            {
                Id cell_id = new_cells[l];
                bool found = cells.find_insert ( cell_id );
                if ( !found )
                {
                    new_cells_inserted = true;
                }
            }

            // rebuild mesh
            if ( new_cells_inserted )
            {
                this->rebuild ( cells.data ( ) );
            }

        }

        void MeshPXest::rebuild ( std::vector<Id>& cells )
        {
            assert ( cells.size ( ) > 0 );

            LOG_DEBUG ( 2, "Rebuild mesh " << this->history_index_ << " with cells " << string_from_range ( cells.begin ( ), cells.end ( ) ) );
            int tdim = this->tdim ( );
            this->attributes_.reset ( new AttributeTable[tdim + 1] );
            this->entities_.reset ( new SortedArray<Id>[tdim + 1] );

            // PERF_TODO check for memory leak
            this->incidence_connections_ = new Connectivity[( tdim + 1 )*( tdim + 1 )];

            this->entities_[tdim] = SortedArray<Id>( cells );
            this->compute_entity_vertex_connectivity ( tdim );

            this->update_cell_id_to_index_map ( );
            this->update_sub_cell_numbers ( );

            // Push back mesh history for cells
            std::vector<int> history_indices ( cells.size ( ), this->history_index_ );
            this->pXest_db_->add_cell_to_mesh_maps ( cells, history_indices );
        }

        // ---------------------------------------------------------------------
        // BUG_TODO: something goes wrong for recursive refinement

        std::vector<Id> MeshPXest::compute_refined_cells_recursive ( const Entity& cell,
                                                                     RefinementTree* tree,
                                                                     std::vector<int>& sub_cell_numbers,
                                                                     std::vector<int>& tree_node_numbers,
                                                                     MeshPXestBuilder* builder,
                                                                     bool leaf_only ) const
        {
            const GDim gd = this->gdim ( );
            const TDim td = this->tdim ( );

            // call mesh::compute_refined_cells to get coordinates and cell sizes
            std::vector<Coordinate> refined_cell_coordinates;
            std::vector<int> refined_cell_sizes;

            ::hiflow::mesh::compute_refined_cells_recursive ( cell,
                                                              *tree,
                                                              refined_cell_coordinates,
                                                              refined_cell_sizes,
                                                              sub_cell_numbers,
                                                              tree_node_numbers,
                                                              leaf_only,
                                                              this->rec_ref_geom_fun_ );

            LOG_DEBUG ( 2, " Parent cell id " << cell.id ( ) );
            LOG_DEBUG ( 2, " Num descendants " << refined_cell_sizes.size ( ) << " coord size " << refined_cell_coordinates.size ( ) );
            LOG_DEBUG ( 2, " Tree node numbers " << string_from_range ( tree_node_numbers.begin ( ), tree_node_numbers.end ( ) ) );
            LOG_DEBUG ( 3, " Refined cell coords size " << string_from_range ( refined_cell_coordinates.begin ( ), refined_cell_coordinates.end ( ) ) );

            // add refined vertices to database
            const EntityCount num_refined_vertices = refined_cell_coordinates.size ( ) / gd;
            std::vector<Id> refined_vertex_ids ( num_refined_vertices );
            for ( int v = 0; v < num_refined_vertices; ++v )
            {
                if ( builder == NULL )
                {
                    refined_vertex_ids[v] = db_->add_vertex ( vec2ptr ( refined_cell_coordinates ) + gd * v );
                }
                else
                {
                    //    refined_vertex_ids[v] = builder->add_vertex ( vec2ptr ( refined_cell_coordinates ) + gd * v );
                }

                LOG_DEBUG ( 2, "added vertex " << v << " with id " << refined_vertex_ids[v]
                            << " and coords " << refined_cell_coordinates[gd * v] << ", " << refined_cell_coordinates[gd * v + 1] << " to DB " );
            }

            // add refined cells to database
            const EntityCount num_children_cells = refined_cell_sizes.size ( );
            assert ( num_children_cells > 1 );

            // get material numbers on cell facets
            std::vector<MaterialNumber> parent_facet_materials;

            if ( td > 1 )
            {
                parent_facet_materials = get_cell_facet_materials ( db_, cell );
            }

            // vector for new cell id:s
            std::vector<Id> new_cells ( num_children_cells );

            LOG_DEBUG ( 2, "refined_vertex_ids size = " << refined_vertex_ids.size ( ) << " and content " <<
                        string_from_range ( refined_vertex_ids.begin ( ), refined_vertex_ids.end ( ) ) );

            int offset = 0;
            for ( int rc = 0; rc < num_children_cells; ++rc )
            {
                // get vertex id:s
                const int cell_size = refined_cell_sizes[rc];
                const std::vector<Id> refined_cell_vertices ( refined_vertex_ids.begin ( ) + offset,
                                                              refined_vertex_ids.begin ( ) + offset + cell_size );
                offset += cell_size;

                LOG_DEBUG ( 2, "Refined cell " << rc << " with size " << cell_size << " and offset " << offset - cell_size << " = "
                            << string_from_range ( refined_cell_vertices.begin ( ), refined_cell_vertices.end ( ) ) );

                // add to database and return vector
                Id cell_id;
                if ( builder == NULL )
                {
                    cell_id = db_->add_entity ( td, refined_cell_vertices );
                }
                else
                {
                    //    cell_id = builder->add_entity ( td, refined_cell_vertices );
                }
                new_cells[rc] = cell_id;

                // inherit cell material number from parent
                if ( builder == NULL )
                {
                    db_->set_material_number ( td, cell_id, get_material_number ( td, cell.index ( ) ) );
                }
                else
                {
                    //    builder->set_material_number ( td, cell_id, get_material_number ( td, cell.index ( ) ) );
                }
                // TODO builder alternative
                // inherit facet material numbers from parent facets (see helper function above)
                if ( !parent_facet_materials.empty ( ) )
                {
                    inherit_facet_material_numbers ( db_,
                                                     cell.cell_type ( ),
                                                     parent_facet_materials,
                                                     sub_cell_numbers[rc],
                                                     refined_cell_vertices,
                                                     tree );
                }
            }
            return new_cells;
        }

        ////////////////////////////////////////////////////////////////////////
        //////////////// MeshPXestBuilder //////////////////////////////////////
        typedef MeshPXestBuilder::VertexHandle VertexHandle;
        typedef MeshPXestBuilder::EntityHandle EntityHandle;

        MeshPXestBuilder::MeshPXestBuilder ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        : MeshDbViewBuilder ( tdim, gdim, 0, period ),
        pXest_db_ ( new MeshPXestDatabase ( tdim, gdim ) )
        {
            db_ = pXest_db_;
            assert ( db_ != 0 );
            assert ( pXest_db_ != 0 );
        }

        MeshPXestBuilder::MeshPXestBuilder ( int tdim, int gdim,
                                             const MeshPXestDatabasePtr db,
                                             const std::vector<EntityHandle>& cells,
                                             std::vector<MasterSlave> period )
        : MeshDbViewBuilder ( tdim, gdim, 0, period )
        {
            this->pXest_db_ = db;
            this->db_ = db;
            this->cells_.resize ( cells.size ( ), 0 );
            this->cells_ = cells;
        }

        void MeshPXestBuilder::add_entity_to_quad_coord ( int tdim, Id entity_id, QuadCoord& coord, int which_mapping )
        {
            this->pXest_db_->add_entity_to_quad_coord ( tdim, entity_id, coord, which_mapping );
        }

        QuadCoord MeshPXestBuilder::get_entity_to_quad_coord ( int tdim, Id entity_id, int which_mapping )
        {
            return this->pXest_db_->get_entity_to_quad_coord ( tdim, entity_id, which_mapping );
        }

        void MeshPXestBuilder::set_conn_data ( int num_vertices, int num_cells, std::vector<double>& vertices, std::vector<int>& tree_to_vertices )
        {
            return this->pXest_db_->set_conn_data ( num_vertices, num_cells, vertices, tree_to_vertices );
        }

        MeshPtr MeshPXestBuilder::build ( )
        {
            if ( !cells_.empty ( ) )
            {
                std::vector<bool> locality;
                return MeshPtr ( new MeshPXest ( this->db_->tdim ( ), this->db_->gdim ( ), this->pXest_db_, this->cells_, 0, locality, this->period_ ) );
            }
            else
            {
                return MeshPtr ( 0 );
            }
        }

    }
} // namespace hiflow
