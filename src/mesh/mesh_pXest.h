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

/// \author Fabian Kissler, Jonathan Schwegler, Philipp Gerstner
// BUGS:
//   recursive refinement with more than 1 levels
//   copy of sub_cell_number attribute

#ifndef HIFLOW_MESH_PXEST_H
#    define HIFLOW_MESH_PXEST_H

#    include <map>
#    include <vector>
#    include <string>

#    include "common/sorted_array.h"
#    include "mesh_db_view.h"
#    include "mesh_database_pXest.h"
#    include "mesh_tools.h"

namespace hiflow
{
    namespace mesh
    {
        class Connectivity;
        class MeshPXestBuilder;

        ///
        /// \brief A view of some of the entities of a MeshDatabase.
        /// This class implements the Mesh interface.
        /// \see Mesh.

        class MeshPXest : public MeshDbView
        {
          public:
            // STRUCT_TODO same constructor interface for MeshPXest and MeshDbView
            /// \brief Constructor
            /// @param[in] tdim topological dimension of mesh
            /// @param[in] gdim geometric dimension of mesh
            /// @param[in] db pointer to underlying database
            /// @param[in] cells array of cell ids that make up the mesh
            /// @param[in] history_index index of mesh in array mesh_history of underlying database db
            /// @param[in] cell_is_local flag indicating for each cell whether it is locally owned or contained in ghost layer. can be empty array
            /// @param[in] period MasterSlave struct describing periodic boundaries
            MeshPXest ( TDim tdim, GDim gdim,
                        MeshPXestDatabasePtr db,
                        const std::vector<Id>& cells,
                        int history_index,
                        const std::vector<bool>& cell_is_local,
                        std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Constructor
            /// @param[in] tdim topological dimension of mesh
            /// @param[in] gdim geometric dimension of mesh
            MeshPXest ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Destructor

            ~MeshPXest ( )
            {
            };

            // STRUCT_TODO -> AdaptionP4est
            /// \brief create a new mesh by refining / coarsening this one. The output is always balanced. <br>
            /// Convention for refinement flags: <br>
            /// flags.size == cells.size <br>
            /// flags[i] = 0: keep cell i <br>
            /// flags[i] = 1: refine cell i <br>
            /// flags[i] < 0: refine family F = {c_1, .. c_n} of cells iff \sum{i = 1}^{n} flags[c_i] <= - n
            /// @param[in] refinements refinement indicator for each cell in current mesh (>= 0: refine one level, <0: do nothing)
            /// @return pointer to newly created refined mesh
            virtual MeshPtr refine ( std::vector<int>& flags ) const;

            /// \brief Get parent cell id of a cell
            /// @param[in] cell_index index of child cell
            /// @return Id of parent cell
            virtual Id get_parent_cell_id ( EntityNumber cell_index ) const;

            /// \brief return entity of parent cell of cell with givewn index
            /// @param[in] cell_index index of child cell
            /// @return entity of parent cell
            virtual Entity get_parent_cell ( EntityNumber cell_index ) const;

            /// \brief check whether cell with given index has a parent
            /// @param[in] cell_index index of child cell
            /// @return true if cell has parent in mesh
            virtual bool cell_has_parent ( EntityNumber cell_index ) const;

            /// \brief Obtain the id:s of the children of a cell entity.
            /// @param[in] cell_index index of parent cell
            /// @return vector of children ids
            virtual std::vector<Id> get_children_cell_ids ( EntityNumber cell_index ) const;

            /// \brief compute sibling cells of a given cell
            /// @param[in] cell_index index of one sibling
            /// @return cell indices of all cells in same family as input cell
            virtual std::vector<EntityNumber> get_sibling_cell_indices ( EntityNumber cell_index ) const;

            /// \brief set pointer to underlying database
            /// @param[in] db pointer to MeshPXestDatabase

            void set_db ( MeshPXestDatabasePtr db )
            {
                this->db_ = db;
                this->pXest_db_ = db;
            }

            /// \brief return pointer to underlying database
            /// @return pointer to MeshPXestDatabase

            MeshPXestDatabasePtr get_db ( ) const
            {
                return boost::static_pointer_cast<MeshPXestDatabase> ( this->db_ );
            }

            // STRUCT_TODO -> AdaptionP4est
            /// \brief refine sequential master mesh uniformly until a minimum number of cells is contained
            /// @param[in] num_ref_steps number of times, the mesh should be refined uniformly
            /// @return pointer to refined mesh
            MeshPtr refine_uniform_seq ( int num_ref_steps ) const;

            /// \brief refine mesh uniformly once
            /// @return pointer to refined mesh
            MeshPtr refine_uniform_seq ( ) const;

            /// \brief compute number of cells in ghost layer of mesh
            /// @return number of cells
            virtual int num_ghost_cells ( ) const;

            /// \brief returns information if a given cell belongs to the locally owned subdomain (i.e. no ghost cell)
            /// @param[in] cell_index index of some cell
            /// @return true if cell is locally owned
            virtual bool cell_is_local ( EntityNumber cell_index ) const;

            virtual void save ( std::string filename, const MPI_Comm& comm ) const;

            virtual void load ( std::string filename, const MPI_Comm& comm );

            /// \brief copy data from reference mesh. Here, the underlying databse is not copied, only its pointer
            /// @param[in] pointer to mesh which should be copied
            virtual void copy_from ( const MeshPtr mesh );

            /// \brief copy data from reference mesh, inlcuding the underlying database
            /// @param[in] pointer to mesh which should be copied
            virtual void deep_copy_from ( const MeshPtr mesh );

            /// \brief set flag for type of connectivity in balance mode. Balance ensures that two cell connected in the specified way differ in at most one level of refinement.
            /// @param[in] mode (2D: 0: Edge connectivity, 1: corner connectivity; 3D: 0: Face connectivity, 1: edge connectivity, 2: corner connectivity)
            virtual void set_connection_mode ( int mode );

            /// \brief set flag for path mode. If set to true, then additional cells are selected for refinement in order to obtain a mesh which is
            /// uniformly coarsenable
            /// @param[in] flag
            virtual void set_patch_mode ( bool flag );

            /// \brief checks if mesh is in such a state, that each family of cells is coarsenable.
            /// @return true if mesh is uniformyl coarsenable
            virtual bool is_uniformly_coarsenable ( ) const;

            virtual int get_history_index ( ) const
            {
                return history_index_;
            }
          protected:
            // STRUCT_TODO -> AdaptionP4est
            /// \brief refines the mesh in order to ensure haviong at most one hanging node per facet.
            /// This is done by calling the p4est_balance_ext routine.
            /// @param[out] add_cells new cells that were created during the balanciing process and which should be added to the mesh
            /// @param[out] rm_cells cells that were refined during the balancing process and which should be removed from the mesh
            /// @param[in] connection_mode type of connectivity
            void balance ( std::vector<Id>& add_cells, std::vector<Id>& rm_cells, int connection_mode ) const;

            /// \brief finds all cells belonging t a family which is ready for coarsening
            /// @param[out] indicex of cells that are not coarsenable
            void find_coarsenable_cells ( SortedArray<EntityNumber>& coarsenable_cells ) const;

            /// \brief yields a mesh which is uniformly coarsenable by refining non coarsenable cells
            /// @return pointer to newly created, refined mesh
            MeshPtr make_mesh_coarsenable ( );

            /// \brief for each local quadrant, store index of corresponding cell in associated QuadData.
            /// Stores -1 in quadrants whose corresponding cell is not in current mesh
            /// @param[in] update_ghost flag indicating whether quadrants in ghost layer should be updated as well
            void update_cell_index_in_quads ( bool update_ghost );

            // STRUCT_TODO -> GhostP4est
            /// \brief add ghost cells and facets to existing mesh
            /// @param[in] comm MPI communicator
            /// @param[in] layer_width width of ghost layer
            /// @return pointer to newly created mesh containing ghost layer
            virtual MeshPtr add_ghost_cells ( MPI_Comm comm, int layer_width ) const;

            /// \brief Access the id of an ancestor of a cell entity
            /// @param[in] cell_id Id of the descendant
            /// @param[in] level_diff number of generations between cell and its ancestor (e.g. for parent: level_diff = 1)
            /// @return Id of ancestor cell
            Id get_ancestor_cell_id ( Id cell_id, int level_diff ) const;

            /// \brief Compute the cell ids of all descendants of a cell entity
            /// @param[in] cell_id Id of the ancestor
            /// @param[in] level_diff number of generations between cell and its descendants (e.g. for children: level_diff = 1)
            /// @return ids of descendant cells
            std::vector<Id> get_descendant_cell_ids ( Id cell_id, int level_diff ) const;

            /// \brief Return ids of current ghost cells
            /// @return ids of ghsot cells
            SortedArray<Id> get_ghost_ids ( ) const;

            // STRUCT_TODO -> GhostP4est
            /// \brief initialze p4est ghost layer -> build p4est_ghost_t
            /// @param[in] layer_width width of ghost layer
            void init_ghost_layer ( int layer_width ) const;

            // STRUCT_TODO -> GhostP4est
            /// \brief Put data that should be exchanged into mirror quads (local quads which are in ghost layer of other processes)
            void init_ghost_data ( ) const;

            // STRUCT_TODO -> GhostP4est
            /// \brief exchange data in ghost layers and put it into builder, sub_domains and remote_indices
            /// @param[out] ids of all cells in local_mesh + ghost_layer
            void exchange_ghost_data ( std::vector< SortedArray<Id> >& new_cells ) const;

            // STRUCT_TODO -> GhostP4est
            /// \brief put cell indices into ghost data after exchanging ghost cells -> communicate indices of newly created ghost cells
            void init_ghost_indices ( ) const;

            /// \brief exchange indices of newly created ghost cells
            // STRUCT_TODO -> GhostP4est
            void exchange_ghost_indices ( ) const;

            // STRUCT_TODO -> GhostP4est
            /// \brief compute subdomain index and remote index for each cell in local mesh
            /// @param[out] sub_domain
            /// @param[out] remote_indices
            void compute_subdomains_and_remote_indices ( std::vector<SubDomainId>& sub_domains, std::vector<EntityNumber>& remote_indices ) const;

            /// \brief setup mapping Id -> Id for cells
            void update_cell_id_to_index_map ( );

            /// \brief Access data from cell_id_to_index
            /// @param[in] cell_id id of some cell
            /// @return index of given cell
            inline EntityNumber get_cell_index ( Id cell_id ) const;

            /// \brief add new cells to the existing mesh and rebuild if necessary
            /// @param[in] new_cells
            void add_cells_and_rebuild ( SortedArray<Id>& new_cells );

            /// \brief rebuild mesh by using existing database and inclusion of cells.
            /// This functions mimics the call of MeshDbView( MeshDbVievBuilder.build() )
            /// @param[in] cells id of cells making up the updated mesh
            void rebuild ( std::vector<Id>& cells );

            /// \brief return sub_cell_number of cell with given Id
            /// @param[in] cell_id
            /// @return sub cell number
            int get_sub_cell_number ( Id cell_id ) const;

            /// \brief updates the sub_cell_number attribute of the mesh
            void update_sub_cell_numbers ( );

            /// \brief set flags indicating for each cell, whether it is locally owned or not
            /// @param[in] flags
            void set_cell_flags ( std::vector<bool>& flags );

            // STRUCT_TODO -> AdaptionDataBase
            /// \brief recursively refine a given cell according to a given tree. This functions adds new cell entities to the
            /// underlying database for each node in the tree, i.e. also for each non-leaf node.
            /// @param[in] cell the cell entity to refine
            /// @param[in] tree tree describing the refinement.
            /// @param[out] sub_cell_numbers the sub cell numbers for each created cell
            /// @param[out] tree_node_numbers the corresponding index in the refined_cell array of the tree for each created cell
            /// @param[in] builder pointer to mesh builder
            /// @param[in] leaf_only if true, only leaf cells are returned
            /// @return Ids of refined cells
            std::vector<Id> compute_refined_cells_recursive ( const Entity& cell,
                                                              RefinementTree* tree,
                                                              std::vector<int>& sub_cell_numbers,
                                                              std::vector<int>& tree_node_numbers,
                                                              MeshPXestBuilder* builder,
                                                              bool leaf_only ) const;

            // STRUCT_TODO -> InterfaceP4est
            /// \brief creates an interface list by using the p4est iterator functionality
            /// @param[out] list
            void create_interface_list ( InterfaceList* list ) const;

            /// \brief mpi rank
            /// @return mpi rank
            int get_rank ( ) const;

            /// \brief pointer to underlying database
            MeshPXestDatabasePtr pXest_db_;

            /// \brief cell id to index mapping
            std::map<Id, EntityNumber> cell_id_to_index_;

            // STRUCT_TODO -> GhostP4est
            /// \brief pointer to temporary data that has to be communicated when exchanging ghost data
            QuadData* ghost_data_;

            // STRUCT_TODO -> put into MeshDbView
            /// \brief index of mesh in mesh_history array of underlying meshPXestDatabase
            int history_index_;

            // STRUCT_TODO -> put into MeshDbView
            /// \brief sub cell number attribute
            AttributePtr sub_cell_number_attr_;

            // STRUCT_TODO -> put into MeshDbView
            /// \brief flag indicating for each cell whether it is local or in the ghost layer
            std::vector<bool> cell_is_local_;

            /// \brief type of connection in balance routine
            int connection_mode_;

            /// \brief flag indicatin whether patch mode is used
            bool patch_mode_;

          private:
            // STRUCT_TODO -> GhostP4est
            mutable std::vector<GhostCommData22*> mirror_data_22_;
            mutable std::vector<GhostCommData33*> mirror_data_33_;
            mutable std::vector<GhostCommData*> mirror_data_;

            friend class MeshP4estBuilder;

            friend MeshPtr hiflow::partition_and_distribute_pXest ( const mesh::MeshPtr master_mesh,
                                                                    const int master_rank,
                                                                    const MPI_Comm& comm,
                                                                    int& uniform_ref_steps );

            friend MeshPtr hiflow::compute_ghost_cells_pXest ( const Mesh& local_mesh,
                                                               const MPI_Comm& comm,
                                                               int layer_width );

        };

        /// \brief A MeshBuilder that creates a MeshDbView
        /// \see MeshBuilder
        /// \see MeshDbView

        class MeshPXestBuilder : public MeshDbViewBuilder
        {
          public:
            typedef MeshBuilder::VertexHandle VertexHandle;
            typedef MeshBuilder::EntityHandle EntityHandle;

            /// \brief Constructor to use when building a Mesh with an existing MeshDatabase
            /// @param[in] db pointer to database
            explicit MeshPXestBuilder ( MeshDatabasePtr db, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Constructor to initialize builder with existing MeshDbView mesh
            /// @param[in] tdim topological dimension of mesh
            /// @param[in] gdim geometric dimension of mesh
            /// @param[in] db pointer to database
            /// @param[in] cells array of cell ids forming the new mesh
            explicit MeshPXestBuilder ( int tdim, int gdim,
                                        const MeshPXestDatabasePtr db,
                                        const std::vector<EntityHandle>& cells,
                                        std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Constructor to use when building a Mesh without an existing MeshDatabase
            /// @param[in] tdim topological dimension of mesh
            /// @param[in] gdim geometric dimension of mesh
            explicit MeshPXestBuilder ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Destructor

            ~MeshPXestBuilder ( )
            {
            }

#    ifdef WITH_P4EST

            /// \brief set p4est conenctivity object
            /// \param[in] conn pointer to connectivity

            void set_p4est_conn ( p4est_connectivity_t* conn )
            {
                pXest_db_->set_p4est_conn ( conn );
            }

            /// \brief set p8est conenctivity object
            /// \param[in] conn pointer to connectivity

            void set_p8est_conn ( p8est_connectivity_t* conn )
            {
                pXest_db_->set_p8est_conn ( conn );
            }

            /// \brief set p4est forest object
            /// \param[in] forest pointer to forest

            void set_p4est ( p4est_t* forest )
            {
                pXest_db_->set_p4est_forest ( forest );
            }

            /// \brief set p8est forest object
            /// \param[in] forest pointer to forest

            void set_p8est ( p8est_t* forest )
            {
                pXest_db_->set_p8est_forest ( forest );
            }
#    endif

            /// \brief add entity (meshdatabase) to quadrant (p4est) mapping
            /// @param[in] tdim top. dimension of entity
            /// @param[in] entity_id id of entity
            /// @param[in] coord struct that describes coordinates of quadrant in forest which corresponds to given entity
            /// @param[in] which_mapping 0: mapping for coarse mesh, 1: mapping for refined mesh
            void add_entity_to_quad_coord ( int tdim, Id entity_id, QuadCoord& coord, int which_mapping );

            /// \brief get some entity to quadrant mapping
            /// @param[in] tdim top. dimension of entity
            /// @param[in] entity_id id of entity
            /// @param[in] which_mapping 0: mapping for coarse mesh, 1: mapping for refined mesh
            /// @return Quadcoord object
            QuadCoord get_entity_to_quad_coord ( int tdim, Id entity_id, int which_mapping );

            /// \brief set data used to build p4est connectivity
            /// @param[in] num_vertices number of vertices in coarse mesh
            /// @param[in] num_cells number of cells in coarse mesh
            /// @param[in] vertices coordinates of vertices
            /// @param[in] tree_to_vertices array containing vertex ids belonging to each cell
            void set_conn_data ( int num_vertices, int num_cells, std::vector<double>& vertices, std::vector<int>& tree_to_vertices );

            /// \brief get ids of cells making up the new mesh
            /// @return ids

            std::vector<Id> get_cells ( )
            {
                return this->cells_;
            }

            /// \brief build the mesh
            /// @return pointer to new mesh
            virtual MeshPtr build ( );

          protected:

            /// \brief build p4est connectivity object
            //void build_connectivity();

            /// pointer to underlying database
            MeshPXestDatabasePtr pXest_db_;
        };
    }
} // namespace hiflow

#endif
