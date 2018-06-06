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

#include "numbering_lagrange.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <sstream>

#include <mpi.h>

#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/king_ordering.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "fe_interface_pattern.h"
#include "dof_interpolation_pattern.h"
#include "fem/felagrange.h"
#include "mesh/iterator.h"
#include "mesh/entity.h"
#include "mesh/types.h"
#include "mesh/cell_type.h"
#include "common/macros.h"
#include "fem/femanager.h"
#include "fem/cell_transformation.h"
#include "degree_of_freedom.h"
#include "common/log.h"
#include "common/sorted_array.h"
#include "space/vector_space.h"

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        const int DEBUG_LEVEL = 0;

        /// creates DoF numbering, i.e. calculates mapping of local to global DoF IDs
        /// and interpolation for hanging DoFs

        template<class DataType>
        void NumberingLagrange<DataType>::number ( DOF_ORDERING order )
        {

            const int verbatim_mode = 0; // 0 -> no output
            // 1 -> some output
            // initial numbering, some DoFs may have multiple DoF IDs
            initial_numbering ( );

            // find DoF IDs that lie on same coordinate and should be identified
            identify_common_dofs ( );

            // Determine equivalence classes and modify numer_ s.t.
            // each numer_ entry maps to one representer of the equivalence class
            //
            // The following algorithm is based on the description given in the book
            // "Numerical Recipes in C", Second edition,
            // W. Press, S. Teukolsky, W. Vetterling, B. P. Flannery, Cambridge University
            // Press, Pages 345 / 346

            size_t cntr = 0;

            while ( cntr < this->dof_interpolation_->dof_identification_list ( ).size ( ) )
            {
                int cntr_val = this->dof_interpolation_->dof_identification_list ( )[cntr].first;

                while ( ( *this->numer_ )[cntr_val] != cntr_val )
                    cntr_val = ( *this->numer_ )[cntr_val];

                int cntr_val2 = this->dof_interpolation_->dof_identification_list ( )[cntr].second;

                while ( ( *this->numer_ )[cntr_val2] != cntr_val2 )
                    cntr_val2 = ( *this->numer_ )[cntr_val2];

                if ( cntr_val != cntr_val2 )
                    ( *this->numer_ )[cntr_val] = cntr_val2;

                ++cntr;
            }

            const size_t e_cntr = this->numer_->size ( );
            for ( cntr = 0; cntr != e_cntr; ++cntr )
                while ( ( *this->numer_ )[cntr] != ( *this->numer_ )[( *this->numer_ )[cntr]] )
                    ( *this->numer_ )[cntr] = ( *this->numer_ )[( *this->numer_ )[cntr]];

            // Need equivalence reduction also in dof constraints
            this->dof_interpolation_->apply_permutation ( *( this->numer_ ) );

            // 2. Consecutive numbering

            std::vector<int> numer_copy ( *( this->numer_ ) );
            std::sort ( numer_copy.begin ( ), numer_copy.end ( ) );

            std::vector<int>::iterator it = std::unique ( numer_copy.begin ( ), numer_copy.end ( ) );
            size_t numer_copy_size = it - numer_copy.begin ( );

            std::vector<int> permutation ( this->numer_->size ( ), -5 );
            for ( size_t i = 0; i != numer_copy_size; ++i )
                permutation[numer_copy[i]] = i;

            LOG_DEBUG ( 2, "Numer size = " << this->numer_->size ( ) );
            LOG_DEBUG ( 2, "Numer = " << string_from_range ( this->numer_->begin ( ), this->numer_->end ( ) ) );
            LOG_DEBUG ( 2, "Numer copy size = " << numer_copy.size ( ) );
            LOG_DEBUG ( 2, "Numer copy = " << string_from_range ( numer_copy.begin ( ), numer_copy.end ( ) ) );
            LOG_DEBUG ( 2, "Permutation size = " << permutation.size ( ) );
            LOG_DEBUG ( 2, "Permutation = " << string_from_range ( permutation.begin ( ), permutation.end ( ) ) );

            this->apply_permutation ( permutation );

            if ( verbatim_mode > 0 )
            {
                std::cout << std::endl;
                std::cout << "INTERPOLATION INFORMATION" << std::endl;
                std::cout << "=========================" << std::endl;
                this->dof_interpolation_->backup ( std::cout );
                std::cout << std::endl;

                std::cout << std::endl;
                std::cout << "NUMER-FIELD" << std::endl;
                std::cout << "===========" << std::endl;
                for ( int i = 0; i < this->numer_->size ( ); ++i )
                    std::cout << i << "\t ->    " << this->numer_->at ( i ) << std::endl;
                std::cout << std::endl;

                if ( this->dof_interpolation_->size ( ) > 1 )
                {
                    std::cout << "# DoFs:              "
                            << this->number_of_dofs_total_ << std::endl;
                    std::cout << "# interpolated DoFs: "
                            << this->dof_interpolation_->size ( ) << std::endl;
                    std::cout << "# real DoFs:         "
                            << this->number_of_dofs_total_ - this->dof_interpolation_->size ( )
                            << std::endl;
                }
                else
                    std::cout << "# DoFs: " << this->number_of_dofs_total_ << std::endl;
            }

            // 3. Apply optimization to ordering according to order
            if ( order != HIFLOW_CLASSIC )
            {
                LOG_INFO ( "DoF reordering strategy", order );

                // Some Typedefs for clearer notation
                using namespace boost;
                using namespace std;
                typedef adjacency_list<vecS, vecS, undirectedS,
                        property<vertex_color_t, default_color_type,
                        property<vertex_degree_t, int> > > Graph;
                typedef graph_traits<Graph>::vertex_descriptor Vertex;
                typedef graph_traits<Graph>::vertices_size_type size_type;

                // Create Boost graph of local dof couplings
                Graph G ( this->number_of_dofs_total_ );
                {
                    // Create graph topology
                    std::vector<int> dof_ind_test, dof_ind_trial;
                    std::vector< SortedArray<int> > raw_struct ( this->number_of_dofs_total_ );

                    // loop over every cell (including ghost cells)
                    typename VectorSpace<DataType>::MeshEntityIterator mesh_it = this->mesh_->begin ( this->mesh_->gdim ( ) );
                    typename VectorSpace<DataType>::MeshEntityIterator e_mesh_it = this->mesh_->end ( this->mesh_->gdim ( ) );
                    while ( mesh_it != e_mesh_it )
                    {
                        // loop over test variables
                        for ( int test_var = 0, tv_e = this->nvar_; test_var != tv_e; ++test_var )
                        {
                            // get dof indices for test variable
                            this->dof_->get_dofs_on_cell ( test_var, ( *mesh_it ).index ( ), dof_ind_test );

                            // loop over trial variables
                            for ( int trial_var = 0, vt_e = this->nvar_; trial_var != vt_e; ++trial_var )
                            {

                                // get dof indices for trial variable
                                this->dof_->get_dofs_on_cell ( trial_var, ( *mesh_it ).index ( ), dof_ind_trial );

                                // detect couplings
                                for ( size_t i = 0, i_e = dof_ind_test.size ( ); i != i_e; ++i )
                                {
                                    for ( size_t j = 0, j_e = dof_ind_trial.size ( ); j != j_e; ++j )
                                    {
                                        raw_struct[dof_ind_test[i]].find_insert ( dof_ind_trial[j] );
                                    } // for (int j=0;...
                                } // for (int i=0;...
                            }
                        }
                        // next cell
                        ++mesh_it;
                    } // while (mesh_it != ...

                    for ( size_t k = 0, k_e = raw_struct.size ( ); k != k_e; ++k )
                    {
                        for ( size_t l = 0, l_e = raw_struct[k].size ( ); l != l_e; ++l )
                        {
                            add_edge ( k, raw_struct[k][l], G );
                        }
                    }
                }

                graph_traits<Graph>::vertex_iterator ui, ui_end;

                property_map<Graph, vertex_degree_t>::type deg = get ( vertex_degree, G );
                for ( boost::tie ( ui, ui_end ) = vertices ( G ); ui != ui_end; ++ui )
                    deg[*ui] = degree ( *ui, G );

                property_map<Graph, vertex_index_t>::type
                index_map = get ( vertex_index, G );

                LOG_INFO ( "Original bandwidth", bandwidth ( G ) );

                std::vector<Vertex> inv_perm ( num_vertices ( G ) );
                std::vector<int> perm ( num_vertices ( G ) );

                switch ( order )
                {
                    case HIFLOW_CLASSIC:
                    {
                        break;
                    }

                    case CUTHILL_MCKEE:
                    {
                        //reverse cuthill_mckee_ordering
                        cuthill_mckee_ordering ( G, inv_perm.rbegin ( ), get ( vertex_color, G ),
                                                 make_degree_map ( G ) );
                        break;
                    }

                    case KING:
                    {
                        //king_ordering
                        king_ordering ( G, inv_perm.rbegin ( ) );
                        break;
                    }

                    default:
                    {
                        break;
                    }
                }

                for ( size_type c = 0; c != inv_perm.size ( ); ++c )
                    perm[index_map[inv_perm[c]]] = static_cast < int > ( c );
                LOG_INFO ( "Bandwidth after reordering",
                           bandwidth ( G, make_iterator_property_map ( &perm[0], index_map, perm[0] ) ) );

                this->apply_permutation ( perm );
            }

            if ( verbatim_mode > 0 )
            {
                std::cout << std::endl;
                std::cout << "INTERPOLATION INFORMATION" << std::endl;
                std::cout << "=========================" << std::endl;
                this->dof_interpolation_->backup ( std::cout );
                std::cout << std::endl;

                std::cout << std::endl;
                std::cout << "NUMER-FIELD" << std::endl;
                std::cout << "===========" << std::endl;
                for ( int i = 0; i < this->numer_->size ( ); ++i )
                    std::cout << i << "\t ->    " << this->numer_->at ( i ) << std::endl;
                std::cout << std::endl;

                if ( this->dof_interpolation_->size ( ) > 1 )
                {
                    std::cout << "# DoFs:              "
                            << this->number_of_dofs_total_ << std::endl;
                    std::cout << "# interpolated DoFs: "
                            << this->dof_interpolation_->size ( ) << std::endl;
                    std::cout << "# real DoFs:         "
                            << this->number_of_dofs_total_ - this->dof_interpolation_->size ( )
                            << std::endl;
                }
                else
                    std::cout << "# DoFs: " << this->number_of_dofs_total_ << std::endl;
            }

        }

        /// create initial DoF numbering ignoring that some DoFs coincide

        template<class DataType>
        void NumberingLagrange<DataType>::initial_numbering ( )
        {
            this->numer_offsets_cell_varloc_->clear ( );
            int numer_size = 0; // size of numer_ field: size = sum_var sum_cells n_local_dofs(var,cell)
            this->numer_offsets_cell_varloc_->resize ( this->nvar_ );

            // loop over vars

            for ( int var = 0; var < this->nvar_; ++var )
            {
                // loop over cells
                this->numer_offsets_cell_varloc_->at ( var ).reserve ( this->mesh_->num_entities ( this->fe_manager_->get_dim ( ) ) );
                for ( mesh::EntityIterator it = this->mesh_->begin ( this->fe_manager_->get_dim ( ) ),
                      e_it = this->mesh_->end ( this->fe_manager_->get_dim ( ) );
                      it != e_it;
                      ++it )
                {
                    const int n_local_dofs = this->fe_manager_->get_fe_on_cell ( it->index ( ), var )->get_nb_dof_on_cell ( );
                    ( *this->numer_offsets_cell_varloc_ )[var].push_back ( numer_size );
                    numer_size += n_local_dofs;
                }

            }

            // initialize numer_

            this->numer_->resize ( numer_size );

            for ( size_t i = 0; i != numer_size; ++i )
            {
                ( *this->numer_ )[i] = i;
            }
        }

        /// common dofs are identified and unified

        template<class DataType>
        void NumberingLagrange<DataType>::identify_common_dofs ( )
        {

            // prepare interface list

            mesh::InterfaceList interface_list = mesh::InterfaceList::create ( this->mesh_ );

            if ( DEBUG_LEVEL > 0 )
            {
                std::cout << interface_list << std::endl;
            }

            // clear old interpolation and identification information

            this->dof_interpolation_->clear_entries ( );

            // loop over interfaces
            for ( mesh::InterfaceList::const_iterator interface = interface_list.begin ( ),
                  e_interface = interface_list.end ( );
                  interface != e_interface;
                  ++interface )
            {
                // loop over vars

                for ( size_t var = 0; var != this->nvar_; ++var )
                {
                    // only identify continuous vars

                    if ( this->fe_manager_->get_ca ( var ) == true )
                    {
                        // calculate interface pattern

                        // FE ansatz of master cell
                        FEType<DataType>* master_type = this->fe_manager_->get_fe_on_cell ( interface->master_index ( ), var );

                        // FE ansatz of slave cells
                        std::vector<FEType<DataType>* > slave_types;
                        slave_types.reserve ( interface->num_slaves ( ) );
                        for ( size_t i = 0, e_i = interface->num_slaves ( ); i != e_i; ++i )
                        {
                            slave_types.push_back ( this->fe_manager_->get_fe_on_cell ( interface->slave_index ( i ), var ) );
                        }

                        // definition of the interface pattern

                        mesh::InterfacePattern pre_pattern = mesh::compute_interface_pattern ( this->mesh_, *interface );
                        FEInterfacePattern<DataType> pattern ( pre_pattern, master_type, slave_types );

                        // only treat interfaces between at least two cells (-> no boundaries)
                        if ( pre_pattern.num_slaves ( ) > 0 )
                        {
                            // add new interpolation definition for new interface patterns if needed
                            typename std::map<FEInterfacePattern<DataType>, DofInterpolationPattern>::const_iterator interpolation_it;
                            interpolation_it = interface_patterns_interpolation_.find ( pattern );

                            if ( interpolation_it == interface_patterns_interpolation_.end ( ) )
                            {
                                DofInterpolationPattern interpolation;
                                compute_interpolation ( pattern, *interface, interpolation );
                                interface_patterns_interpolation_.insert ( std::make_pair ( pattern, interpolation ) );
                                interpolation_it = interface_patterns_interpolation_.find ( pattern );
                            }
                            DofInterpolationPattern const& interpolation = interpolation_it->second;

                            // apply interpolation rules to DoFs

                            // -> get DoF IDs of master and slave cells

                            std::vector<DofID> dof_master;
                            this->dof_->get_dofs_on_cell ( var, interface->master_index ( ), dof_master );

                            std::vector<std::vector<DofID> > dof_slaves;
                            dof_slaves.resize ( pattern.num_slaves ( ) );
                            for ( size_t s = 0, e_s = pattern.num_slaves ( ); s != e_s; ++s )
                                this->dof_->get_dofs_on_cell ( var, interface->slave_index ( s ), dof_slaves[s] );

                            // -> identification of DoFs (Slave <-> Master)

                            for ( size_t s = 0, e_s = pattern.num_slaves ( ); s != e_s; ++s )
                            {
                                std::vector<std::pair<DofID, DofID> > const& identification_list =
                                        interpolation.interpolation_slave ( s ).dof_identification_list ( );
                                for ( size_t i = 0, e_i = identification_list.size ( ); i != e_i; ++i )
                                {
                                    DofID dof_id_slave = dof_slaves[s][identification_list[i].first];
                                    DofID dof_id_master = dof_master[identification_list[i].second];
                                    this->dof_interpolation_->insert_dof_identification ( dof_id_slave, dof_id_master );
                                }
                            }

                            // TODO: refactor master-master and slave-master translation into common function.
                            // -> interpolation of DoFs (Master <-> Master)

                            DofInterpolation const& int_master = interpolation.interpolation_master ( );
                            for ( DofInterpolation::const_iterator it = int_master.begin ( ), e_it = int_master.end ( ); it != e_it; ++it )
                            {
                                // get DoF ID of interpolated DoF
                                DofID dof_id_master = dof_master[it->first];

                                // get DoF IDs of interpolating DoFs
                                std::vector<std::pair<DofID, double> > sum = it->second;
                                for ( size_t md = 0, e_md = sum.size ( ); md != e_md; ++md )
                                    sum[md].first = dof_master[sum[md].first];

                                // add interpolation definition to DofInterpolation
                                bool status;
                                status = this->dof_interpolation_->push_back ( make_pair ( dof_id_master, sum ) );
                                // TODO (Staffan): Is the following assert correct? (see below)

                                assert ( status );
                            }

                            // -> interpolation of DoFs (Slave <-> Master)

                            for ( size_t s = 0, e_s = pattern.num_slaves ( ); s != e_s; ++s )
                            {
                                DofInterpolation const& int_slave = interpolation.interpolation_slave ( s );
                                for ( DofInterpolation::const_iterator it = int_slave.begin ( ), e_it = int_slave.end ( ); it != e_it; ++it )
                                {
                                    // get DoF ID of interpolated slave DoF
                                    DofID dof_id_slave = dof_slaves[s][it->first];

                                    // get DoF IDs of interpolating master DoFs
                                    std::vector<std::pair<DofID, double> > sum = it->second;
                                    for ( size_t md = 0, e_md = sum.size ( ); md != e_md; ++md )
                                        sum[md].first = dof_master[sum[md].first];

                                    // add interpolation definition to DofInterpolation
                                    bool status;
                                    status = this->dof_interpolation_->push_back ( make_pair ( dof_id_slave, sum ) );
                                    // Note (staffan): This assert is not correct, at least
                                    // in 3d, since constraints for dofs on hanging edges
                                    // will be created at least twice.

                                    // assert ( status );
                                }
                            } // for (int s=0; s<pattern.num_slaves(); ++s)
                        } // if (pre_pattern.num_slaves() >= 1)
                    } // if (fe_manager_->get_ca(var) == true)
                } // for (int var=0; var<fe_manager_->get_nb_var(); ++var)
            } // for (interface = interface_list.begin();...)

        }

        template<class DataType>
        void NumberingLagrange<DataType>::print_interface_patterns ( ) const
        {
            typename std::map<FEInterfacePattern<DataType>, DofInterpolationPattern>::const_iterator it;

            std::cout << "Interface-Modes:" << std::endl;
            for ( it = interface_patterns_interpolation_.begin ( );
                  it != interface_patterns_interpolation_.end ( );
                  ++it )
            {
                std::cout << it->first << std::endl;
                std::cout << it->second << std::endl;
            }
        }

        namespace
        {
            // Helper functions for compute_interpolation().

            template<class DataType>
            void log_debug_fe_pattern_info ( const mesh::Entity& entity, const FEType<DataType>& fe_type, int gdim )
            {
                LOG_DEBUG ( 2, "Cell index = " << entity.index ( ) << ", fe degree = " << fe_type.get_fe_deg ( ) );
                std::vector<DataType> coord;
                entity.get_coordinates ( coord );
                for ( int v = 0; v < entity.num_vertices ( ); ++v )
                {
                    LOG_DEBUG ( 2, "Vertex " << v << " = "
                                << string_from_pointer_range ( &coord[gdim * v], &coord[gdim * ( v + 1 )] ) );
                }
            }

            template<class DataType>
            void log_compute_interpolation_start ( const FEInterfacePattern<DataType>& pattern,
                                                   const mesh::Interface& interface,
                                                   const mesh::Mesh* mesh )
            {
                LOG_DEBUG ( 1, "\nCOLLECTING INFORMATION FOR THE INTERPOLATION\n"
                            << "============================================\n"
                            << "FEInterfacePattern = \n" << pattern << "\n" );

                const mesh::Entity& master = mesh->get_entity ( mesh->tdim ( ), interface.master_index ( ) );

                LOG_DEBUG ( 2, "Master Info" );
                log_debug_fe_pattern_info ( master, pattern.fe_type_master ( ), mesh->gdim ( ) );

                const std::vector<FEType<DataType>* >& fe_type_slaves = pattern.fe_type_slaves ( );

                for ( int s = 0; s < pattern.num_slaves ( ); ++s )
                {
                    const mesh::Entity& slave = mesh->get_entity ( mesh->tdim ( ), interface.slave_index ( s ) );

                    LOG_DEBUG ( 2, "Slave " << s << " Info" );
                    log_debug_fe_pattern_info ( slave, *fe_type_slaves[s], mesh->gdim ( ) );
                }
            }

            template<class DataType>
            void log_interface_matrix ( int level, const typename NumberingLagrange<DataType>::InterfaceMatrix& interface_matrix )
            {
                std::stringstream sstr;
                sstr << std::setprecision ( 4 );

                for ( typename NumberingLagrange<DataType>::InterfaceMatrix::const_iterator it = interface_matrix.begin ( ), end = interface_matrix.end ( );
                      it != end; ++it )
                {
                    sstr << "\t" << it->first << "\t->";

                    for ( int i = 0; i < it->second.size ( ); ++i )
                    {
                        sstr << "\t" << it->second[i];
                    }

                    const double sum = std::accumulate ( it->second.begin ( ), it->second.end ( ), 0. );

                    sstr << "\tRow sum = " << sum << "\n";
                }
                LOG_DEBUG ( level, "\n" << sstr.str ( ) );
            }

            template<class DataType>
            void extract_interface_matrix ( const typename NumberingLagrange<DataType>::InterpolationWeights& interpolation_weights,
                                            const std::vector<DofID>& interface_dofs,
                                            typename NumberingLagrange<DataType>::InterfaceMatrix& interface_matrix )
            {
                // Copy rows corresponding to interface dofs.
                for ( size_t i = 0, end = interface_dofs.size ( ); i != end; ++i )
                {
                    const DofID dof = interface_dofs[i];
                    interface_matrix.insert ( std::make_pair ( dof, interpolation_weights[dof] ) );
                }

#if 0
                // Filter small entries (not needed if interpolation_weights are already filtered?)
                for ( InterpolationMap::iterator it = interface_matrix.begin ( ), end = interface_matrix.end ( ); it != end; ++it )
                {
                    for ( size_t i = 0, e_i = it->second.size ( ); i != e_i; ++i )
                    {
                        if ( std::abs ( it->second[i] ) < COEFFICIENT_EPS )
                        {
                            it->second[i] = 0.0;
                        }
                    }
                }
#endif
            }

            template<class DataType>
            void multiply_interface_matrices ( const typename NumberingLagrange<DataType>::InterfaceMatrix& A,
                                               const typename NumberingLagrange<DataType>::InterfaceMatrix& B,
                                               typename NumberingLagrange<DataType>::InterfaceMatrix& C )
            {

                // It is assumed that the sizes of the matrices are as follows: A
                // = M x N; B = P x Q; where M is the number of interface dofs of
                // the constrained element, N is the number of cell dofs of the
                // constraining element, P is the number of interface dofs of the
                // constraining element, and Q is the number of cell dofs of the
                // constrained element. This means that the keys of B all lie in
                // the range 0..N-1, and the keys of A lie in the range 0..Q-1 .

                // Copy A into C, since they will have the same row indices.
                C = A;

                // Determine num columns of C = num columns of B (same for all rows).
                const int num_cols = B.begin ( )->second.size ( );

                // Resize all columns of C.
                for ( typename NumberingLagrange<DataType>::InterfaceMatrix::iterator it = C.begin ( ), end = C.end ( ); it != end; ++it )
                {
                    it->second.clear ( );
                    it->second.resize ( num_cols, 0. );
                }

                // Matrix multiplication: triple loop. In index notation: C_ij =
                // \sum{k}{A_ik * B_kj} . Here, i = rowC->first, k = rowB->first
                // and j loops through columns of B.
                // TODO: not clear if this is an efficient matrix-matrix multiplication.
                for ( typename NumberingLagrange<DataType>::InterfaceMatrix::iterator rowC = C.begin ( ), endRowC = C.end ( ); rowC != endRowC; ++rowC )
                {
                    for ( typename NumberingLagrange<DataType>::InterfaceMatrix::const_iterator rowB = B.begin ( ), endRowB = B.end ( ); rowB != endRowB; ++rowB )
                    {
                        const DataType A_coef = A.find ( rowC->first )->second[rowB->first];

                        for ( size_t j = 0; j != num_cols; ++j )
                        {
                            rowC->second[j] += A_coef * rowB->second[j];
                        }
                    }
                }

#if 0    // old implementation
                std::map<int, std::vector<DataType> > map_m2m = map_m2v;
                for ( std::map<int, std::vector<DataType> >::iterator it = map_m2m.begin ( ); it != map_m2m.end ( ); ++it )
                {
                    it->second.resize ( master_ansatz.get_nb_dof_on_cell ( ) );
                    for ( int i = 0; i < it->second.size ( ); ++i )
                        it->second[i] = 0.;

                    for ( int vindex = 0; vindex < map_m2v[it->first].size ( ); ++vindex )
                    {
                        DataType weight = map_m2v[it->first][vindex];
                        if ( std::abs ( weight ) > COEFFICIENT_EPS )
                        {
                            for ( int mindex = 0; mindex < map_v2m[vindex].size ( ); ++mindex )
                            {
                                it->second[mindex] += weight * map_v2m[vindex][mindex];
                            }
                        }
                    }
                }
#endif
            }

            template<class DataType>
            void find_constrained_dofs ( const typename NumberingLagrange<DataType>::InterfaceMatrix& matrix,
                                         DataType tol,
                                         std::vector<int>& constrained )
            {
                constrained.clear ( );
                constrained.reserve ( matrix.size ( ) );

                for ( typename NumberingLagrange<DataType>::InterfaceMatrix::const_iterator it = matrix.begin ( ), e_it = matrix.end ( ); it != e_it; ++it )
                {
                    for ( size_t i = 0, end = it->second.size ( ); i != end; ++i )
                    {
                        const DataType coef = it->second[i];

                        // A dof is unconstrained if the corresponding row
                        // contains exactly one entry == 1. Since the row sum is
                        // guaranteed to be 1., we only check if each
                        // coefficient is either 0 or 1 here.
                        if ( ( std::abs ( coef ) > tol ) && ( std::abs ( coef - 1. ) > tol ) )
                        {
                            constrained.push_back ( it->first );
                            break;
                        }
                    }
                }
            }

            template<class DataType>
            void correct_number_unconstrained_master_dofs ( typename NumberingLagrange<DataType>::InterfaceMatrix& matrix,
                                                            DataType tol,
                                                            int target_num_unconstrained )
            {
                std::vector<int> constrained_dofs;
                find_constrained_dofs<DataType>( matrix, tol, constrained_dofs );

                const int num_unconstrained = matrix.size ( ) - constrained_dofs.size ( );

                LOG_DEBUG ( 2, "Current number unconstrained dofs = " << num_unconstrained
                            << "\nNeeded number unconstrained dofs = " << target_num_unconstrained );

                assert ( num_unconstrained <= target_num_unconstrained );

                // Sort the constrained dofs so that they can be eliminated in order of their ID:s.
                std::sort ( constrained_dofs.begin ( ), constrained_dofs.end ( ) );
                const int num_dofs_to_relax = target_num_unconstrained - num_unconstrained;

                for ( std::vector<int>::const_iterator it = constrained_dofs.begin ( ),
                      end = constrained_dofs.begin ( ) + num_dofs_to_relax; it != end; ++it )
                {
                    typename NumberingLagrange<DataType>::InterfaceMatrix::iterator row = matrix.find ( *it );
                    assert ( row != matrix.end ( ) );

                    // Replace dependencies with identity row.
                    row->second.assign ( row->second.size ( ), 0. );
                    row->second[row->first] = 1.;
                }
            }

            template<class DataType>
            void add_master_dof_interpolation_pattern ( const typename NumberingLagrange<DataType>::InterfaceMatrix& matrix,
                                                        DataType tol,
                                                        DofInterpolationPattern& pattern )
            {
                // Only add interpolations for constrained dofs.
                std::vector<int> constrained_dofs;
                find_constrained_dofs<DataType>( matrix, tol, constrained_dofs );

                std::vector< std::pair<int, double> > dependencies;

                for ( std::vector<int>::const_iterator it = constrained_dofs.begin ( ),
                      end = constrained_dofs.end ( ); it != end; ++it )
                {
                    dependencies.clear ( );

                    typename NumberingLagrange<DataType>::InterfaceMatrix::const_iterator row_it = matrix.find ( *it );
                    assert ( row_it != matrix.end ( ) );

                    const int row = row_it->first;
                    const DataType row_coef = row_it->second[row];

                    for ( size_t col = 0, end_col = row_it->second.size ( ); col != end_col; ++col )
                    {
                        const DataType col_coef = row_it->second[col];

                        // Add dependencies from columns with non-zero
                        // weights, skipping the column corresponding to the constrained dof itself.
                        if ( col != row &&
                             ( std::abs ( col_coef ) > tol ) )
                        {
                            // Normalize weight with 1 - row_dof_coef, since
                            // the weight for the left-hand-side constrained
                            // dof might not be zero.
                            // TODO: should it not always be zero, so that weight is always one?
                            const double weight = static_cast < double > ( col_coef / ( 1. - row_coef ) );
                            dependencies.push_back ( std::make_pair ( col, weight ) );
                        }
                    }
                    pattern.insert_interpolation_master ( std::make_pair ( row, dependencies ) );
                }
            }

            template<class DataType>
            void add_slave_interpolation_pattern ( int slave,
                                                   const typename NumberingLagrange<DataType>::InterfaceMatrix& matrix,
                                                   DataType tol,
                                                   DofInterpolationPattern& pattern )
            {
                // On slave, all dofs are constrained, so no need to extract
                // them specifically.

                std::vector< std::pair<int, double> > dependencies;

                for ( typename NumberingLagrange<DataType>::InterfaceMatrix::const_iterator row_it = matrix.begin ( ), end_it = matrix.end ( );
                      row_it != end_it; ++row_it )
                {
                    dependencies.clear ( );

                    const int row = row_it->first;

                    for ( size_t col = 0, end_col = row_it->second.size ( ); col != end_col; ++col )
                    {
                        const double col_coef = static_cast < double > ( row_it->second[col] );
                        // TODO: it is assumed here that the entry for which
                        // col == row will be zero, and so automatically
                        // skipped. No normalization is necessary here.
                        if ( std::abs ( col_coef ) > tol )
                        {
                            dependencies.push_back ( std::make_pair ( col, col_coef ) );
                        }
                    }
                    pattern.insert_interpolation_slave ( slave, std::make_pair ( row, dependencies ) );
                }
            }

            template void log_debug_fe_pattern_info<double>( const mesh::Entity& entity, const FEType<double>& fe_type, int gdim );
            template void log_debug_fe_pattern_info<float>( const mesh::Entity& entity, const FEType<float>& fe_type, int gdim );

            template void log_compute_interpolation_start<double>( const FEInterfacePattern<double>& pattern,
                    const mesh::Interface& interface,
                    const mesh::Mesh* mesh );
            template void log_compute_interpolation_start<float>( const FEInterfacePattern<float>& pattern,
                    const mesh::Interface& interface,
                    const mesh::Mesh* mesh );

            template void log_interface_matrix<double>( int level, const NumberingLagrange<double>::InterfaceMatrix& interface_matrix );
            template void log_interface_matrix<float>( int level, const NumberingLagrange<float>::InterfaceMatrix& interface_matrix );

            template void extract_interface_matrix<double>( const NumberingLagrange<double>::InterpolationWeights& interpolation_weights,
                    const std::vector<DofID>& interface_dofs,
                    NumberingLagrange<double>::InterfaceMatrix& interface_matrix );
            template void extract_interface_matrix<float>( const NumberingLagrange<float>::InterpolationWeights& interpolation_weights,
                    const std::vector<DofID>& interface_dofs,
                    NumberingLagrange<float>::InterfaceMatrix& interface_matrix );

            template void multiply_interface_matrices<double>( const NumberingLagrange<double>::InterfaceMatrix& A,
                    const NumberingLagrange<double>::InterfaceMatrix& B,
                    NumberingLagrange<double>::InterfaceMatrix& C );
            template void multiply_interface_matrices<float>( const NumberingLagrange<float>::InterfaceMatrix& A,
                    const NumberingLagrange<float>::InterfaceMatrix& B,
                    NumberingLagrange<float>::InterfaceMatrix& C );

            template void find_constrained_dofs<double>( const NumberingLagrange<double>::InterfaceMatrix& matrix,
                    double tol,
                    std::vector<int>& constrained );
            template void find_constrained_dofs<float>( const NumberingLagrange<float>::InterfaceMatrix& matrix,
                    float tol,
                    std::vector<int>& constrained );

            template void correct_number_unconstrained_master_dofs<double>( NumberingLagrange<double>::InterfaceMatrix& matrix,
                    double tol,
                    int target_num_unconstrained );
            template void correct_number_unconstrained_master_dofs<float>( NumberingLagrange<float>::InterfaceMatrix& matrix,
                    float tol,
                    int target_num_unconstrained );

            template void add_master_dof_interpolation_pattern<double>( const NumberingLagrange<double>::InterfaceMatrix& matrix,
                    double tol,
                    DofInterpolationPattern& pattern );
            template void add_master_dof_interpolation_pattern<float>( const NumberingLagrange<float>::InterfaceMatrix& matrix,
                    float tol,
                    DofInterpolationPattern& pattern );

            template void add_slave_interpolation_pattern<double>( int slave,
                    const NumberingLagrange<double>::InterfaceMatrix& matrix,
                    double tol,
                    DofInterpolationPattern& pattern );
            template void add_slave_interpolation_pattern<float>( int slave,
                    const NumberingLagrange<float>::InterfaceMatrix& matrix,
                    float tol,
                    DofInterpolationPattern& pattern );
        } // helper function namespace

        /// \brief computes the interpolation description for given interface mode
        /// \details the interpolation description is given by ...
        /// \param[in] pattern the FEInterfacePattern that describes the neighbouring
        ///                    cells at the considered interface
        /// \return the interface mode interpolation
        /// \see DofInterpolationPattern

        template<class DataType>
        void NumberingLagrange<DataType>::compute_interpolation ( const FEInterfacePattern<DataType>& pattern,
                                                                  const mesh::Interface& interface,
                                                                  DofInterpolationPattern& interpolation )
        {
            // Log general information for this interpolation.
            log_compute_interpolation_start<DataType>( pattern, interface, this->mesh_ );

            const mesh::TDim tdim = this->mesh_->tdim ( );

            // TODO: check reasonable tolerances
            /*DataType eps;
            if ( typeid (DataType ) == typeid (double ) )
            {
                eps = 1.0e-14;
            }
            else if ( typeid (DataType ) == typeid (float ) )
            {
                eps = 1.0e-6f;
            }
            else
            {
                LOG_ERROR ( "Unknown data type: neither double nor float." );
                exit ( -1 );
             */

            const DataType COEFFICIENT_EPS = 1.e3 * std::numeric_limits<DataType>::epsilon ( );

            // Prepare output argument.
            interpolation.set_number_slaves ( pattern.num_slaves ( ) );

            // Find degree and ansatz of interface.
            int which_slave = -1;
            const int interface_deg = pattern.get_interface_degree ( &which_slave );
            const FEType<DataType>* virtual_ansatz;
            if ( which_slave == -1 )
            {
                virtual_ansatz = &pattern.fe_type_master ( );
            }
            else
            {
                virtual_ansatz = pattern.fe_type_slaves ( )[which_slave];
            }

            LOG_DEBUG ( 2, "Interface degree = " << interface_deg );
            if ( which_slave == -1 )
            {
                LOG_DEBUG ( 2, "Dominating cell: master" );
            }
            else
            {
                LOG_DEBUG ( 2, "Dominating cell: slave " << which_slave );
            }
            LOG_DEBUG ( 2, "Virtual FE Ansatz = " << virtual_ansatz->get_name ( ) );

            // ***************************************************************************
            // Master <-> Master interpolation
            // ***************************************************************************

            ////////////////
            // 1. Compute master -> master interpolation matrix, as product of
            // master -> virtual and virtual -> master interface matrices.

            const mesh::Entity& master_cell = this->mesh_->get_entity ( tdim, interface.master_index ( ) );
            const FEType<DataType>& master_ansatz = pattern.fe_type_master ( );

            // Master to virtual mapping.
            InterpolationWeights weights_m2v;
            compute_weights ( master_cell,
                              *virtual_ansatz,
                              master_cell,
                              master_ansatz,
                              weights_m2v );

            // Get master DoFs lying on the interface.
            const std::vector<DofID>& master_dofs_on_interface =
                    master_ansatz.get_dof_on_subentity ( tdim - 1, pattern.master_facet_number ( ) );

            LOG_DEBUG ( 3, "Master dofs on interface = "
                        << string_from_range ( master_dofs_on_interface.begin ( ),
                                               master_dofs_on_interface.end ( ) ) );

            InterfaceMatrix m2v;
            extract_interface_matrix<DataType>( weights_m2v, master_dofs_on_interface, m2v );

            LOG_DEBUG ( 2, "m2v =" );
            log_interface_matrix<DataType>( 2, m2v );

            // Virtual to master mapping.
            InterpolationWeights weights_v2m;
            compute_weights ( master_cell,
                              master_ansatz,
                              master_cell,
                              *virtual_ansatz,
                              weights_v2m );

            // Get virtual DoFs lying on the interface.
            const std::vector<DofID>& virtual_dofs_on_interface =
                    virtual_ansatz->get_dof_on_subentity ( tdim - 1, pattern.master_facet_number ( ) );

            LOG_DEBUG ( 3, "Virtual dofs on interface = "
                        << string_from_range ( virtual_dofs_on_interface.begin ( ),
                                               virtual_dofs_on_interface.end ( ) ) );

            InterfaceMatrix v2m;
            extract_interface_matrix<DataType>( weights_v2m, virtual_dofs_on_interface, v2m );

            LOG_DEBUG ( 2, "v2m =" );
            log_interface_matrix<DataType>( 2, v2m );

            // Multiply m2v and v2m to get m2m (master -> master interpolation).
            InterfaceMatrix m2m;
            multiply_interface_matrices<DataType>( m2v, v2m, m2m );

            LOG_DEBUG ( 2, "m2m =" );
            log_interface_matrix<DataType>( 2, m2m );

            ////////////////
            // 2. Correct the number of unconstrained "real" DoFs on interface.
            const int needed_num_unconstrained_dofs =
                    virtual_ansatz->get_nb_dof_on_subentity ( tdim - 1, pattern.master_facet_number ( ) );
            correct_number_unconstrained_master_dofs<DataType>( m2m, COEFFICIENT_EPS, needed_num_unconstrained_dofs );

            LOG_DEBUG ( 2, "corrected m2m = " );
            log_interface_matrix<DataType>( 2, m2m );

            ////////////////
            // 3. Insert master interpolation into output argument.
            add_master_dof_interpolation_pattern<DataType>( m2m, COEFFICIENT_EPS, interpolation );

            // ***************************************************************************
            // Master -> Slave interpolation
            // ***************************************************************************

            for ( size_t s = 0, e_s = pattern.num_slaves ( ); s != e_s; ++s )
            {
                ////////////////
                // 1. Compute slave -> master interpolation matrix as product
                // of slave -> virtual and virtual -> master interface matrices.

                const mesh::Entity& slave_cell = this->mesh_->get_entity ( tdim,
                                                                           interface.slave_index ( s ) );
                const FEType<DataType>& slave_ansatz = *pattern.fe_type_slaves ( )[s];

                // Slave to virtual mapping.
                std::vector<std::vector<DataType> > weights_s2v;
                compute_weights ( master_cell,
                                  *virtual_ansatz,
                                  slave_cell,
                                  slave_ansatz,
                                  weights_s2v );

                // Get slave DoFs lying on the interface.
                std::vector<int> slave_dofs_on_interface =
                        slave_ansatz.get_dof_on_subentity ( tdim - 1, pattern.slave_facet_number ( s ) );

                LOG_DEBUG ( 3, "Slave " << s << " dofs on interface = "
                            << string_from_range ( slave_dofs_on_interface.begin ( ),
                                                   slave_dofs_on_interface.end ( ) ) );

                InterfaceMatrix s2v;
                extract_interface_matrix<DataType>( weights_s2v, slave_dofs_on_interface, s2v );

                LOG_DEBUG ( 2, "Slave " << s << " s2v =" );
                log_interface_matrix<DataType>( 2, s2v );

                // Multiply s2v and v2m to get s2m (slave -> master interpolation.
                InterfaceMatrix s2m;
                multiply_interface_matrices<DataType>( s2v, v2m, s2m );

                LOG_DEBUG ( 2, "Slave " << s << " s2m =" );
                log_interface_matrix<DataType>( 2, s2m );

                ////////////////
                // 2. Insert slave interpolation into output argument.
                add_slave_interpolation_pattern<DataType>( s, s2m, COEFFICIENT_EPS, interpolation );

            }

            LOG_DEBUG ( 2, "Interpolation pattern computed: \n" << interpolation );
        }

        namespace
        {

            template<class DataType>
            std::vector<DataType> map_coords ( int gdim,
                                               const CellTransformation<DataType>& T,
                                               const std::vector<DataType>& pt )
            {
                std::vector<DataType> res ( gdim, 0. );

                switch ( gdim )
                {
                    case 1:
                        res[0] = T.x ( pt );
                        break;
                    case 2:
                        res[0] = T.x ( pt );
                        res[1] = T.y ( pt );
                        break;
                    case 3:
                        res[0] = T.x ( pt );
                        res[1] = T.y ( pt );
                        res[2] = T.z ( pt );
                        break;
                    default:
                        throw "Dimension not supported\n";
                }

                return res;
            }

            template<class DataType>
            std::vector<DataType> inv_map_coords ( int gdim,
                                                   const CellTransformation<DataType>& T,
                                                   const std::vector<DataType>& pt )
            {
                std::vector<DataType> res ( gdim, 0. );

                switch ( gdim )
                {
                    case 1:
                        T.inverse ( pt[0], res[0] );
                        break;
                    case 2:
                        T.inverse ( pt[0], pt[1], res[0], res[1] );
                        break;
                    case 3:
                        T.inverse ( pt[0], pt[1], pt[2], res[0], res[1], res[2] );
                        break;
                    default:
                        throw "Dimension not supported\n";
                }

                return res;
            }

            template std::vector<double> map_coords ( int gdim,
                                                      const CellTransformation<double>& T,
                                                      const std::vector<double>& pt );
            template std::vector<float> map_coords ( int gdim,
                                                     const CellTransformation<float>& T,
                                                     const std::vector<float>& pt );

            template std::vector<double> inv_map_coords ( int gdim,
                                                          const CellTransformation<double>& T,
                                                          const std::vector<double>& pt );
            template std::vector<float> inv_map_coords ( int gdim,
                                                         const CellTransformation<float>& T,
                                                         const std::vector<float>& pt );
        }

        /// for all DoFs of (cellB, ansatzB) the weights of DoFs of (cellA, ansatzA)
        /// are calculated

        template<class DataType>
        void NumberingLagrange<DataType>::compute_weights ( mesh::Entity const& cellA,
                                                            doffem::FEType<DataType> const& ansatzA,
                                                            mesh::Entity const& cellB,
                                                            doffem::FEType<DataType> const& ansatzB,
                                                            InterpolationWeights& weights ) const
        {
            // Absolute tolerance for setting a weight to 0 or 1.

            // TODO: check reasonable tolerances
            /*DataType eps;
            if ( typeid (DataType ) == typeid (double ) )
            {
                eps = 1.0e-10;
            }
            else if ( typeid (DataType ) == typeid (float ) )
            {
                eps = 1.0e-5f;
            }
            else
            {
                LOG_ERROR ( "Unknown data type: neither double nor float." );
                exit ( -1 );
            }*/

            const DataType COEFFICIENT_EPS = 1.e3 * std::numeric_limits<DataType>::epsilon ( );

            const CellTransformation<DataType>* transA = this->fe_manager_->get_cell_transformation ( cellA.index ( ) );
            const CellTransformation<DataType>* transB = this->fe_manager_->get_cell_transformation ( cellB.index ( ) );

            const int gdim = this->mesh_->gdim ( );

            const int row_length = ansatzA.get_nb_dof_on_cell ( );

            // Cell B's dof points (reference cell).
            const std::vector< std::vector<DataType> >& dof_coords_B = ansatzB.get_coord_on_cell ( );
            weights.resize ( dof_coords_B.size ( ) );

            for ( size_t i = 0, e_i = dof_coords_B.size ( ); i != e_i; ++i )
            {
                LOG_DEBUG ( 3, "Dof point " << i << " reference cell B = "
                            << string_from_range ( dof_coords_B[i].begin ( ), dof_coords_B[i].end ( ) ) );

                // Compute coordinates of cell B's DoFs (physical cell).
                const std::vector<DataType> mapped_dof_pt = map_coords ( gdim, *transB, dof_coords_B[i] );

                LOG_DEBUG ( 3, "Dof point " << i << " physical cell B = "
                            << string_from_range ( mapped_dof_pt.begin ( ), mapped_dof_pt.end ( ) ) );

                // Map coordinates back to reference cell of A.
                const std::vector<DataType> inv_mapped_dof_pt = inv_map_coords ( gdim, *transA, mapped_dof_pt );

                LOG_DEBUG ( 3, "Dof point " << i << " ref cell A = "
                            << string_from_range ( inv_mapped_dof_pt.begin ( ), inv_mapped_dof_pt.end ( ) ) );

                // Calculate weights by evaluating A:s shape functions at the
                // inversely mapped pt. The rows in the matrix correspond to
                // the points, and the columns to the shape functions.
                weights[i].resize ( row_length, 0. );
                ansatzA.N ( inv_mapped_dof_pt, weights[i] );
            }

            // Filter weights -- make sure we get zeros and ones exactly correct.
            for ( size_t i = 0, e_i = weights.size ( ); i != e_i; ++i )
            {
                for ( size_t j = 0, e_j = row_length; j != e_j; ++j )
                {
                    DataType& w = weights[i][j];
                    if ( std::abs ( w ) < COEFFICIENT_EPS )
                    {
                        w = 0.;
                    }
                    else if ( std::abs ( w - 1. ) < COEFFICIENT_EPS )
                    {
                        w = 1.;
                    }
                }
            }
        }

        template class NumberingLagrange<double>;
        template class NumberingLagrange<float>;

    }
}
