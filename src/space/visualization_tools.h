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

#ifndef HIFLOW_SPACE_VISUALIZATION_TOOLS
#    define HIFLOW_SPACE_VISUALIZATION_TOOLS

/// \author Staffan Ronnas, Martin Baumann, Teresa Beck

#    include <map>
#    include <string>
#    include <vector>

#    include <boost/function.hpp>
#    include <mpi.h>
#    include "common/permutation.h"
#    include "space/pp_vector.h"
#    include "space/vector_space.h"
#    include "mesh/mesh.h"
#    include "mesh/entity.h"
#    include "linear_algebra/la_descriptor.h"
#    include "dof/dof_fem_types.h"
#    include "common/vector_algebra.h"
#    include "assembly/assembly_assistant_values.h"

using namespace hiflow::mesh;

namespace hiflow
{

    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    template<class LAD>
    class EvalFeFunction
    {
        typedef hiflow::doffem::DofID DofID;

        // Type of function for evaluation.
        typedef boost::function3<void,
        const mesh::Entity&, // cell
        const std::vector<typename LAD::DataType>&, // reference coordinates
        std::vector<typename LAD::DataType>& // values of function at the points
        > EvalFunction;

      public:

        EvalFeFunction ( const VectorSpace<typename LAD::DataType>& space,
                         const typename LAD::VectorType& fun, int var = 0 )
        : space_ ( space ), var_ ( var ), fun_ ( fun )
        {
            // sort (id,val) by ordering of id once,
            // but only if we are in parallel mode
            static int num_procs;
            MPI_Comm_size ( space_.get_mpi_comm ( ), &num_procs );
            if ( num_procs != 1 )
            {
                // get data
                std::vector<DofID> id_tmp;
                std::vector<typename LAD::DataType> val_tmp;

                fun_.GetAllDofsAndValues ( id_tmp, val_tmp );

                // calculate permutation
                std::vector<int> permutation;
                compute_sorting_permutation ( id_tmp, permutation );

                // permute
                permute_vector ( permutation, id_tmp, id_ );
                permute_vector ( permutation, val_tmp, val_ );
            }
        }

        void operator() ( const Entity& cell,
                const std::vector<typename LAD::DataType>& ref_coords,
                std::vector<typename LAD::DataType>& values ) const
        {
            const int gdim = space_.mesh ( ).gdim ( );
            const int num_points = ref_coords.size ( ) / gdim;

            std::vector<int> global_dof_ids;
            space_.GetDofIndices ( var_, cell, &global_dof_ids );
            const int num_dofs = global_dof_ids.size ( );

            std::vector < std::vector<typename LAD::DataType> > shape_fun ( num_points, std::vector<typename LAD::DataType > ( num_dofs, 1.e13 ) );
            std::vector<typename LAD::DataType > pt ( gdim, 0. );

            int k = 0;
            for ( int i = 0; i < num_points; ++i )
            {
                for ( int c = 0; c < gdim; ++c )
                    pt[c] = ref_coords[k++];
                space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), var_ )->N ( pt, shape_fun[i] );
            }

            std::vector<typename LAD::DataType > dof_values ( num_dofs, 1.e25 );

            static int num_procs;
            MPI_Comm_size ( space_.get_mpi_comm ( ), &num_procs );

            if ( num_procs == 1 )
            {
                // in sequential world, DofIDs are already sorted
                fun_.GetValues ( &( global_dof_ids[0] ), global_dof_ids.size ( ), &( dof_values[0] ) );
            }
            else
            {
                // in parallel world, DofIDs are not sorted
                // -> id and val fields need to be sorted and accesses are related to
                //    a seek through the data

                std::vector<DofID>::const_iterator it;
                for ( int i = 0; i < num_dofs; ++i )
                {
                    it = std::lower_bound ( id_.begin ( ), id_.end ( ), global_dof_ids[i] );
                    const int index = it - id_.begin ( );
                    dof_values[i] = val_[index];
                }

                // slow version
                //fun_.GetValues(&global_dof_ids[0], num_dofs, &dof_values[0]);
            }

            values.clear ( );
            values.resize ( num_points, 0. );
            for ( int i = 0; i < num_points; ++i )
            {
                values[i] = 0.;
                for ( int j = 0; j < num_dofs; ++j )
                    values[i] += shape_fun[i][j] * dof_values[j];
            }
        }

      private:
        const VectorSpace<typename LAD::DataType>& space_;
        int var_;
        const typename LAD::VectorType& fun_;

        // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
        std::vector<DofID> id_;
        std::vector<typename LAD::DataType> val_;
    };

    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    template<class LAD, int DIM>
    class EvalDerivativeFeFunction
    {
        typedef hiflow::doffem::DofID DofID;

        // Type of function for evaluation.
        typedef boost::function3<void,
        const mesh::Entity&, // cell
        const std::vector<typename LAD::DataType>&, // reference coordinates
        std::vector< typename LAD::DataType >& // values of function at the points
        > EvalFunction;

      public:

        EvalDerivativeFeFunction ( const VectorSpace<typename LAD::DataType>& space,
                                   const typename LAD::VectorType& fun, int var = 0, int deriv_dir = 0 )
        : space_ ( space ), var_ ( var ), fun_ ( fun ), dir_ ( deriv_dir )
        {
            // sort (id,val) by ordering of id once,
            // but only if we are in parallel mode
            static int num_procs;
            MPI_Comm_size ( space_.get_mpi_comm ( ), &num_procs );
            if ( num_procs != 1 )
            {
                // get data
                std::vector<DofID> id_tmp;
                std::vector< typename LAD::DataType > val_tmp;

                fun_.GetAllDofsAndValues ( id_tmp, val_tmp );

                // calculate permutation
                std::vector<int> permutation;
                compute_sorting_permutation ( id_tmp, permutation );

                // permute
                permute_vector ( permutation, id_tmp, id_ );
                permute_vector ( permutation, val_tmp, val_ );
            }
        }

        void operator() ( const Entity& cell,
                const std::vector<typename LAD::DataType>& ref_coords,
                std::vector< typename LAD::DataType >& values ) const
        {
            const int gdim = space_.mesh ( ).gdim ( );
            const int num_points = ref_coords.size ( ) / gdim;

            std::vector<int> global_dof_ids;
            space_.GetDofIndices ( var_, cell, &global_dof_ids );
            const int num_dofs = global_dof_ids.size ( );

            // store reference points
            std::vector < Vec< DIM, typename LAD::DataType > > q_points ( num_points );
            int k = 0;
            for ( int i = 0; i < num_points; ++i )
            {
                for ( int c = 0; c < gdim; ++c )
                {
                    q_points[i][c] = ref_coords[k++];
                }
            }

            // get FE type
            std::vector< const FEType< typename LAD::DataType > *> fe_types;
            fe_types.push_back ( space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), var_ ) );

            // compute shape function gradients
            FunctionValues < std::vector < Vec< DIM, typename LAD::DataType > > > grad_phi_hat;
            grad_phi_hat.compute ( q_points, EvalShapeFunctionGradients<DIM, typename LAD::DataType > ( fe_types ) );

            // compute jacobians on physical cell
            FunctionValues < Mat< DIM, DIM, typename LAD::DataType > > J;
            const CellTransformation< typename LAD::DataType >& cell_transform = this->space_.GetCellTransformation ( cell );

            J.compute ( q_points, EvalPhysicalJacobian<DIM, typename LAD::DataType > ( cell_transform ) );

            // compute inverse-transpose of jacobians
            FunctionValues < Mat< DIM, DIM, typename LAD::DataType > > JinvT;
            JinvT.compute ( J, EvalInvTranspose<DIM, typename LAD::DataType > ( ) );

            // compute JinvT * grad_phi_hat
            FunctionValues < std::vector < Vec< DIM, typename LAD::DataType > > > grad_phi;
            grad_phi.compute ( JinvT, EvalMappedShapeFunctionGradients<DIM, typename LAD::DataType > ( grad_phi_hat ) );

            assert ( grad_phi.size ( ) == num_points );

            // get dofs
            std::vector<typename LAD::DataType > dof_values ( num_dofs, 1.e25 );

            static int num_procs;
            MPI_Comm_size ( space_.get_mpi_comm ( ), &num_procs );

            if ( num_procs == 1 )
            {
                // in sequential world, DofIDs are already sorted
                fun_.GetValues ( &( global_dof_ids[0] ), global_dof_ids.size ( ), &( dof_values[0] ) );
            }
            else
            {
                // in parallel world, DofIDs are not sorted
                // -> id and val fields need to be sorted and accesses are related to
                //    a seek through the data

                std::vector<DofID>::const_iterator it;
                for ( int i = 0; i < num_dofs; ++i )
                {
                    it = std::lower_bound ( id_.begin ( ), id_.end ( ), global_dof_ids[i] );
                    const int index = it - id_.begin ( );
                    dof_values[i] = val_[index];
                }

                // slow version
                //fun_.GetValues(&global_dof_ids[0], num_dofs, &dof_values[0]);
            }

            values.clear ( );
            values.resize ( num_points, 0. );
            for ( int i = 0; i < num_points; ++i )
            {
                values[i] = 0.;
                for ( int j = 0; j < num_dofs; ++j )
                {
                    values[i] += grad_phi[i][j][dir_] * dof_values[j];
                }
            }
        }

      private:
        const VectorSpace<typename LAD::DataType>& space_;
        int var_;
        int dir_;
        //hiflow::la::CoupledVector<typename LAD::DataType>& fun_;
        const typename LAD::VectorType& fun_;

        // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
        std::vector<DofID> id_;
        std::vector< typename LAD::DataType > val_;
    };

    template<class LAD>
    class EvalFeFunction2
    {
        typedef hiflow::doffem::DofID DofID;

        // Type of function for evaluation.
        typedef boost::function3<void,
        const mesh::Entity&, // cell
        const std::vector<typename LAD::DataType>&, // reference coordinates
        std::vector<typename LAD::DataType>& // values of function at the points
        > EvalFunction;

      public:

        EvalFeFunction2 ( const VectorSpace<typename LAD::DataType>& space,

                          const hiflow::la::PpVector<LAD>& fun, int var = 0 )
        : space_ ( space ), var_ ( var ), fun_ ( fun )
        {
            // sort (id,val) by ordering of id once,
            // but only if we are in parallel mode
            static int num_procs;
            MPI_Comm_size ( space_.get_mpi_comm ( ), &num_procs );
            if ( num_procs != 1 )
            {
                // get data
                std::vector<DofID> id_tmp;
                std::vector<typename LAD::DataType> val_tmp;
                fun.GetDofsAndValues ( id_tmp, val_tmp );

                // calculate permutation
                std::vector<int> permutation;
                compute_sorting_permutation ( id_tmp, permutation );

                // permute
                permute_vector ( permutation, id_tmp, id_ );
                permute_vector ( permutation, val_tmp, val_ );
            }
        }

        void operator() ( const Entity& cell,
                const std::vector<typename LAD::DataType>& ref_coords,
                std::vector<typename LAD::DataType>& values ) const
        {
            const int gdim = space_.mesh ( ).gdim ( );
            const int num_points = ref_coords.size ( ) / gdim;

            std::vector<int> global_dof_ids;
            space_.GetDofIndices ( var_, cell, &global_dof_ids );
            const int num_dofs = global_dof_ids.size ( );

            std::vector < std::vector<typename LAD::DataType> > shape_fun ( num_points, std::vector<typename LAD::DataType > ( num_dofs, 1.e13 ) );
            std::vector<typename LAD::DataType > pt ( gdim, 0. );

            int k = 0;
            for ( int i = 0; i < num_points; ++i )
            {
                for ( int c = 0; c < gdim; ++c )
                    pt[c] = ref_coords[k++];
                space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), var_ )->N ( pt, shape_fun[i] );
            }
            std::vector<typename LAD::DataType > dof_values ( global_dof_ids.size ( ), 1.e25 );

            static int num_procs;
            MPI_Comm_size ( space_.get_mpi_comm ( ), &num_procs );
            if ( num_procs == 1 )
            {
                // in sequential world, DofIDs are already sorted
                // -> direct accesses are possible

                for ( int i = 0; i < global_dof_ids.size ( ); ++i )
                    dof_values[i] = fun_[global_dof_ids[i]];
            }
            else
            {
                // in parallel world, DofIDs are not sorted
                // -> id and val fields need to be sorted and accesses are related to
                //    a seek through the data

                std::vector<DofID>::const_iterator it;
                for ( int i = 0; i < num_dofs; ++i )
                {
                    it = std::lower_bound ( id_.begin ( ), id_.end ( ), global_dof_ids[i] );
                    const int index = it - id_.begin ( );
                    dof_values[i] = val_[index];
                }

                // slow version
                //fun_.GetValues(&global_dof_ids[0], num_dofs, &dof_values[0]);
            }

            values.clear ( );
            values.resize ( num_points, 0. );

            for ( int i = 0; i < num_points; ++i )
            {
                values[i] = 0.;
                for ( int j = 0; j < num_dofs; ++j )
                    values[i] += shape_fun[i][j] * dof_values[j];
            }
        }

      private:
        const VectorSpace<typename LAD::DataType>& space_;
        int var_;
        const hiflow::la::PpVector<LAD>& fun_;

        // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
        std::vector<DofID> id_;
        std::vector<typename LAD::DataType> val_;
    };

}

#endif
