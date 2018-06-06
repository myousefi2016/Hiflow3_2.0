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

#ifndef __FEM_FEMANAGER_H_
#    define __FEM_FEMANAGER_H_

#    include <cassert>
#    include <vector>
#    include "common/macros.h"

#    include "mesh/mesh.h"

#    include "dof/dof_fem_types.h"
#    include "fem/fetype.h"
#    include "fem/feinstance.h"
#    include "fem/cell_transformation.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class FEManager femanager.h
        /// \brief Ancestor Manager of information about the finite element ansatz
        /// \author Michael Schick<br>Martin Baumann
        ///

        template<class DataType>
        class FEManager
        {
          public:

            typedef std::vector<DataType> Coord;

            /// Constructor with parameters geometrical dimension and nb of variables
            explicit FEManager ( int dim, int nb_var );

            /// Destructor
            virtual ~FEManager ( );

            /// Initialize reference to given mesh and automatically storing cell transformations
            void set_mesh ( const mesh::Mesh& mesh );

            /// Information if a continuous (true) or discontinuous (false) ansatz is used for a specific variable
            inline void set_ca ( int var, bool mode );

            /// Get dimension
            inline int get_dim ( ) const;
            /// Get total number of variables
            inline int get_nb_var ( ) const;

            /// Information if a continuous (true) or discontinuous (false) ansatz is used for a specific variable
            inline bool get_ca ( int var ) const;

            /// Information about the used FEType on a given mesh cell for a specific variable
            inline FEType<DataType>* get_fe_on_cell ( int cell_index, int var ) const;

            /// Get the Id of a Finite Element on a specific cell_index and variable
            inline FETypeID get_fetype_id ( int cell_index, int var ) const;

            /// Get cell transformation for a specific mesh cell
            inline CellTransformation<DataType>* get_cell_transformation ( int cell_index ) const;

            /// Status information
            void get_status ( ) const;
            /// Get number of initialized Finite Elements

            int nfe ( ) const
            {
                return fe_inst_.nfe ( );
            }

            /// \brief Initialize FE Tank by given fe ansatz and parameters
            ///        for all variables (all cells get same ansatz and parameters)
            virtual void init_fe_tank ( const typename FEType<DataType>::FEAnsatz& ansatz, const std::vector<std::vector<int> >& param );
            /// \brief Initialize FE Tank by given fe ansatz and parameters
            ///        for a specific variable (all cells get same ansatz)
            virtual void init_fe_tank ( int var, const typename FEType<DataType>::FEAnsatz& ansatz, const std::vector<int>& param );

          protected:

            /// Initialize Linear Cell Transformations by given coordinates of mesh on cells
            void init_cell_transformation ( );

            /// Dimension
            int dim_;

            /// Total number of variables
            int nb_var_;

            /// This vector stores for every variable if fe_tank was initialized
            std::vector<bool> initialized_;

            /// Continuous (true) or discontinuous (false) ansatz for each variable
            std::vector<bool> ca_;

            /// Const pointer to given mesh
            const mesh::Mesh* mesh_;
            /// Number of cells on the mesh
            int num_entities_;

            /// FE Tank, which holds pointers to every FEType for every cell and variable
            std::vector<FEType<DataType>* > fe_tank_;

            /// Instance holder for every needed FEType
            FEInstance<DataType> fe_inst_;

            /// Cell transformations for every mesh cell
            std::vector<CellTransformation<DataType>* > cell_transformation_;
        };

        //------------ INLINE FUNCTIONS FOR FEMANAGER ---------------

        template<class DataType>
        void FEManager<DataType>::set_ca ( int var, bool mode )
        {
            interminable_assert ( var >= 0 && var < nb_var_ );
            ca_[var] = mode;
        }

        template<class DataType>
        int FEManager<DataType>::get_dim ( ) const
        {
            return dim_;
        }

        template<class DataType>
        int FEManager<DataType>::get_nb_var ( ) const
        {
            return nb_var_;
        }

        template<class DataType>
        bool FEManager<DataType>::get_ca ( int var ) const
        {
            interminable_assert ( var >= 0 && var < nb_var_ );
            return ca_[var];
        }

        template<class DataType>
        CellTransformation<DataType>* FEManager<DataType>::get_cell_transformation ( int cell_index ) const
        {
            assert ( cell_index >= 0 && cell_index < cell_transformation_.size ( ) );
            return cell_transformation_[cell_index];
        }

        /// \details Returns the finite element ansatz for the desired cell and variable

        template<class DataType>
        FEType<DataType>* FEManager<DataType>::get_fe_on_cell ( int cell_index, int var ) const
        {
            // Calculate fe_tank_ index for given cell_id and var
            int idx = var * num_entities_ + cell_index;

            // Get FE ansatz on cell for given var
            return fe_tank_[idx];
        }

        /// \details The ID is stored in FEType \see FEType::FiniteElement

        template<class DataType>
        FETypeID FEManager<DataType>::get_fetype_id ( int cell_index, int var ) const
        {
            return get_fe_on_cell ( cell_index, var )->get_my_id ( );
        }

    } // namespace doffem
} // namespace hiflow
#endif
