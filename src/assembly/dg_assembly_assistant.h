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

#ifndef _DG_ASSEMBLY_ASSISTANT_H_
#    define _DG_ASSEMBLY_ASSISTANT_H_

#    include "assembly/assembly_assistant.h"
#    include "assembly/dg_assembly.h"

namespace hiflow
{
    using doffem::FEType;

    /// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok

    /// \brief Provides assembly functionalities needed for DG, i.e., local
    /// assembly over interfaces where function values are needed from 
    /// several adjacent cells.

    template<int DIM, class DataType>
    class DGAssemblyAssistant
    {
      public:
        typedef typename DGGlobalAssembler<DataType>::InterfaceSide InterfaceSide;
        typedef la::SeqDenseMatrix<DataType> LocalMatrix;
        typedef std::vector<DataType> LocalVector;

        /// \brief Standard constructor
        DGAssemblyAssistant ( );

        /// \brief Initialize an interface between two elements (trial_elem, test_elem).
        /// The physical points of the two quadratures must be the same.
        /// \param trial_elem Trial space finite element
        /// \param test_elem Test space finite element
        /// \param trial_quad Quadrature on trial space finite element
        /// \param test_quad Quadrature on test space finite element
        /// \param trial_facet_number Facet number of trial space finite element
        /// \param test_facet_number Facet number of test space finite element
        /// \param trial_if_side Kind of interface of trial space finite element
        /// \param test_if_side Kind of interface of test space finite element
        void initialize_for_interface ( const Element<DataType>& trial_elem, const Element<DataType>& test_elem,
                                        const Quadrature<DataType>& trial_quad, const Quadrature<DataType>& test_quad,
                                        int trial_facet_number, int test_facet_number,
                                        InterfaceSide trial_if_side, InterfaceSide test_if_side );
        // may change names: left = trial, right = test
        /// \brief Get AssemblyAssistant corresponding to trial space finite element

        AssemblyAssistant<DIM, DataType>& trial ( )
        {
            assert ( this->trial_aa_ != NULL );
            return *( this->trial_aa_ );
        };
        /// \brief Get AssemblyAssistant corresponding to test space finite element

        AssemblyAssistant<DIM, DataType>& test ( )
        {
            assert ( this->test_aa_ != NULL );
            return *( this->test_aa_ );
        };

        // AssemblyAssistant functions that are the same for trial and test function
        /// \brief Get number of quadrature points on current interface; 
        /// Provided by test space finite element

        inline int num_quadrature_points ( ) const
        {
            assert ( this->test_aa_ != NULL );
            return this->test_aa_->num_quadrature_points ( );
        };
        /// \brief Get coordinates of a quadrature point on current interface; 
        /// Provided by test space finite element
        /// \param q Number of quadrature point

        inline Vec<DIM, DataType> x ( int q ) const
        {
            assert ( this->test_aa_ != NULL );
            return this->test_aa_->x ( q );
        };
        /// \brief Get quadrature weight at a quadrature point on current interface; 
        /// Provided by test space finite element
        /// \param q Number of quadrature point

        inline DataType w ( int q ) const
        {
            assert ( this->test_aa_ != NULL );
            return this->test_aa_->w ( q );
        };
        /// \brief Get surface element at a quadrature point on current interface; 
        /// Provided by test space finite element
        /// \param q Number of quadrature point

        inline DataType ds ( int q ) const
        {
            assert ( this->test_aa_ != NULL );
            assert ( this->trial_aa_ != NULL );
            if (
                 ( this->right_if_side_ == DGGlobalAssembler<DataType>::INTERFACE_SLAVE )
                 || ( this->right_if_side_ == DGGlobalAssembler<DataType>::INTERFACE_BOUNDARY )
                 || ( this->right_if_side_ == this->left_if_side_ )
                 )
            {
                return this->test_aa_->ds ( q );
            }
            else
            {
                assert ( this->left_if_side_ == DGGlobalAssembler<DataType>::INTERFACE_SLAVE );
                return this->trial_aa_->ds ( q );
            }
        }
        /// \brief Get number of variables 
        /// Provided by test space finite element

        inline int num_vars ( ) const
        {
            assert ( this->test_aa_ != NULL );
            return this->test_aa_->num_vars ( );
        };
      protected:
        /// AssemblyAssistant corresponding to "master" finite element
        AssemblyAssistant<DIM, DataType> master_aa_;
        /// AssemblyAssistant corresponding to "slave" finite element
        AssemblyAssistant<DIM, DataType> slave_aa_;
        /// AssemblyAssistant corresponding to trial space finite element
        AssemblyAssistant<DIM, DataType>* trial_aa_;
        /// AssemblyAssistant corresponding to test space finite element
        AssemblyAssistant<DIM, DataType>* test_aa_;

        /// Flag if current interface is a boundary facet
        bool is_boundary_;

        /// Interface sides of neighbouring cells
        InterfaceSide left_if_side_, right_if_side_;
    };

    template<int DIM, class DataType>
    DGAssemblyAssistant<DIM, DataType>::DGAssemblyAssistant ( )
    {
    }

    template<int DIM, class DataType>
    void DGAssemblyAssistant<DIM, DataType>::initialize_for_interface (
                                                                        const Element<DataType>& left_elem, const Element<DataType>& right_elem,
                                                                        const Quadrature<DataType>& left_quad, const Quadrature<DataType>& right_quad,
                                                                        int left_facet_number, int right_facet_number,
                                                                        InterfaceSide left_if_side, InterfaceSide right_if_side )
    {

        this->left_if_side_ = left_if_side;
        this->right_if_side_ = right_if_side;
        this->is_boundary_ = ( right_if_side == DGGlobalAssembler<DataType>::INTERFACE_BOUNDARY );
        // Convention: left -- trial  functions, right -- test functions. TODO: Maybe change this vice versa
        this->master_aa_.initialize_for_facet ( right_elem, right_quad, right_facet_number );
        this->test_aa_ = &( this->master_aa_ );
        // only init both facets if they are different.
        if ( this->is_boundary_ || left_if_side == right_if_side )
        { //if left elem == right elem
            this->trial_aa_ = this->test_aa_;
        }
        else
        {
            this->slave_aa_.initialize_for_facet ( left_elem, left_quad, left_facet_number );
            this->trial_aa_ = &( this->slave_aa_ );
        }

        assert ( this->trial_aa_->num_quadrature_points ( ) == this->test_aa_->num_quadrature_points ( ) );
    }
} //namespace hiflow
#endif
