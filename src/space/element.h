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

#ifndef HIFLOW_ELEMENT_H
#    define HIFLOW_ELEMENT_H

#    include "fem/cell_transformation.h"
#    include "mesh/entity.h"
#    include "fem/fetype.h"
#    include "space/vector_space.h"

/// @author Staffan Ronnas

namespace hiflow
{

    template<class T> class Quadrature;

    /// \brief Class representing one element of a finite element space.
    ///
    /// \details This class provides a local view of one element of
    /// the finite element space, and brings together information from
    /// the mesh, dof and fem data structures.

    template<class DataType>
    class Element
    {
      public:

        /// \brief Construct an element on a given cell in a space.
        ///
        /// \param[in] space         finite element space to which element belongs
        /// \param[in] cell_index    index of cell for which element should be created

        Element ( const VectorSpace<DataType>& space, int cell_index )
        : space_ ( &space ),
        cell_ ( space.mesh ( ).get_entity ( space.get_dim ( ), cell_index ) )
        {
        };

        /// \return the cell transformation associated with the element

        const doffem::CellTransformation<DataType>* get_cell_transformation ( ) const
        {
            return &this->space_->GetCellTransformation ( this->cell_ );
        }

        /// \brief Accesss the finite element ansatz for a variable on the cell.
        /// \param[in] var    the number of the variable
        /// \return a pointer to the finite element ansatz for variable var on the cell.

        const doffem::FEType<DataType>* get_fe_type ( int var ) const
        {
            return this->space_->fe_manager ( ).get_fe_on_cell ( this->cell_.index ( ), var );
        }

        /// \brief Access the dof indices associated with the element.
        /// \param[out] indices   vector of dof indices associated with the element.

        void get_dof_indices ( std::vector<int>& indices ) const
        {
            this->space_->GetDofIndices ( this->cell_, &indices );
        }

        /// \brief Access the dof indices on the boundary for a given variable associated with the element.
        /// \param[in]  var       number of the variable
        /// \param[out] indices   vector of dof indices for variable @p var associated with the element.

        void get_dof_indices_on_facet ( int var, int tdim, int sindex, std::vector<int>& indices ) const
        {
            this->space_->GetDofIndicesFacet ( var, this->cell_, tdim, sindex, &indices );
        }

        /// \brief Access the dof indices on the boundary associated with the element.
        /// \param[out] indices   vector of dof indices associated with the element.

        void get_dof_indices_on_facet ( int tdim, int sindex, std::vector<int>& indices ) const
        {
            this->space_->GetDofIndicesFacet ( this->cell_, tdim, sindex, &indices );
        }

        /// \brief Access the dof indices for a given variables associated with the element.
        /// \param[in]  var       number of the variable
        /// \param[out] indices   vector of dof indices for variable @p var associated with the element.

        void get_dof_indices ( int var, std::vector<int>& indices ) const
        {
            this->space_->GetDofIndices ( var, this->cell_, &indices );
        }

        /// \return the number of variables in the finite element space

        int get_num_variables ( ) const
        {
            return this->space_->fe_manager ( ).get_nb_var ( );
        }

        /// \brief Access the number of dofs for a given variable associated with the element.
        /// \param[in]  var   the number of the variable
        /// \return the number of dofs associated with variable var

        int get_num_dofs ( int var ) const
        {
            const doffem::FEType<DataType>* fe_type = get_fe_type ( var );
            assert ( fe_type != 0 );
            return fe_type->get_nb_dof_on_cell ( );
        }

        /// \return the cell index of the element.

        mesh::EntityNumber get_cell_index ( ) const
        {
            return this->cell_.index ( );
        }

        /// \return a reference to the cell entity associated with the element.

        const mesh::Entity& get_cell ( ) const
        {
            return this->cell_;
        }

        /// \return true if and only if the element is adjacent to the boundary of the mesh.

        bool is_boundary ( ) const
        {
            const mesh::TDim facet_dim = this->space_->mesh ( ).tdim ( ) - 1;
            for ( mesh::IncidentEntityIterator it = this->cell_.begin_incident ( facet_dim );
                  it != this->cell_.end_incident ( facet_dim ); ++it )
            {
                bool is_boundary = this->space_->mesh ( ).is_boundary_facet ( it->index ( ) );
                if ( is_boundary )
                {
                    return true;
                }
            }
            return false;
        }

        /// \return the local facet numbers of the element which lie on the boundary.

        std::vector<int> boundary_facet_numbers ( ) const
        {
            const mesh::TDim facet_dim = this->space_->mesh ( ).tdim ( ) - 1;
            std::vector<int> facet_numbers ( 0 );
            if ( is_boundary ( ) )
            {
                facet_numbers.reserve ( 6 );
                int i = 0;
                for ( mesh::IncidentEntityIterator it = this->cell_.begin_incident ( facet_dim );
                      it != this->cell_.end_incident ( facet_dim ); ++it )
                {
                    if ( this->space_->mesh ( ).is_boundary_facet ( it->index ( ) ) )
                    {
                        facet_numbers.push_back ( i );
                    }
                    ++i;
                }
            }
            return facet_numbers;
        }

        /// \return true if and only if the element belongs to the local subdomain.

        bool is_local ( ) const
        {
            const mesh::TDim tdim = this->space_->mesh ( ).tdim ( );
            if ( this->space_->mesh ( ).has_attribute ( "_sub_domain_", tdim ) )
            {
                int cell_sub_subdomain;
                this->cell_.get ( "_sub_domain_", &cell_sub_subdomain );
                return cell_sub_subdomain == this->space_->dof ( ).get_my_subdomain ( );
            }
            else
            {
                // assume it is true if we have no other information
                return true;
            }
        }

        /// \brief Evaluate a discrete finite element solution function on the element.
        /// \param[in] var            variable which should be evaluated
        /// \param[in] pt             coordinates of the point at which the solution should be evaluated
        /// \param[in] dof_values     solution values of the dofs on the element
        /// \return the solution evaluated at point @p pt.

        DataType evaluate_fe_solution ( int var,
                                        const std::vector<DataType>& pt,
                                        const std::vector<DataType>& dof_values ) const
        {
            return this->space_->get_solution_value ( var, this->cell_, pt, dof_values );
        }

        /// \return a reference to the finite element space to which the element belongs.

        const VectorSpace<DataType>& space ( ) const
        {
            return *this->space_;
        }

      private:
        const VectorSpace<DataType>* space_;
        mesh::Entity cell_;
    };

} // namespace hiflow

#endif
