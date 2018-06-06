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

#ifndef HIFLOW_VECTOR_SPACE_H_
#    define HIFLOW_VECTOR_SPACE_H_

#    include <vector>

#    include "mesh/entity.h"
#    include "fem/cell_transformation.h"
#    include "dof/degree_of_freedom.h"
#    include "dof/dof_partition.h"
#    include "fem/femanager.h"
#    include "mesh/iterator.h"
#    include "mesh/mesh.h"
#    include "mesh/types.h"

namespace hiflow
{

    /// @brief A HiFlow vector space.
    /// @author Martin Baumann, Chandramowli Subramanian, Michael Schick
    ///
    /// Contains FEM and DOF data.

    template<class DataType>
    class VectorSpace
    {
      public:
        typedef mesh::Mesh Mesh;
        typedef mesh::Entity MeshEntity;
        typedef mesh::EntityIterator MeshEntityIterator;
        typedef mesh::EntityNumber MeshEntityNumber;
        typedef mesh::Coordinate Coordinate;

        typedef doffem::DegreeOfFreedom<DataType> DegreeOfFreedom;
        typedef doffem::DofPartition<DataType> DofPartition;
        typedef doffem::FEManager<DataType> FEManager;
        typedef doffem::FEType<DataType> FEType;
        typedef doffem::CellTransformation<DataType> CellTransformation;

        VectorSpace ( );
        VectorSpace ( const MPI_Comm& comm );
        virtual ~VectorSpace ( );

        inline void GetDofIndicesFacet ( int var, const MeshEntity& cell,
                                         int tdim, int sindex, std::vector<int>* inddof ) const;
        inline void GetDofIndicesFacet ( const MeshEntity& cell,
                                         int tdim, int sindex, std::vector<int>* inddof ) const;
        inline void GetDofIndices ( int var, const MeshEntity& cell,
                                    std::vector<int>* inddof ) const;
        inline void GetDofIndices ( const MeshEntity& cell,
                                    std::vector<int>* inddof ) const;
        inline int GetNumberDofsOnFacet ( int tdim, const MeshEntity& cell ) const;
        inline int GetNumberDofsOnCell ( const MeshEntity& cell ) const;
        inline void GetCoordOfCell ( const MeshEntity& cell,
                                     std::vector<DataType>* coord ) const;

        inline const CellTransformation& GetCellTransformation
        ( const MeshEntity& cell ) const;

        inline const MPI_Comm& get_mpi_comm ( ) const
        {
            return comm_;
        }
        /// @return dof

        const DofPartition& dof ( ) const
        {
            return *( this->dof_ );
        }

        DofPartition& dof ( )
        {
            return *( this->dof_ );
        }

        /// @return fe manager

        const FEManager& fe_manager ( ) const
        {
            return *( this->fe_manager_ );
        }

        FEManager& fe_manager ( )
        {
            return *( this->fe_manager_ );
        }

        /// @return mesh

        const Mesh& mesh ( ) const
        {
            return *( this->mesh_ );
        }

        mesh::MeshPtr meshPtr ( ) const
        {
            return this->mesh_;
        }

        /// Get mapping of dofs to finite element function for variables defined in
        /// vars. Output vector is ordered by global indices of the specific dofs.
        /// This functionality is needed/useful if a system of PDE is solved and the
        /// BoomerAMG preconditioner of the hypre software package is used.
        /// @param[in] vars Variables/solution function of which the mapping global-
        /// dof-id -> var is needed
        std::vector<int> get_dof_func ( std::vector<int> vars = std::vector<int>( 0 ) ) const;

        /// Initialize scalar space (DoF, Fem, ...)
        void Init ( int degree, Mesh& mesh, hiflow::doffem::DOF_ORDERING order = hiflow::doffem::HIFLOW_CLASSIC );

        /// Initialize vector valued space (DoF, Fem, ...)
        void Init ( std::vector<int> degree, Mesh& mesh, std::vector<bool> is_cg = std::vector<bool>( 0 ), hiflow::doffem::DOF_ORDERING order = hiflow::doffem::HIFLOW_CLASSIC );

        /// Initialize p-refinement ready vector space (DoF, Fem, ...)
        void Init_p ( std::vector<std::vector<int> > degree, Mesh& mesh, hiflow::doffem::DOF_ORDERING order = hiflow::doffem::HIFLOW_CLASSIC );

        /// Clear space (DoF, Fem, ...)
        void Clear ( );

        /// @return dimension

        int get_dim ( ) const
        {
            return this->fe_manager_->get_dim ( );
        }

        /// @return number of variables

        int get_nb_var ( ) const
        {
            return this->fe_manager_->get_nb_var ( );
        }

        /// \brief Get the value for a specific variable (solution) for
        ///        given coordinates and and Cell, which contains the
        ///        coordinates
        template<class VectorType>
        DataType get_solution_value ( int var, const MeshEntity& cell,
                                      const std::vector<DataType>& coord,
                                      const VectorType& sol ) const;

        /// Print status information
        void Print ( ) const;

      private:
        // no implementation of copy constructor or assignement operator
        VectorSpace ( const VectorSpace& );
        VectorSpace operator= ( const VectorSpace& );

        // DegreeOfFreedom*  dof_;
        DofPartition* dof_;
        FEManager* fe_manager_;

        Mesh* mesh_; // only a reference
        /// MPI Communicator
        MPI_Comm comm_;
    };

    /// INLINE FUNCTIONS

    /// Get dof indices on the boundary of a cell of one variable.
    /// @param var variable number
    /// @param cell cell
    /// @param tdim dimension of subentity
    /// @param sindex number of subentity
    /// @param inddof dof indices

    template<class DataType>
    inline void VectorSpace<DataType>::GetDofIndicesFacet ( int var, const MeshEntity& cell,
                                                            int tdim, int sindex, std::vector<int>* inddof ) const
    {
        this->dof_->get_dofs_on_subentity ( var, cell.index ( ), tdim, sindex, ( *inddof ) );
    }

    /// Get all boundary dof indices on a cell.
    /// @param cell cell
    /// @param inddof dof indices

    template<class DataType>
    inline void VectorSpace<DataType>::GetDofIndicesFacet ( const MeshEntity& cell,
                                                            int tdim, int sindex, std::vector<int>* inddof ) const
    {
        // TODO high level function in degree_of_freedom
        // now doing this manually

        inddof->resize ( this->GetNumberDofsOnFacet ( tdim, cell ) );
        std::vector<int> inddof_var;

        // loop on every var
        std::vector<int>::iterator iter = inddof->begin ( );

        for ( int var = 0; var<this->get_nb_var ( ); ++var )
        {
            this->GetDofIndicesFacet ( var, cell, tdim, sindex, &inddof_var );
            // copy values to inddof
            copy ( inddof_var.begin ( ), inddof_var.end ( ), iter );
            iter += inddof_var.size ( );
        }

        assert ( iter == inddof->end ( ) );
    }

    /// Get dof indices on a cell of one variable.
    /// @param var variable number
    /// @param cell cell
    /// @param inddof dof indices

    template<class DataType>
    inline void VectorSpace<DataType>::GetDofIndices ( int var, const MeshEntity& cell,
                                                       std::vector<int>* inddof ) const
    {
        this->dof_->get_dofs_on_cell ( var, cell.index ( ), ( *inddof ) );
    }

    /// Get all dof indices on a cell.
    /// @param cell cell
    /// @param inddof dof indices

    template<class DataType>
    inline void VectorSpace<DataType>::GetDofIndices ( const MeshEntity& cell,
                                                       std::vector<int>* inddof ) const
    {
        // TODO high level function in degree_of_freedom
        // now doing this manually

        inddof->resize ( this->GetNumberDofsOnCell ( cell ) );
        std::vector<int> inddof_var;

        // loop on every var
        std::vector<int>::iterator iter = inddof->begin ( );

        for ( int var = 0; var<this->get_nb_var ( ); ++var )
        {
            this->GetDofIndices ( var, cell, &inddof_var );
            // copy values to inddof
            copy ( inddof_var.begin ( ), inddof_var.end ( ), iter );
            iter += inddof_var.size ( );
        }

        assert ( iter == inddof->end ( ) );
    }

    /// @return number of DOFs on facet of cell
    /// @param cell cell

    template<class DataType>
    inline int VectorSpace<DataType>::GetNumberDofsOnFacet ( int tdim, const MeshEntity& cell ) const
    {
        return this->dof_->get_nb_dofs_on_subentity ( tdim, cell.index ( ) );
    }

    /// @return number of DOFs on @em cell
    /// @param cell cell

    template<class DataType>
    inline int VectorSpace<DataType>::GetNumberDofsOnCell ( const MeshEntity& cell ) const
    {
        return this->dof_->get_nb_dofs_on_cell ( cell.index ( ) );
    }

    /// Get coordinates of a cell w.r.t. dof numbering.
    /// @param cell cell
    /// @param coord coordinates of this cell stored in a vector of vectors

    template<class DataType>
    void VectorSpace<DataType>::GetCoordOfCell ( const MeshEntity& cell,
                                                 std::vector<DataType>* coord ) const
    {
        coord->clear ( );
        cell.get_coordinates<DataType>( *coord );
    }

    /// @return cell transformation of @em cell
    ///         only dependent on mesh

    template<class DataType>
    inline const doffem::CellTransformation<DataType>& VectorSpace<DataType>::GetCellTransformation
    ( const MeshEntity& cell ) const
    {
        return *( this->fe_manager ( ).get_cell_transformation ( cell.index ( ) ) );
    }

} // namespace hiflow

#endif  // HIFLOW_VECTOR_SPACE_H_
