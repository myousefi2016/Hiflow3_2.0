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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef LEVEL_CONNECTION_H
#    define LEVEL_CONNECTION_H

#    include "vector_transfer.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Manages all objects that exist between two levels in the Hierarchy:
            /// Sets up DataTransferInformation, data transfer, DofIdentification.
            /// The BasicConnection does not setup any transfer objects.

            template<class LAD>
            class BasicConnection
            {
              public:

                /// Initializes the object. Collective on the communicator of the finer
                /// grid.
                /// @param coarse_level A coarser grid in the MultiLevelHierarchy
                /// @param fine_level A finer grid in the MultiLevelHierarchy.
                /// The communicator of the coarser grid must contain a subset of the processes
                /// in the communicator of the finer grid.
                BasicConnection ( BasicLevel<LAD>& coarse_level,
                                  BasicLevel<LAD>& fine_level );

                virtual ~BasicConnection ( )
                {
                }

                /// @return The DoF identification object

                boost::shared_ptr<const DofIdentification<LAD> > get_dof_identification ( ) const
                {
                    return dof_ident_;
                }

              protected:
                boost::shared_ptr<DataTransferInformation<LAD> > info_;
                boost::shared_ptr<DofIdentification<LAD> > dof_ident_;

                boost::shared_ptr<DataTransferCoarseToFine<LAD> > transfer_to_fine_;
                boost::shared_ptr<DataTransferFineToCoarse<LAD> > transfer_to_coarse_;
            };

            /// An connection that allows the transfers of solution vectors

            template<class LAD>
            class SolutionTransferringConnection : public BasicConnection<LAD>
            {
              public:
                using BasicConnection<LAD>::dof_ident_;
                using BasicConnection<LAD>::info_;
                using BasicConnection<LAD>::transfer_to_coarse_;
                using BasicConnection<LAD>::transfer_to_fine_;

                /// Initializes the object. Collective on the communicator of the finer
                /// grid.
                /// @param coarse_level A coarser grid in the MultiLevelHierarchy
                /// @param fine_level A finer grid in the MultiLevelHierarchy.
                /// The communicator of the coarser grid must contain a subset of the processes
                /// in the communicator of the finer grid.

                SolutionTransferringConnection ( BasicLevel<LAD>& coarse_level,
                                                 BasicLevel<LAD>& fine_level );

                virtual ~SolutionTransferringConnection ( )
                {
                }

                /// @return A pointer to the vector transfer object that enables
                /// the transfer of the content between a vector on the coarse grid and
                /// a vector on the fine grid.

                boost::shared_ptr<VectorTransfer<LAD> > get_solution_transfer ( ) const
                {
                    return solution_transfer_;
                }

              private:
                boost::shared_ptr<VectorTransfer<LAD> > solution_transfer_;

            };

            /// An connection that allows the transfers of solution vectors

            template<class LAD>
            class GMGConnection : public SolutionTransferringConnection<LAD>
            {
              public:
                typedef GMGLevel<LAD, GMGConnection<LAD> > LevelType;

                using BasicConnection<LAD>::dof_ident_;
                using BasicConnection<LAD>::info_;
                using BasicConnection<LAD>::transfer_to_coarse_;
                using BasicConnection<LAD>::transfer_to_fine_;

                /// Initializes the object. Collective on the communicator of the finer
                /// grid.
                /// @param coarse_level A coarser grid in the MultiLevelHierarchy
                /// @param fine_level A finer grid in the MultiLevelHierarchy.
                /// The communicator of the coarser grid must contain a subset of the processes
                /// in the communicator of the finer grid.

                GMGConnection ( LevelType& coarse_level,
                                LevelType& fine_level );

                virtual ~GMGConnection ( )
                {
                }

                /// Transfers the residual of the fine grid to the right hand side
                /// of the coarse grid

                void transfer_fine_res_to_coarse_rhs ( void )
                {
                    assert ( fine_res_to_coarse_rhs_transfer_ != NULL );
                    fine_res_to_coarse_rhs_transfer_->transfer_to_coarse ( );
                }

                /// Transfers the solution of the coarse grid to the residual of the fine grid.

                void transfer_coarse_sol_to_fine_res ( void )
                {
                    assert ( coarse_sol_to_fine_res_transfer_ != NULL );
                    coarse_sol_to_fine_res_transfer_->transfer_to_fine ( );
                }

                /// Transfers the solution of the coarse grid to the solution vector on
                /// the fine grid.

                void transfer_coarse_sol_to_fine_sol ( )
                {
                    assert ( this->get_solution_transfer ( ) != NULL );
                    this->get_solution_transfer ( )->transfer_to_fine ( );
                }

              private:
                boost::shared_ptr<VectorTransfer<LAD> > fine_res_to_coarse_rhs_transfer_;
                boost::shared_ptr<VectorTransfer<LAD> > coarse_sol_to_fine_res_transfer_;

            };

        }
    }
}

#endif
