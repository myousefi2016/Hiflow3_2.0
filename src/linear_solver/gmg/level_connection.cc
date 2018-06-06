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

#include "level_connection.h"
#include "basic_hierarchy.h"
#include "gmg_level_impl.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            template<class LAD>
            BasicConnection<LAD>::BasicConnection ( BasicLevel<LAD>& coarse_level,
                                                    BasicLevel<LAD>& fine_level )
            {
                info_ = boost::shared_ptr<DataTransferInformation<LAD> >(
                        new DataTransferInformation<LAD>( &coarse_level, &fine_level ) );

                transfer_to_fine_ = boost::shared_ptr<DataTransferCoarseToFine<LAD> >(
                        new DataTransferCoarseToFine<LAD>( *info_ ) );

                transfer_to_coarse_ = boost::shared_ptr<DataTransferFineToCoarse<LAD> >(
                        new DataTransferFineToCoarse<LAD>( *info_ ) );

                dof_ident_ = boost::shared_ptr<DofIdentification<LAD> >(
                        new DofIdentification<LAD>( fine_level, coarse_level, info_, transfer_to_coarse_ ) );

                dof_ident_->identify_dofs ( );
            }

            template class BasicConnection<LADescriptorCoupledD>;
            template class BasicConnection<LADescriptorCoupledS>;

            template class ConnectedLevel<LADescriptorCoupledD, BasicConnection<LADescriptorCoupledD> >;
            template class ConnectedLevel<LADescriptorCoupledS, BasicConnection<LADescriptorCoupledS> >;

            template class BasicHierarchy<ConnectedLevel<LADescriptorCoupledD, BasicConnection<LADescriptorCoupledD> > >;
            template class BasicHierarchy<ConnectedLevel<LADescriptorCoupledS, BasicConnection<LADescriptorCoupledS> > >;

            template<class LAD>
            SolutionTransferringConnection<LAD>::SolutionTransferringConnection ( BasicLevel<LAD>& coarse_level,
                                                                                  BasicLevel<LAD>& fine_level )
            : BasicConnection<LAD>( coarse_level, fine_level )
            {

                solution_transfer_ = boost::shared_ptr<VectorTransfer<LAD> >(
                        new VectorTransfer<LAD>(
                        dof_ident_,
                        info_,
                        fine_level.sol ( ),
                        coarse_level.sol ( ),
                        transfer_to_fine_,
                        transfer_to_coarse_ ) );
            }

            template class SolutionTransferringConnection<LADescriptorCoupledD>;
            template class SolutionTransferringConnection<LADescriptorCoupledS>;

            template class ConnectedLevel<LADescriptorCoupledD, SolutionTransferringConnection<LADescriptorCoupledD> >;
            template class ConnectedLevel<LADescriptorCoupledS, SolutionTransferringConnection<LADescriptorCoupledS> >;

            template class BasicHierarchy<ConnectedLevel<LADescriptorCoupledD, SolutionTransferringConnection<LADescriptorCoupledD> > >;
            template class BasicHierarchy<ConnectedLevel<LADescriptorCoupledS, SolutionTransferringConnection<LADescriptorCoupledS> > >;

            template<class LAD>
            GMGConnection<LAD>::GMGConnection ( LevelType& coarse_level,
                                                LevelType& fine_level )
            : SolutionTransferringConnection<LAD>( coarse_level, fine_level )
            {

                fine_res_to_coarse_rhs_transfer_ = boost::shared_ptr<VectorTransfer<LAD> >(
                        new VectorTransfer<LAD>(
                        dof_ident_,
                        info_,
                        fine_level.res ( ),
                        coarse_level.rhs ( ),
                        transfer_to_fine_,
                        transfer_to_coarse_ ) );

                coarse_sol_to_fine_res_transfer_ = boost::shared_ptr<VectorTransfer<LAD> >(
                        new VectorTransfer<LAD>(
                        dof_ident_,
                        info_,
                        fine_level.res ( ),
                        coarse_level.sol ( ),
                        transfer_to_fine_,
                        transfer_to_coarse_ ) );
            }

            template class GMGConnection<LADescriptorCoupledD>;
            template class GMGConnection<LADescriptorCoupledS>;

            template class GMGLevel<LADescriptorCoupledD, GMGConnection<LADescriptorCoupledD> >;
            template class GMGLevel<LADescriptorCoupledS, GMGConnection<LADescriptorCoupledS> >;

            template class BasicHierarchy<GMGLevel<LADescriptorCoupledD, GMGConnection<LADescriptorCoupledD> > >;
            template class BasicHierarchy<GMGLevel<LADescriptorCoupledS, GMGConnection<LADescriptorCoupledS> > >;

            template void visualize_multilevel_solutions<GMGLevel<LADescriptorCoupledD, GMGConnection<LADescriptorCoupledD> > >( const BasicHierarchy<GMGLevel<LADescriptorCoupledD, GMGConnection<LADescriptorCoupledD> > >&, const std::string& );
            template void visualize_multilevel_solutions<GMGLevel<LADescriptorCoupledS, GMGConnection<LADescriptorCoupledS> > >( const BasicHierarchy<GMGLevel<LADescriptorCoupledS, GMGConnection<LADescriptorCoupledS> > >&, const std::string& );

        } // namespace gmg
    } // namespace la
} // namespace hiflow
