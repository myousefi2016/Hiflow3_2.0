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

#ifndef __FEM_FECELLWISEMANAGER_H_
#    define __FEM_FECELLWISEMANAGER_H_

#    include "fem/femanager.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class FECellwiseManager fecellwisemanager.h
        /// \brief Derived class of FEManager which is needed in case of p refinement
        /// \author Michael Schick<br>Martin Baumann
        ///

        template<class DataType>
        class FECellwiseManager : public FEManager<DataType>
        {
          public:

            typedef std::vector<DataType> Coord;

            /// Use this constructor with dimension and number of variables
            explicit FECellwiseManager ( int dim, int nb_var );

            /// Default destructor

            virtual ~FECellwiseManager ( )
            {
            }

            /// \brief Initialize FE Tank in p-refinement ready version
            virtual void init_fe_tank ( int var,
                                        const typename FEType<DataType>::FEAnsatz &ansatz,
                                        const std::vector<int>& param );

            /// After once (and only once!) using init_fe_tank, use this for p refinements TODO
            /// But this function is not implemented yet and the interface needs to be overthought
            void reinit_fe_tank ( const std::vector<std::vector<int> >& fe_data );
        };

    }
} // namespace hiflow

#endif
