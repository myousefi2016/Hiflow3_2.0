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

#ifndef _DOF_FE_INTERFACE_PATTERN_H_
#    define _DOF_FE_INTERFACE_PATTERN_H_

#    include "mesh/interface.h"

#    include <iostream>

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        class FEType;

        /// \author Michael Schick<br>Martin Baumann

        template<class DataType>
        class FEInterfacePattern
        {
          public:

            FEInterfacePattern ( mesh::InterfacePattern, FEType<DataType>*, std::vector<FEType<DataType>* > );
            ~FEInterfacePattern ( );

            int num_slaves ( ) const
            {
                return interface_pattern_.num_slaves ( );
            }

            mesh::InterfacePattern const& interface_pattern ( ) const
            {
                return interface_pattern_;
            }

            //   mesh::EntityNumber&       master_representer()       { return master_representer_; }
            //   mesh::EntityNumber const& master_representer() const { return master_representer_; }
            // 
            //   std::vector<mesh::EntityNumber>&       slave_representer()       { return slave_representer_; }
            //   std::vector<mesh::EntityNumber> const& slave_representer() const { return slave_representer_; }
            //   mesh::EntityNumber const& slave_representer(int i) const { return slave_representer_[i]; }

            FEType<DataType>& fe_type_master ( )
            {
                return *fe_type_master_;
            }

            FEType<DataType> const& fe_type_master ( ) const
            {
                return *fe_type_master_;
            }

            std::vector<FEType<DataType>* >& fe_type_slaves ( )
            {
                return fe_type_slaves_;
            }

            std::vector<FEType<DataType>* > const& fe_type_slaves ( ) const
            {
                return fe_type_slaves_;
            }

            int master_facet_number ( ) const
            {
                return interface_pattern_.master_facet_number ( );
            }

            int slave_facet_number ( int i ) const
            {
                return interface_pattern_.slave_facet_number ( i );
            }

            /// \brief Compute FE degree of interface, and find FEType of cell
            /// with this degree.
            int get_interface_degree ( int* which_slave ) const;

            bool operator== ( const FEInterfacePattern<DataType>& );
            bool operator< ( const FEInterfacePattern<DataType>& ) const;

            /// overloaded out stream operator
            friend std::ostream& operator<< ( std::ostream&, const FEInterfacePattern<double>& );
            friend std::ostream& operator<< ( std::ostream&, const FEInterfacePattern<float>& );

          private:

            /// \brief represents the geometrical description of the interface
            mesh::InterfacePattern interface_pattern_;

            /// \brief finite element ansatz of master cell
            FEType<DataType>* fe_type_master_;

            /// \brief finite element ansatz of each slave cell
            std::vector<FEType<DataType>* > fe_type_slaves_;

            /// \brief holds index of cell A of first occurrence of interface mode, 
            ///        temporally used to calculate interpolation
            //   mesh::EntityNumber              master_representer_;

            /// \brief holds index of cell B of first occurrence of interface mode, 
            ///        temporally used to calculate interpolation
            //   std::vector<mesh::EntityNumber> slave_representer_;

        };

    }
} // namespace hiflow

#endif
