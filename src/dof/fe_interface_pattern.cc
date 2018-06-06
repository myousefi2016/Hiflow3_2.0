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

#include "dof/fe_interface_pattern.h"
#include "fem/fetype.h"
#include "common/macros.h"

/// \author Michael Schick<br>Martin Baumann

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FEInterfacePattern<DataType>::FEInterfacePattern ( mesh::InterfacePattern interface_pattern,
                                                           FEType<DataType>* fe_type_master,
                                                           std::vector<FEType<DataType>* > fe_type_slaves )
        {
            interface_pattern_ = interface_pattern;
            fe_type_master_ = fe_type_master;
            fe_type_slaves_ = fe_type_slaves;
        }

        template<class DataType>
        FEInterfacePattern<DataType>::~FEInterfacePattern ( )
        {

        }

        template<class DataType>
        int FEInterfacePattern<DataType>::get_interface_degree ( int* which_slave ) const
        {
            int min_degree = fe_type_master_->get_fe_deg ( );

            if ( which_slave )
            {
                *which_slave = -1;
            }

            for ( size_t i = 0, end = num_slaves ( ); i != end; ++i )
            {
                const int deg = fe_type_slaves_[i]->get_fe_deg ( );
                if ( deg < min_degree )
                {
                    min_degree = deg;
                    if ( which_slave )
                    {
                        *which_slave = i;
                    }
                }
            }

            return min_degree;
        }

        template<class DataType>
        bool FEInterfacePattern<DataType>::operator== ( const FEInterfacePattern<DataType>& test )
        {
            if ( interface_pattern_ == test.interface_pattern_ &&
                 fe_type_master_ == test.fe_type_master_ &&
                 fe_type_slaves_ == test.fe_type_slaves_ )
                return true;
            return false;
        }

        /// first check InterfacePattern, as FEInterfacePattern is a specialization

        template<class DataType>
        bool FEInterfacePattern<DataType>::operator< ( const FEInterfacePattern<DataType>& test ) const
        {
            if ( interface_pattern_ < test.interface_pattern_ )
            {
                return true;
            }
            else if ( interface_pattern_ == test.interface_pattern_ )
            {
                if ( fe_type_master_ < test.fe_type_master_ )
                {
                    return true;
                }
                else if ( fe_type_master_ == test.fe_type_master_ )
                {
                    return fe_type_slaves_ < test.fe_type_slaves_;
                }
            }
            return false;
        }

        std::ostream& operator<< ( std::ostream &s, const FEInterfacePattern<double>& pattern )
        {
            s << pattern.interface_pattern ( );
            s << "Ansatz Master: " << pattern.fe_type_master ( ).get_name ( ) << std::endl;
            if ( pattern.num_slaves ( ) == 1 )
                s << "Ansatz Slave:  " << pattern.fe_type_slaves ( )[0]->get_name ( ) << std::endl;
            else
            {
                for ( size_t i = 0, e_i = pattern.num_slaves ( ); i != e_i; ++i )
                    s << "Ansatz Slave " << i << ": " << pattern.fe_type_slaves ( )[i]->get_name ( ) << std::endl;
            }
            return s;
        }

        std::ostream& operator<< ( std::ostream &s, const FEInterfacePattern<float>& pattern )
        {
            s << pattern.interface_pattern ( );
            s << "Ansatz Master: " << pattern.fe_type_master ( ).get_name ( ) << std::endl;
            if ( pattern.num_slaves ( ) == 1 )
                s << "Ansatz Slave:  " << pattern.fe_type_slaves ( )[0]->get_name ( ) << std::endl;
            else
            {
                for ( size_t i = 0, e_i = pattern.num_slaves ( ); i != e_i; ++i )
                    s << "Ansatz Slave " << i << ": " << pattern.fe_type_slaves ( )[i]->get_name ( ) << std::endl;
            }
            return s;
        }

        // template instantiation
        template class FEInterfacePattern<double>;
        template class FEInterfacePattern<float>;

    }
} // namespace hiflow
