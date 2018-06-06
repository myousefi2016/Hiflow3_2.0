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

/// @author Dimitar Lukarski

#include "lvector.h"
#include "lmp_mem.h"
#include "init_vec_mat.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "lmp_log.h"

namespace hiflow
{
    namespace la
    {

        // Class lVector

        template <typename ValueType>
        lVector<ValueType>::lVector ( )
        {
            this->size_ = 0;
            this->name_ = "";
            this->platform_name_ = "";
            this->platform_id_ = CPU;
            this->implementation_name_ = "";
            this->implementation_id_ = NAIVE;
        }

        template <typename ValueType>
        lVector<ValueType>::~lVector ( )
        {
        }

        template <typename ValueType>
        lVector<ValueType> *lVector<ValueType>::CloneWithoutContent ( ) const
        {
            std::string cloned_name;
            cloned_name = "clone from ";
            cloned_name.append ( this->name_ );

            lVector<ValueType> *new_vector = init_vector<ValueType>( this->get_size ( ),
                    cloned_name,
                    this->get_platform ( ),
                    this->get_implementation ( ) );
            new_vector->Zeros ( );
            return new_vector;
        }

        template <typename ValueType>
        void lVector<ValueType>::print ( std::ostream &out ) const
        {

            LOG_INFO ( "lVector",
                       "name='" << this->get_name ( )
                       << "', Elements=" << this->get_size ( )
                       << ", Precision=" << sizeof (ValueType )*8 << "bit"
                       << ", Platform:" << this->get_platform_name ( )
                       << ", Implementation:" << this->get_implementation_name ( )
                       //             << ", L2 norm:" << this->Norm2()
                       );

        }

        template <typename ValueType>
        std::string lVector<ValueType>::get_name ( void ) const
        {
            return this->name_;
        }

        template <typename ValueType>
        enum IMPLEMENTATION lVector<ValueType>::get_implementation ( void ) const
        {
            return implementation_id_;
        }

        template <typename ValueType>
        enum PLATFORM lVector<ValueType>::get_platform ( void ) const
        {
            return platform_id_;
        }

        template <typename ValueType>
        std::string lVector<ValueType>::get_implementation_name ( void ) const
        {
            return this->implementation_name_;
        }

        template <typename ValueType>
        std::string lVector<ValueType>::get_platform_name ( void ) const
        {
            return this->platform_name_;
        }

        template <typename ValueType>
        void lVector<ValueType>::Sync ( void ) const
        {
            // do nothing
        }

        template <typename ValueType>
        void lVector<ValueType>::CopyFromIndexset ( const lVector<ValueType>& vec )
        {

            int *tmp_buff;

            tmp_buff = ( int * ) malloc ( vec.get_indexset_size ( ) * sizeof (int ) );
            assert ( tmp_buff != NULL );

            vec.get_indexset ( tmp_buff );
            this->set_indexset ( tmp_buff, vec.get_indexset_size ( ) );

            free ( tmp_buff );

        }

        template <typename ValueType>
        void lVector<ValueType>::CopyStructureFrom ( const lVector<ValueType>& vec2 )
        {
            if ( this != &vec2 )
            {
                // just init the structure, no data copy
                this->Init ( vec2.get_size ( ), vec2.get_name ( ) );
            }
        }

        template <typename ValueType>
        void lVector<ValueType>::CloneFrom ( const lVector<ValueType>& other )
        {
            if ( this != &other )
            {
                this->Clear ( );
                this->CopyStructureFrom ( other );
                this->CopyFrom ( other );
            }
        }

        template <typename ValueType>
        void lVector<ValueType>::ReadFile ( const char* filename )
        {
            LOG_ERROR ( "lVector::ReadFile() does not support this vector format" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lVector<ValueType>::WriteFile ( const char* filename )
        {
            LOG_ERROR ( "lVector::ReadFile() does not support this vector format" );
            this->print ( );
            exit ( -1 );
        }

        template class lVector<double>;
        template class lVector<float>;

    } // namespace la
} // namespace hiflow
