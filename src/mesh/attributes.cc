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

#include "attributes.h"

namespace hiflow
{
    namespace mesh
    {

        //////////////// AttributeTypeException ////////////////

        const char* AttributeTypeException::what ( ) const throw ( )
        {
            // TODO (Staffan) better error message
            return "Erroneous attribute type";
        }

        //////////////// ImmutableAttributeException ////////////////

        const char* ImmutableAttributeException::what ( ) const throw ( )
        {
            return "Immutable attribute";
        }

        //////////////// AttributeBase ////////////////

        bool AttributeBase::get_bool_value ( int index ) const
        {
            throw AttributeTypeException ( );
        }

        int AttributeBase::get_int_value ( int index ) const
        {
            throw AttributeTypeException ( );
        }

        double AttributeBase::get_double_value ( int index ) const
        {
            throw AttributeTypeException ( );
        }

        void AttributeBase::set_bool_value ( int index, bool value )
        {
            throw AttributeTypeException ( );
        }

        void AttributeBase::set_int_value ( int index, int value )
        {
            throw AttributeTypeException ( );
        }

        void AttributeBase::set_double_value ( int index, double value )
        {
            throw AttributeTypeException ( );
        }

        void AttributeBase::set_all_to_zero ( )
        {
            throw AttributeTypeException ( );
        }

        void AttributeBase::set_all_to_one ( )
        {
            throw AttributeTypeException ( );
        }

        int AttributeBase::size ( )
        {
            throw AttributeTypeException ( );
        }

        //////////////// BoolAttribute ////////////////

        BoolAttribute::BoolAttribute ( )
        {
        }

        BoolAttribute::BoolAttribute ( const std::vector<bool>& data )
        : data_ ( data )
        {
        }

        bool BoolAttribute::get_bool_value ( int index ) const
        {
            assert ( index >= 0 );
            assert ( index < static_cast < int > ( data_.size ( ) ) );
            return data_[index];
        }

        void BoolAttribute::set_bool_value ( int index, bool value )
        {
            assert ( index >= 0 );
            assert ( index < static_cast < int > ( data_.size ( ) ) );
            data_[index] = value;
        }

        void BoolAttribute::set_all_to_zero ( )
        {
            for ( std::vector<bool>::iterator it = data_.begin ( ), e_it = data_.end ( ); it != e_it; ++it )
            {
                *it = 0;
            }
        }

        void BoolAttribute::set_all_to_one ( )
        {
            for ( std::vector<bool>::iterator it = data_.begin ( ); it != data_.end ( ); ++it )
            {
                *it = 1;
            }
        }

        void BoolAttribute::resize ( int new_size, bool value )
        {
            data_.resize ( new_size, value );
        }

        int BoolAttribute::size ( )
        {
            return data_.size ( );
        }

        //////////////// IntAttribute ////////////////

        IntAttribute::IntAttribute ( )
        {
        }

        IntAttribute::IntAttribute ( const std::vector<int>& data )
        : data_ ( data )
        {
        }

        int IntAttribute::get_int_value ( int index ) const
        {
            assert ( index >= 0 );
            assert ( index < static_cast < int > ( data_.size ( ) ) );
            return data_[index];
        }

        void IntAttribute::set_int_value ( int index, int value )
        {
            assert ( index >= 0 );
            assert ( index < static_cast < int > ( data_.size ( ) ) );
            data_[index] = value;
        }

        void IntAttribute::set_all_to_zero ( )
        {
            for ( std::vector<int>::iterator it = data_.begin ( ), e_it = data_.end ( ); it != e_it; ++it )
            {
                *it = 0;
            }
        }

        void IntAttribute::set_all_to_one ( )
        {
            for ( std::vector<int>::iterator it = data_.begin ( ), e_it = data_.end ( ); it != e_it; ++it )
            {
                *it = 1;
            }
        }

        void IntAttribute::resize ( int new_size, int value )
        {
            data_.resize ( new_size, value );
        }

        int IntAttribute::size ( )
        {
            return data_.size ( );
        }

        //////////////// DoubleAttribute ////////////////

        DoubleAttribute::DoubleAttribute ( )
        {
        }

        DoubleAttribute::DoubleAttribute ( const std::vector<double>& data )
        : data_ ( data )
        {
        }

        double DoubleAttribute::get_double_value ( int index ) const
        {
            assert ( index >= 0 );
            assert ( index < static_cast < int > ( data_.size ( ) ) );
#if 1
            // NaN check
            if ( data_[index] != data_[index] )
            {
                return -0.999E30; // arbitrary choice that is unlikely to turn up otherwise
            }
#endif
            return data_[index];
        }

        void DoubleAttribute::set_double_value ( int index, double value )
        {
            assert ( index >= 0 );
            assert ( index < static_cast < int > ( data_.size ( ) ) );
            data_[index] = value;
        }

        void DoubleAttribute::set_all_to_zero ( )
        {
            for ( std::vector<double>::iterator it = data_.begin ( ), e_it = data_.end ( ); it != e_it; ++it )
            {
                *it = 0;
            }
        }

        void DoubleAttribute::set_all_to_one ( )
        {
            for ( std::vector<double>::iterator it = data_.begin ( ), e_it = data_.end ( ); it != e_it; ++it )
            {
                *it = 1;
            }
        }

        void DoubleAttribute::resize ( int new_size, double value )
        {
            data_.resize ( new_size, value );
        }

        int DoubleAttribute::size ( )
        {
            return data_.size ( );
        }

        //////////////// InheritedAttribute ////////////////

        InheritedAttribute::InheritedAttribute ( )
        {
        }

        InheritedAttribute::InheritedAttribute ( AttributePtr data, AttributePtr parent_index )
        : data_ ( data ), parent_index_ ( parent_index )
        {
        }

        bool InheritedAttribute::get_bool_value ( int index ) const
        {
            return data_->get_bool_value ( get_parent_index ( index ) );
        }

        int InheritedAttribute::get_int_value ( int index ) const
        {
            return data_->get_int_value ( get_parent_index ( index ) );
        }

        double InheritedAttribute::get_double_value ( int index ) const
        {
            return data_->get_double_value ( get_parent_index ( index ) );
        }

        void InheritedAttribute::set_bool_value ( int index, bool value )
        {
            throw ImmutableAttributeException ( );
        }

        void InheritedAttribute::set_int_value ( int index, int value )
        {
            throw ImmutableAttributeException ( );
        }

        void InheritedAttribute::set_double_value ( int index, double value )
        {
            throw ImmutableAttributeException ( );
        }

        int InheritedAttribute::get_parent_index ( int index ) const
        {
            assert ( parent_index_ != 0 );
            return parent_index_->get_int_value ( index );
        }

        void InheritedAttribute::set_all_to_zero ( )
        {
            // do nothing
        }

        void InheritedAttribute::set_all_to_one ( )
        {
            // do nothing
        }

        int InheritedAttribute::size ( )
        {
            return data_->size ( );
        }

        //////////////// MissingAttributeException ////////////////

        const char* MissingAttributeException::what ( ) const throw ( )
        {
            // TODO (Staffan) better message
            return "Missing attribute";
        }

        //////////////// AttributeTable ////////////////

        void AttributeTable::add_attribute ( const std::string& name, AttributePtr attribute )
        {
            if ( attribute != 0 )
            {
                attributes_[name] = attribute;
            }
        }

        AttributePtr AttributeTable::get_attribute ( const std::string& name ) const
        {
            std::map< std::string, AttributePtr >::const_iterator it = attributes_.find ( name );
            if ( it == attributes_.end ( ) )
            {
                std::cerr << "Attribute " << name << " is missing!\n";
                throw MissingAttributeException ( );
            }
            return it->second;
        }

        bool AttributeTable::has_attribute ( const std::string& name ) const
        {
            return attributes_.find ( name ) != attributes_.end ( );
        }

        std::vector<std::string> AttributeTable::get_attribute_names ( ) const
        {
            std::vector<std::string> names;
            for ( std::map<std::string, AttributePtr>::const_iterator it = attributes_.begin ( ); it != attributes_.end ( ); ++it )
            {
                names.push_back ( it->first );
            }
            return names;
        }

        void AttributeTable::get ( const std::string& name, int index, int* value ) const
        {
            assert ( value != 0 );
            *value = get_attribute ( name )->get_int_value ( index );
        }

        void AttributeTable::get ( const std::string& name, int index, double* value ) const
        {
            assert ( value != 0 );
            *value = get_attribute ( name )->get_double_value ( index );
        }

        void AttributeTable::set ( const std::string& name, int index, int value )
        {
            return get_attribute ( name )->set_int_value ( index, value );
        }

        void AttributeTable::set ( const std::string& name, int index, double value )
        {
            return get_attribute ( name )->set_double_value ( index, value );
        }

    }
} // namespace hiflow
