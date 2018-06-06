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

/// \author Thomas Gengenbach and Staffan Ronnas

#ifndef HIFLOW_MESH_ATTRIBUTES_H
#    define HIFLOW_MESH_ATTRIBUTES_H

#    include "mesh/types.h"

#    include <exception>
#    include <map>
#    include <string>

namespace hiflow
{
    namespace mesh
    {

        class Attribute
        {
          public:

            virtual ~Attribute ( )
            {
            }
            virtual bool get_bool_value ( int index ) const = 0;
            virtual int get_int_value ( int index ) const = 0;
            virtual double get_double_value ( int index ) const = 0;

            virtual void set_bool_value ( int index, bool value ) = 0;
            virtual void set_int_value ( int index, int value ) = 0;
            virtual void set_double_value ( int index, double value ) = 0;

            virtual void set_all_to_zero ( ) = 0;
            virtual void set_all_to_one ( ) = 0;

            virtual int size ( ) = 0;
        };

        typedef SharedPtr<Attribute>::Type AttributePtr;

        struct AttributeTypeException : public std::exception
        {
            virtual const char* what ( ) const throw ( );
        };

        struct ImmutableAttributeException : public std::exception
        {
            virtual const char* what ( ) const throw ( );
        };

        class AttributeBase : public Attribute
        {
          public:
            virtual bool get_bool_value ( int index ) const;
            virtual int get_int_value ( int index ) const;
            virtual double get_double_value ( int index ) const;
            virtual void set_bool_value ( int index, bool value );
            virtual void set_int_value ( int index, int value );
            virtual void set_double_value ( int index, double value );
            virtual int size ( );

            virtual void set_all_to_zero ( );
            virtual void set_all_to_one ( );
        };

        class BoolAttribute : public AttributeBase
        {
          public:
            BoolAttribute ( );
            BoolAttribute ( const std::vector<bool>& data );

            virtual bool get_bool_value ( int index ) const;
            virtual void set_bool_value ( int index, bool value );

            virtual void set_all_to_zero ( );
            virtual void set_all_to_one ( );

            virtual int size ( );

            void resize ( int new_size, bool value = false );

          private:
            std::vector<bool> data_;
        };

        class IntAttribute : public AttributeBase
        {
          public:
            IntAttribute ( );
            IntAttribute ( const std::vector<int>& data );

            virtual int get_int_value ( int index ) const;
            virtual void set_int_value ( int index, int value );

            virtual void set_all_to_zero ( );
            virtual void set_all_to_one ( );

            virtual int size ( );

            void resize ( int new_size, int value = 0 );

          private:
            std::vector<int> data_;
        };

        class DoubleAttribute : public AttributeBase
        {
          public:
            DoubleAttribute ( );
            DoubleAttribute ( const std::vector<double>& data );

            virtual double get_double_value ( int index ) const;
            virtual void set_double_value ( int index, double value );

            virtual void set_all_to_zero ( );
            virtual void set_all_to_one ( );

            virtual int size ( );

            void resize ( int new_size, double value = 0.0 );

          private:
            std::vector<double> data_;
        };

        class InheritedAttribute : public Attribute
        {
          public:
            InheritedAttribute ( );
            InheritedAttribute ( AttributePtr data, AttributePtr parent_index );

            virtual bool get_bool_value ( int index ) const;
            virtual int get_int_value ( int index ) const;
            virtual double get_double_value ( int index ) const;

            virtual void set_bool_value ( int index, bool value );
            virtual void set_int_value ( int index, int value );
            virtual void set_double_value ( int index, double value );

            virtual void set_all_to_zero ( );
            virtual void set_all_to_one ( );

            virtual int size ( );

            AttributePtr data_attribute ( ) const
            {
                return data_;
            }

          private:
            int get_parent_index ( int index ) const;

            /// \brief Attribute with inherited data
            AttributePtr data_;

            /// \brief Attribute for parent index
            AttributePtr parent_index_;

        };

        struct MissingAttributeException : public std::exception
        {
            virtual const char* what ( ) const throw ( );
        };

        class AttributeTable
        {
          public:

            // Need to use output parameter instead of return value since
            // template member functions cannot be specialized.
            template<typename T>
            void get ( const std::string& name, int index, T* value ) const;
            void get ( const std::string& name, int index, int* value ) const;
            void get ( const std::string& name, int index, double* value ) const;

            template<typename T>
            void set ( const std::string& name, int index, const T& value );
            void set ( const std::string& name, int index, int value );
            void set ( const std::string& name, int index, double value );

            void add_attribute ( const std::string& name, AttributePtr attribute );
            AttributePtr get_attribute ( const std::string& name ) const;
            bool has_attribute ( const std::string& name ) const;

            std::vector<std::string> get_attribute_names ( ) const;

          private:

            std::map< std::string, AttributePtr > attributes_;

        };

        //////////////// AttributeTable::get() ////////////////

        template<typename T>
        void AttributeTable::get ( const std::string& name, int index, T* value ) const
        {
            assert ( false );
        }

        //////////////// AttributeTable::set() ////////////////

        template<typename T>
        void AttributeTable::set ( const std::string& name, int index, const T& value )
        {
            assert ( false );
        }

    }
} // namespace hiflow

#endif /* _ATTRIBUTES_H_ */
