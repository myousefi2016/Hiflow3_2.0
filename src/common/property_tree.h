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

#ifndef _PROPERTY_TREE_H_
#    define _PROPERTY_TREE_H_

#    include <string>
#    include <vector>
#    include <iostream>
#    include <sstream>
#    include <algorithm>
#    include <mpi.h>
#    include "log.h"
#    include "macros.h"

class TiXmlNode;

namespace hiflow
{

    /// @brief Describes a tree of properties, i.e a hierarchical
    /// collection of key-value pairs. The keys are of type
    /// std::string, and the values of the template parameter type
    /// ValueType. ValueType must be default-constructible.
    /// @author Jonas Latt, Jonas Fietz, Mathias Krause, Tobias Hahn

    class PropertyTree
    {
      public:

        typedef std::map<std::string, PropertyTree*>::const_iterator const_iterator;
        typedef std::map<std::string, PropertyTree*>::iterator iterator;

        /// Constructs a new PropertyTree from an XML file
        PropertyTree ( const std::string& fName, int master_rank, const MPI_Comm& comm );

        /// Copy constructor
        PropertyTree ( const PropertyTree& srcTree );

        /// Destructor
        ~PropertyTree ( );

        /// Assignment
        PropertyTree & operator= ( const PropertyTree & srcTree );

        /// Return the subtree with given name
        PropertyTree const& operator[] ( std::string name ) const;

        /// Standard read function of node value.
        template <typename T> bool read ( T& value ) const;
        template <typename T> bool read ( std::vector<T>& value ) const;

        /// Direct retrieval of node value.
        template <typename T> T get ( ) const;
        /// Direct retrieval of node value with the possibility to set an default value if the value was not found.    
        template <typename T> T get ( T def ) const;

        /// Merges the tree with another, preferring the other values in case of same keys.
        void merge ( const PropertyTree& tree );

        /// Returns true if the tree contains child with specified key
        bool contains ( const std::string& key ) const;

        /// Returns true if this tree is empty, i.e. the tree is a leaf.
        bool isEmpty ( ) const;

        /// Returns the number of keys in this level.
        int size ( ) const;

        friend std::ostream& operator<< ( std::ostream& os, const PropertyTree& tree );

        /// Prints out the XML structure read in, mostly for debugging purposes
        void print ( int indent ) const;

        /// Return an iterator for this level at the tree
        const_iterator begin_children ( ) const;
        iterator begin_children ( );

        /// Return an iterator end element 
        const_iterator end_children ( ) const;
        iterator end_children ( );

        /// Return the name of the element
        std::string getName ( ) const;
        /// Return the text of the element
        std::string getText ( ) const;

      private:
        PropertyTree ( );

        /// Construct a XML node.
        PropertyTree ( TiXmlNode* pParent, int master_rank, const MPI_Comm& comm );

        /// Construct a new tree node with given name and value.
        PropertyTree ( const std::string& name, const std::string& text, int master_rank, const MPI_Comm& comm );

        /// Set value of node.
        template <typename T> void set ( const T& value );

        /// Add a child node.
        void add ( const std::string& key );
        template <typename T> void add ( const std::string& key, const T& value );

        /// Add a subtree
        void add ( const PropertyTree& tree );

        /// Remove a subtree.
        bool remove ( const std::string& key );

        void mainProcessorIni ( TiXmlNode* pParent, int master_rank, const MPI_Comm& comm );
        void slaveProcessorIni ( int master_rank, const MPI_Comm& comm );
        void bCastString ( std::string* sendbuf, int master_rank, const MPI_Comm& comm );
      private:
        std::string name;
        std::string text;
        std::map<std::string, PropertyTree*> children;
        static PropertyTree notFound;
    };

    template <typename T>
    bool PropertyTree::read ( T& value ) const
    {
        std::stringstream valueStr ( text );
        T tmp = T ( );
        if ( !( valueStr >> tmp ) )
        {
            if ( this->isEmpty ( ) ) //Only log errors for leafs of the tree 
                LOG_ERROR ( "Cannot read value from property element " << name );
            return false;
        }
        value = tmp;
        return true;
    }

    template <>
    inline bool PropertyTree::read ( bool& value ) const
    {
        std::stringstream valueStr ( text );
        std::string word;
        valueStr >> word;
        // Transform to lower-case, so that "true" and "false" are case-insensitive.
        std::transform ( word.begin ( ), word.end ( ), word.begin ( ), ::tolower );
        if ( word == "true" )
        {
            value = true;
            return true;
        }
        else if ( word == "false" )
        {
            value = false;
            return true;
        }
        else
        {
            LOG_ERROR ( "Cannot read boolean value from XML element " << name );
        }
        return false;
    }

    template <typename T>
    inline bool PropertyTree::read ( std::vector<T>& values ) const
    {
        std::stringstream multiValueStr ( text );
        std::string word;
        std::vector<T> tmp ( values );
        while ( multiValueStr >> word )
        {
            std::stringstream valueStr ( word );
            T value;
            if ( !( valueStr >> value ) )
            {
                if ( this->isEmpty ( ) ) //Only log errors for leafs of the tree 
                    LOG_ERROR ( "Cannot read value array from property element. " << name );
                return false;
            }
            tmp.push_back ( value );
        }
        values.swap ( tmp );
        return true;
    }

    template <typename T>
    inline T PropertyTree::get ( ) const
    {
        std::stringstream valueStr ( text );
        T tmp = T ( );
        if ( !( valueStr >> tmp ) )
        {
            LOG_ERROR ( "Cannot read value from property element. " << name << "." );
            exit ( -1 );
        }
        return tmp;
    }

    template <typename T>
    inline T PropertyTree::get ( T def ) const
    {
        std::stringstream valueStr ( text );
        T tmp = T ( );
        if ( !( valueStr >> tmp ) )
        {
            LOG_INFO ( "Property tree element not found. Using default", def );
            return def;
        }
        return tmp;
    }

    template <>
    inline bool PropertyTree::read ( std::string& entry ) const
    {
        if ( name == "XML node not found" )
        {
            return false;
        }
        entry = text;
        return true;
    }

    template <typename T>
    inline void PropertyTree::set ( const T &value )
    {
        std::stringstream valueStr;
        valueStr << value;
        text = valueStr.str ( );
    }

} // namespace hiflow

#endif  // _PROPERTY_TREE_H_
