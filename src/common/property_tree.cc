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

/** \file
 * Input/Output in XML format -- non-generic code.
 */

#include "property_tree.h"
#include <cctype>
#include <stack>
#include <iostream>

#define TIXML_USE_STL
#include <tinyxml.h>

/// @author Jonas Latt, Jonas Fietz, Mathias Krause, Tobias Hahn

namespace hiflow
{

    PropertyTree PropertyTree::notFound;

    PropertyTree::PropertyTree ( TiXmlNode* pParent, int master_rank, const MPI_Comm& comm )
    {
        int my_rank = -1;
        MPI_Comm_rank ( comm, &my_rank );
        if ( my_rank == master_rank )
        {
            mainProcessorIni ( pParent, master_rank, comm );
        }
        else
        {
            slaveProcessorIni ( master_rank, comm );
        }
    }

    PropertyTree::PropertyTree ( const std::string& fName, int master_rank, const MPI_Comm& comm )
    {
        TiXmlDocument* doc = 0;
        int loadOK = false;
        int my_rank = -1;
        MPI_Comm_rank ( comm, &my_rank );
        if ( my_rank == master_rank )
        {
            doc = new TiXmlDocument ( fName.c_str ( ) );
            loadOK = doc->LoadFile ( );
            if ( !loadOK )
            {
                std::cout << std::string ( "Problem processing input XML file " ) << fName << std::endl;
            }
        }

        if ( my_rank == master_rank )
        {
            mainProcessorIni ( doc, master_rank, comm );
            delete doc;
        }
        else
        {
            slaveProcessorIni ( master_rank, comm );
        }
    }

    void PropertyTree::mainProcessorIni ( TiXmlNode* pParent, int master_rank, const MPI_Comm& comm )
    {
        assert ( pParent->Type ( ) == TiXmlNode::DOCUMENT || pParent->Type ( ) == TiXmlNode::ELEMENT );

        if ( pParent->Type ( ) == TiXmlNode::DOCUMENT )
        {
            // ignore the surrounding PARAM-block
            pParent = pParent->FirstChild ( );
        }

        name = pParent->ValueStr ( );
        bCastString ( &name, master_rank, comm );

        TiXmlNode * pChild;

        int type = 0;
        for ( pChild = pParent->FirstChild ( ); pChild != 0; pChild = pChild->NextSibling ( ) )
        {
            type = pChild->Type ( );
            MPI_Bcast ( static_cast < void* > ( &type ), 1, MPI_INT, master_rank, comm );
            if ( type == TiXmlNode::ELEMENT )
            {
                PropertyTree* new_child = new PropertyTree ( pChild, master_rank, comm );
                children[new_child->getName ( )] = new_child;
            }
            else if ( type == TiXmlNode::TEXT )
            {
                text = pChild->ToText ( )->ValueStr ( );
                bCastString ( &text, master_rank, comm );
            }
        }
        type = TiXmlNode::UNKNOWN;
        MPI_Bcast ( static_cast < void* > ( &type ), 1, MPI_INT, master_rank, comm );
    }

    void PropertyTree::slaveProcessorIni ( int master_rank, const MPI_Comm& comm )
    {
        bCastString ( &name, master_rank, comm );

        int type = 0;
        do
        {
            MPI_Bcast ( static_cast < void* > ( &type ), 1, MPI_INT, master_rank, comm );
            if ( type == TiXmlNode::ELEMENT )
            {
                PropertyTree* new_child = new PropertyTree ( 0, master_rank, comm );
                children[new_child->getName ( )] = new_child;
            }
            else if ( type == TiXmlNode::TEXT )
            {
                bCastString ( &text, master_rank, comm );
            }
        }
        while ( type != TiXmlNode::UNKNOWN );
    }

    void PropertyTree::bCastString ( std::string* sendBuf, int master_rank, const MPI_Comm& comm )
    {
        int length = ( int ) sendBuf->size ( );
        MPI_Bcast ( static_cast < void* > ( &length ), 1, MPI_INT, master_rank, comm );
        char* buffer = new char[length + 1];
        int rank = -1;
        MPI_Comm_rank ( comm, &rank );
        if ( rank == master_rank )
        {
            std::copy ( sendBuf->c_str ( ), sendBuf->c_str ( ) + length + 1, buffer );
        }
        MPI_Bcast ( static_cast < void* > ( buffer ), length + 1, MPI_CHAR, master_rank, comm );
        if ( rank != master_rank )
        {
            *sendBuf = buffer;
        }
        delete [] buffer;
    }

    PropertyTree::PropertyTree ( )
    {
        name = "XML node not found";
    }

    PropertyTree::PropertyTree ( const PropertyTree& srcTree )
    : name ( srcTree.name ), text ( srcTree.text )
    {
        if ( !srcTree.isEmpty ( ) )
        {
            PropertyTree::const_iterator src_end = srcTree.end_children ( );
            for ( PropertyTree::const_iterator it = srcTree.begin_children ( ); it != src_end; ++it )
                add ( *( it->second ) );
        }
    }

    PropertyTree::~PropertyTree ( )
    {
        if ( !isEmpty ( ) )
        {
            PropertyTree::iterator end_tree = this->end_children ( );
            for ( PropertyTree::iterator it = this->begin_children ( ); it != end_tree; ++it )
                delete it->second;
        }
    }

    void PropertyTree::print ( int indent ) const
    {
        std::string indentStr ( indent, ' ' );
        std::cout << indentStr << "[" << name << "]" << std::endl;
        if ( !text.empty ( ) )
        {
            std::cout << indentStr << "  " << text << std::endl;
        }
        if ( !isEmpty ( ) )
        {
            PropertyTree::const_iterator end_tree = this->end_children ( );
            for ( PropertyTree::const_iterator it = this->begin_children ( ); it != end_tree; ++it )
                it->second->print ( indent + 2 );
        }
    }

    PropertyTree & PropertyTree::operator= ( const PropertyTree & srcTree )
    {
        if ( this != &srcTree )
        {
            if ( !isEmpty ( ) )
            {
                PropertyTree::iterator end_tree = this->end_children ( );
                for ( PropertyTree::iterator it = this->begin_children ( ); it != end_tree; ++it )
                    remove ( it->first );
            }
            name = srcTree.getName ( );
            text = srcTree.getText ( );

            if ( !srcTree.isEmpty ( ) )
            {
                PropertyTree::const_iterator src_end = srcTree.end_children ( );
                for ( PropertyTree::const_iterator it = srcTree.begin_children ( ); it != src_end; ++it )
                    add ( *( it->second ) );
            }
        }

        return *this;
    }

    PropertyTree const& PropertyTree::operator[] ( std::string name ) const
    {
        std::map<std::string, PropertyTree*>::const_iterator element = children.find ( name );
        if ( element != children.end ( ) )
        {
            return *( element->second );
        }
        else
        {
            LOG_ERROR ( "Child " << name << " in tree " << this->name << " not found!" );
            return notFound;
        }
    }

    std::map<std::string, PropertyTree*>::const_iterator PropertyTree::begin_children ( ) const
    {
        return children.begin ( );
    }

    std::map<std::string, PropertyTree*>::iterator PropertyTree::begin_children ( )
    {
        return children.begin ( );
    }

    std::map<std::string, PropertyTree*>::const_iterator PropertyTree::end_children ( ) const
    {
        return children.end ( );
    }

    std::map<std::string, PropertyTree*>::iterator PropertyTree::end_children ( )
    {
        return children.end ( );
    }

    std::string PropertyTree::getName ( ) const
    {
        return name;
    }

    std::string PropertyTree::getText ( ) const
    {
        return text;
    }

    void PropertyTree::add ( const std::string& key )
    {
        // TODO 
    }

    void PropertyTree::add ( const PropertyTree& tree )
    {
        PropertyTree::iterator has_element = children.find ( tree.getName ( ) );
        if ( has_element == children.end ( ) )
        {
            PropertyTree* newTree = new PropertyTree ( tree );
            children[tree.getName ( )] = newTree;
        }
        else
        {
            has_element->second->merge ( tree );
        }
    }

    bool PropertyTree::remove ( const std::string& key )
    {
        std::map<std::string, PropertyTree*>::iterator element = children.find ( key );
        if ( element != children.end ( ) )
        {
            delete element->second;
            children.erase ( element );
            return true;
        }
        else
        {
            return false;
        }
    }

    void PropertyTree::merge ( const PropertyTree& tree )
    {
        if ( name == tree.getName ( ) )
        {
            //text=tree.get<std::string>();
            text = tree.getText ( );
        }
        if ( !tree.isEmpty ( ) )
        {
            for ( PropertyTree::const_iterator it = tree.begin_children ( ), end = tree.end_children ( ); it != end; ++it )
            {
                add ( *( it->second ) );
            }
        }
    }

    bool PropertyTree::contains ( const std::string& key ) const
    {
        std::map<std::string, PropertyTree*>::const_iterator element = children.find ( key );
        if ( element == children.end ( ) )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    bool PropertyTree::isEmpty ( ) const
    {
        switch ( children.size ( ) )
        {
            case 0:
                return true;
            default:
                return false;
        }
    }

    int PropertyTree::size ( ) const
    {
        return children.size ( );
    }

    std::ostream& operator<< ( std::ostream& os, const PropertyTree& tree )
    {
        std::stack< std::pair<std::string, const PropertyTree* > > s;
        s.push ( std::make_pair ( tree.getName ( ), &tree ) );

        while ( !s.empty ( ) )
        {
            const std::pair<std::string, const PropertyTree*> node = s.top ( );
            s.pop ( );

            if ( !node.second->text.empty ( ) )
                os << node.first << " -> " << node.second->get<std::string>( ) << "\n";
            else
                os << node.first << "\n";
            for ( PropertyTree::const_iterator it = node.second->begin_children ( ), end = node.second->end_children ( ); it != end; ++it )
            {
                s.push ( std::make_pair ( node.first + std::string ( "->" ) + it->second->getName ( ), it->second ) );
            }

        }
        return os;
    }

} // namespace hiflow
