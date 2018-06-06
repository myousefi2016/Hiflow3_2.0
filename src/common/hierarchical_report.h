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

#ifndef HIFLOW_COMMON_HIERARCHICAL_REPORT_H
#    define HIFLOW_COMMON_HIERARCHICAL_REPORT_H

/// \author Staffan Ronnas

namespace hiflow
{

    template<class Data>
    class HierarchicalReport
    {
      public:
        inline HierarchicalReport ( );
        virtual inline ~HierarchicalReport ( );

        // Start a new section.
        inline Data* begin_section ( const std::string& name, const Data& data = Data ( ) );

        // End a section.
        inline Data* end_section ( );

        // Traverse sections depth-first with callback functor
        // visitor. Visitor must define functions enter and exit both
        // taking as parameters the name and a pointer to the
        // section's data. For each section, enter() will be called
        // before its children are visited, and exit() will be called
        // after they have been visited.
        template<class Visitor>
        inline void traverse_depth_first ( Visitor& visitor ) const;

      private:
        typedef int SectionId;

        template<class Visitor>
        void recursive_visit ( Visitor& visitor, SectionId id ) const;

        HierarchicalReport ( const HierarchicalReport& );
        HierarchicalReport& operator= ( const HierarchicalReport& );

        std::vector< SectionId > parent_;
        std::vector< std::vector<SectionId> > children_;
        std::vector< std::string > names_;
        std::vector< Data* > data_;
        SectionId curr_;
    };

    template<class Data>
    HierarchicalReport<Data>::HierarchicalReport ( )
    {
        // Create root node.
        parent_.push_back ( -1 );
        children_.push_back ( std::vector<SectionId>( ) );
        names_.push_back ( std::string ( "root" ) );
        data_.push_back ( 0 );
        curr_ = 0;
    }

    template<class Data>
    HierarchicalReport<Data>::~HierarchicalReport ( )
    {
        for ( typename std::vector< Data* >::iterator it = data_.begin ( ), end_it = data_.end ( );
              it != end_it; ++it )
        {
            delete *it;
        }
    }

    template<class Data>
    Data* HierarchicalReport<Data>::begin_section ( const std::string& name,
                                                    const Data& data )
    {
        // Create new section.
        SectionId child_id = parent_.size ( );
        parent_.push_back ( curr_ );
        children_.push_back ( std::vector<SectionId>( ) );
        names_.push_back ( name );
        data_.push_back ( new Data ( data ) );

        // Update children information.
        children_[curr_].push_back ( child_id );

        // Update curr_.
        curr_ = child_id;

        return data_.back ( );
    }

    template<class Data>
    Data* HierarchicalReport<Data>::end_section ( )
    {
        if ( curr_ == -1 )
        {
            // TODO: how should this error (mismatched begin/end) be handled?
            assert ( false );
            return 0;
        }

        Data* data = data_[curr_];
        curr_ = parent_[curr_];
        return data;
    }

    template<class Data>
    template<class Visitor>
    void HierarchicalReport<Data>::traverse_depth_first ( Visitor& visitor ) const
    {
        recursive_visit ( visitor, 0 );
    }

    template<class Data>
    template<class Visitor>
    void HierarchicalReport<Data>::recursive_visit ( Visitor& visitor, SectionId id ) const
    {
        visitor.enter ( names_[id], data_[id] );

        for ( std::vector<SectionId>::const_iterator it = children_[id].begin ( ),
              end_it = children_[id].end ( ); it != end_it; ++it )
        {
            recursive_visit ( visitor, *it );
        }

        visitor.exit ( names_[id], data_[id] );
    }

}

#endif
