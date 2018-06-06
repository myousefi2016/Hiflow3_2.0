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

#include "periodic_boundary_conditions.h"
#include "common/macros.h"

/// @author Martin Baumann

namespace hiflow
{

    template<class DataType>
    void PeriodicBoundaryConditions<DataType>::handle_multiple_correspondences ( )
    {
        bool debug_print = false;

        std::map<doffem::DofID, std::vector<doffem::DofID> >::iterator it;
        it = corresponding_dof_.begin ( );
        while ( it != corresponding_dof_.end ( ) )
        {
            std::map<doffem::DofID, std::vector<doffem::DofID> >::iterator it2;
            it2 = corresponding_dof_.begin ( );
            while ( it2 != corresponding_dof_.end ( ) )
            {
                if ( it2 != it )
                {
                    bool remove_it2 = false;

                    // key value of it2 matches key value of it
                    if ( it2->first == it->first )
                    {
                        for ( int i = 0; i < it2->second.size ( ); ++i )
                            it->second.push_back ( it2->second[i] );
                        remove_it2 = true;
                        if ( debug_print )
                            std::cout << "MODE1: " << it->first
                                << "\t" << it2->first << std::endl;
                    }

                    // one of the values in the vector of it2 matches key value of it
                    bool found_matching = false;
                    for ( int i = 0; i < it2->second.size ( ); ++i )
                    {
                        if ( it2->second[i] == it->first )
                        {
                            found_matching = true;
                            break;
                        }
                        if ( debug_print )
                            std::cout << "MODE2: " << it->first
                                << "\t" << it2->second.at ( 0 )
                            << ", SIZE = " << it2->second.size ( ) << std::endl;
                    }
                    if ( found_matching == true )
                    {
                        it->second.push_back ( it2->first );
                        for ( int i = 0; i < it2->second.size ( ); ++i )
                            if ( it2->second[i] != it->first )
                                it->second.push_back ( it2->second[i] );
                        remove_it2 = true;
                    }

                    // key value of it2 matches one of the values in it
                    for ( int i = 0; i < it->second.size ( ); ++i )
                    {
                        if ( it2->first == it->second[i] )
                        {
                            for ( int j = 0; j < it2->second.size ( ); ++j )
                                it->second.push_back ( it2->second[j] );
                            remove_it2 = true;
                            if ( debug_print )
                                std::cout << "MODE3" << std::endl;
                            break;
                        }
                    }

                    // value in it2 matches value in it
                    for ( int i = 0; i < it->second.size ( ); ++i )
                    {
                        bool found_matching = false;
                        for ( int j = 0; j < it2->second.size ( ); ++j )
                        {
                            if ( it->second[i] == it2->second[j] )
                                found_matching = true;
                        }
                        if ( found_matching == true )
                        {
                            it->second.push_back ( it2->first );
                            for ( int j = 0; j < it2->second.size ( ); ++j )
                                if ( it2->second[j] != it->second[i] )
                                    it->second.push_back ( it2->second[j] );
                            remove_it2 = true;
                            if ( debug_print )
                                std::cout << "MODE4" << std::endl;
                            break;
                        }
                    }

                    // remove it2

                    if ( remove_it2 == true )
                        corresponding_dof_.erase ( it2 );
                    else
                        if ( debug_print )
                        std::cout << "--------- NOTHING ------- " << std::endl;
                }
                ++it2;
            }
            ++it;
        }

        // cleaning

        it = corresponding_dof_.begin ( );
        while ( it != corresponding_dof_.end ( ) )
        {
            // erase entry if DoF in vector is same as DoF in key
            std::vector<doffem::DofID>::iterator values;
            for ( values = it->second.begin ( ); values != it->second.end ( ); ++values )
                if ( *values == it->first )
                    values = it->second.erase ( values );

            // erase DoF in vector if it occurs multiply
            for ( values = it->second.begin ( ); values != it->second.end ( ); ++values )
            {
                // first sort vector than erase double entries
                sort ( it->second.begin ( ), it->second.end ( ) );
                std::vector<doffem::DofID>::iterator last_unique;
                last_unique = std::unique ( it->second.begin ( ), it->second.end ( ) );
                it->second.resize ( last_unique - it->second.begin ( ) );
            }

            ++it;
        }

        // print correspondence list

        if ( debug_print )
        {
            std::cout << "PeriodicBoundary Correspondence List" << std::endl;
            std::map<doffem::DofID, std::vector<doffem::DofID> >::const_iterator it;
            for ( it = corresponding_dof_.begin ( ); it != corresponding_dof_.end ( ); ++it )
            {
                std::cout << "  " << it->first << "\t <-> \t";
                std::vector<doffem::DofID>::const_iterator itc;
                for ( itc = it->second.begin ( ); itc != it->second.end ( ); ++itc )
                    std::cout << *itc << "\t";
                std::cout << std::endl;
            }
        }

    }

    template<class DataType>
    bool PeriodicBoundaryConditions<DataType>::insert_dof_id_entry ( MaterialNumber mat,
                                                                     int var,
                                                                     doffem::DofID dof_id )
    {
        std::map< std::pair<MaterialNumber, int>, std::vector<doffem::DofID> >::iterator it;
        it = dof_id_map_.find ( std::make_pair ( mat, var ) );

        // insert new map entry for tuple (mat,var) if none existing
        if ( it == dof_id_map_.end ( ) )
        {
            std::vector<doffem::DofID> ids;
            dof_id_map_.insert ( std::make_pair ( std::make_pair ( mat, var ), ids ) );
            it = dof_id_map_.find ( std::make_pair ( mat, var ) );
        }

        // is dof_id already in vector of DoFs?
        std::vector<doffem::DofID>& dof_vector = it->second;
        if ( find ( dof_vector.begin ( ), dof_vector.end ( ), dof_id ) == dof_vector.end ( ) )
        {
            dof_id_map_.find ( std::make_pair ( mat, var ) )->second.push_back ( dof_id );
            return true;
        }

        //std::cout << "ALREADY INSERTED: " << mat << "\t" << var << "\t" << dof_id << std::endl;

        return false;
    }

    template<class DataType>
    void PeriodicBoundaryConditions<DataType>::apply_boundary_conditions ( VectorSpace<DataType>& space )
    {
        // calculates corresponding DoFs
        compute_conditions ( space );

        // adapt DoFs
        change_dofs ( space );
    }

    template<class DataType>
    void PeriodicBoundaryConditions<DataType>::add_boundary_tuple ( MaterialNumber a, MaterialNumber b )
    {
        boundary_descriptor_.insert ( std::make_pair ( a, b ) );
        boundary_descriptor_.insert ( std::make_pair ( b, a ) );
    }

    /// there are boundary conditions for a boundary with material number 'mat' if there is an entry in boundary_descriptor_

    template<class DataType>
    bool PeriodicBoundaryConditions<DataType>::is_periodic_boundary ( MaterialNumber mat ) const
    {
        std::map<MaterialNumber, MaterialNumber>::const_iterator it;
        it = boundary_descriptor_.find ( mat );
        if ( it != boundary_descriptor_.end ( ) )
            return true;

        return false;
    }

    template<class DataType>
    void PeriodicBoundaryConditions<DataType>::change_dofs ( VectorSpace<DataType>& space )
    {
        // Calculate permutation: Blind DoFs should be replaced by real DoFs

        std::vector<int> perm;
        perm.resize ( space.dof ( ).get_nb_dofs ( ) );
        for ( int i = 0; i < perm.size ( ); ++i )
            perm[i] = i;

        std::map<doffem::DofID, std::vector<doffem::DofID> > & cd = corresponding_dof_;
        std::map<doffem::DofID, std::vector<doffem::DofID> >::const_iterator cd_it = cd.begin ( );
        while ( cd_it != cd.end ( ) )
        {
            // cd.first  -> master dof
            // cd.second -> vector of slave dofs
            for ( int slave = 0; slave < cd_it->second.size ( ); ++slave )
            {
                assert ( cd_it->second.at ( slave ) < perm.size ( ) );
                perm[cd_it->second[slave]] = cd_it->first;
            }
            ++cd_it;
        }

        // Make permutation consecutive

        std::vector<bool> is_used ( space.dof ( ).get_nb_dofs ( ), true );
        cd_it = cd.begin ( );
        while ( cd_it != cd.end ( ) )
        {
            for ( int slave = 0; slave < cd_it->second.size ( ); ++slave )
                is_used[cd_it->second[slave]] = false;
            ++cd_it;
        }
        std::vector<int> consecutive_perm;
        consecutive_perm.resize ( is_used.size ( ), -1 );
        int counter = -1;
        for ( int i = 0; i < is_used.size ( ); ++i )
        {
            if ( is_used[i] == true )
                ++counter;
            consecutive_perm[i] = counter;
        }

        // Mix permutations together

        for ( int i = 0; i < perm.size ( ); ++i )
            perm[i] = consecutive_perm[perm[i]];

        // Apply permutation

        space.dof ( ).apply_permutation ( perm, "PeriodicBoundaryConditions" );
    }

    template<class DataType>
    std::ostream& operator<< ( std::ostream& s, const PeriodicBoundaryConditions<DataType>& cond )
    {
        s << "PeriodicBoundaryConditions" << std::endl;
        s << "  Corresponding Boundaries:" << std::endl;
        typename std::map<typename PeriodicBoundaryConditions<DataType>::MaterialNumber, typename PeriodicBoundaryConditions<DataType>::MaterialNumber>::const_iterator it2;
        it2 = cond.boundary_descriptor_.begin ( );
        while ( it2 != cond.boundary_descriptor_.end ( ) )
        {
            s << "    " << it2->first << " <-> " << it2->second << std::endl;
            ++it2;
        }
        s << "  Boundary DoFs on (mat, var):" << std::endl;
        typename std::map < std::pair<typename PeriodicBoundaryConditions<DataType>::MaterialNumber, int>, std::vector<doffem::DofID> >::const_iterator it = cond.dof_id_map_.begin ( );
        while ( it != cond.dof_id_map_.end ( ) )
        {
            s << "    (" << it->first.first << ", " << it->first.second << "): " << it->second.size ( ) << std::endl;
            ++it;
        }
        s << "  # corresponding DoFs: " << cond.corresponding_dof_.size ( ) << std::endl;

        return s;
    }

    template class PeriodicBoundaryConditions<double>;
    template class PeriodicBoundaryConditions<float>;

} // namespace hiflow
