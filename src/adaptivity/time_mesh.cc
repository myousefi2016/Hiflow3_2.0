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

#include "time_mesh.h"
#include "common/log.h"

#include <cmath>
#include <utility>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/function.hpp>

const int DEBUG_LEVEL = 1;

namespace hiflow
{

    template<class DataType>
    TimeMesh<DataType>::TimeMesh ( )
    {
        SortedArray<DataType> tmp;
        tmp.insert ( 0. );
        this->times_.push_back ( tmp );
    }

    template<class DataType>
    void TimeMesh<DataType>::clear ( )
    {
        this->times_.clear ( );
        SortedArray<DataType> tmp;
        tmp.insert ( 0. );
        this->times_.push_back ( tmp );
    }

    template<class DataType>
    void TimeMesh<DataType>::add_times ( std::vector<DataType>& add_times, int counter )
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        for ( int t = 0; t < add_times.size ( ); ++t )
        {
            this->times_[counter].insert ( add_times[t] );
        }
    }

    template<class DataType>
    void TimeMesh<DataType>::add_time ( DataType add_time, int counter )
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        LOG_DEBUG ( 2, "add " << add_time << " to curent times: " << string_from_range ( this->times_[counter].data ( ).begin ( ), this->times_[counter].data ( ).end ( ) ) );
        this->times_[counter].insert ( add_time );
    }

    template<class DataType>
    void TimeMesh<DataType>::set_time ( int time_step, DataType add_time, int counter )
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        this->times_[counter].data ( ).at ( time_step ) = add_time;
    }

    template<class DataType>
    DataType TimeMesh<DataType>::delta_t ( int time_step, int counter ) const
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        if ( time_step < 1 )
        {
            return 0.;
        }

        if ( time_step >= this->times_[counter].size ( ) )
        {
            return 0.;
        }
        return (this->time ( time_step, counter ) - this->time ( time_step - 1, counter ) );
    }

    template<class DataType>
    std::vector<DataType> TimeMesh<DataType>::get_all_times ( int counter ) const
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        return this->times_[counter].data ( );
    }

    template<class DataType>
    int TimeMesh<DataType>::num_intervals ( int counter ) const
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        return this->times_[counter].size ( ) - 1;
    }

    template<class DataType>
    int TimeMesh<DataType>::num_levels ( ) const
    {
        return this->times_.size ( );
    }

    template<class DataType>
    int TimeMesh<DataType>::num_intervals ( ) const
    {
        int counter = this->times_.size ( ) - 1;
        return this->num_intervals ( counter );
    }

    template<class DataType>
    DataType TimeMesh<DataType>::time ( int time_step, int counter ) const
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        if ( time_step < this->times_[counter].size ( ) && time_step >= 0 )
        {
            return this->times_[counter].data ( ).at ( time_step );
        }
        return -1.;
    }

    template<class DataType>
    DataType TimeMesh<DataType>::time ( int time_step ) const
    {
        int counter = this->times_.size ( ) - 1;
        return this->time ( time_step, counter );
    }

    template<class DataType>
    DataType TimeMesh<DataType>::start ( int counter ) const
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        return this->time ( 0, counter );
    }

    template<class DataType>
    DataType TimeMesh<DataType>::start ( ) const
    {
        int counter = this->times_.size ( ) - 1;
        return this->start ( counter );
    }

    template<class DataType>
    DataType TimeMesh<DataType>::end ( int counter ) const
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        return this->time ( times_[counter].size ( ) - 1, counter );
    }

    template<class DataType>
    DataType TimeMesh<DataType>::end ( ) const
    {
        int counter = this->times_.size ( ) - 1;

        return this->end ( counter );
    }

    template<class DataType>
    void TimeMesh<DataType>::refine ( const std::vector<int>& markers, int counter )
    {
        int num_intervals = markers.size ( );
        assert ( num_intervals == this->num_intervals ( counter ) );

        int old_end = this->times_.size ( );

        this->times_.resize ( counter + 2 );
        for ( int l = old_end; l<this->times_.size ( ); ++l )
        {
            SortedArray<DataType> tmp_time;
            tmp_time.insert ( 0. );
            this->times_[l] = tmp_time;
        }

        int index = counter + 1;
        assert ( index < this->times_.size ( ) );

        // loop over old time intervals
        int interval = 0;
        while ( interval < num_intervals )
        {
            int c_l = markers[interval];
            int c_r = markers[interval + 1];

            DataType t_p = this->time ( interval, counter );
            DataType t_c = this->time ( interval + 1, counter );
            DataType t_n = this->time ( interval + 2, counter );

            // t_p __________ t_c ____________ t_n
            //      interval       interval+1

            // check if coarsening is possible
            bool coarsen = true;
            if ( c_l > 0 || c_r > 0 ) coarsen = false; // no coarsening if one of the intervals is marked for refinement
            if ( c_l + c_r > -2 ) coarsen = false; // no coarsening if number of assigned coarsening markers is not high enough

            if ( coarsen )
            {
                this->add_time ( t_n, index );
            }
            else
            {
                if ( c_l <= 0 ) // keep left interval
                {
                    this->add_time ( t_c, index );
                }
                if ( c_l > 0 ) // refine left interval
                {
                    this->add_time ( 0.5 * ( t_p + t_c ), index );
                    this->add_time ( t_c, index );
                }
                if ( c_r <= 0 ) // keep right interval
                {
                    this->add_time ( t_n, index );
                }
                if ( c_r > 0 ) // refine righ tinterval
                {
                    this->add_time ( 0.5 * ( t_c + t_n ), index );
                    this->add_time ( t_n, index );
                }
            }
            interval += 2;
        }
        this->make_regular ( index );
    }

    template<class DataType>
    void TimeMesh<DataType>::make_regular ( int counter )
    {
        assert ( this->times_.size ( ) > counter );
        assert ( counter >= 0 );

        int num_intervals = this->num_intervals ( counter );
        if ( num_intervals % 2 == 0 )
        {
            return;
        }

        DataType t_prev_end = this->time ( num_intervals - 1, counter );
        DataType t_end = this->time ( num_intervals, counter );

        LOG_DEBUG ( 2, "Time mesh before regularization " << string_from_range ( this->times_[counter].begin ( ), this->times_[counter].end ( ) ) );

        this->set_time ( num_intervals, 0.5 * ( t_prev_end + t_end ), counter );
        this->add_time ( t_end, counter );

        LOG_DEBUG ( 2, "Time mesh after regularization " << string_from_range ( this->times_[counter].begin ( ), this->times_[counter].end ( ) ) );
    }

    template<class DataType>
    void TimeMesh<DataType>::load_all ( const std::string& prefix, int num_mesh )
    {
        assert ( num_mesh > 0 );

        this->times_.clear ( );
        this->times_.resize ( num_mesh );
        for ( int l = 0; l < num_mesh; ++l )
        {
            this->load ( prefix, l );
            LOG_DEBUG ( 1, "Loaded time mesh on level " << l << ": " << string_from_range ( this->times_[l].begin ( ), this->times_[l].end ( ) ) );
        }
    }

    template<class DataType>
    void TimeMesh<DataType>::load ( const std::string& prefix, int counter )
    {

        std::stringstream pre;
        pre << prefix << "." << counter << ".csv";
        std::string filename = pre.str ( );

        this->times_.resize ( counter + 1 );
        this->times_[counter].clear ( );

        std::ifstream in;
        in.open ( filename.c_str ( ) );
        std::vector < std::vector<DataType> > data;

        const std::string& delim = "";
        std::string line;
        std::string strnum;

        while ( getline ( in, line ) )
        {
            data.push_back ( std::vector<DataType>( ) );

            for ( std::string::const_iterator i = line.begin ( ); i != line.end ( ); ++i )
            {
                // If i is not a delim, then append it to strnum
                if ( delim.find ( *i ) == std::string::npos )
                {
                    strnum += *i;
                    if ( i + 1 != line.end ( ) ) // If it's the last char, do not continue
                        continue;
                }
            }

            // if strnum is still empty, it means the previous char is also a
            // delim (several delims appear together). Ignore this char.
            if ( strnum.empty ( ) )
                continue;

            // If we reach here, we got a number. Convert it to DataType.
            DataType number;

            std::istringstream ( strnum ) >> number;
            data.back ( ).push_back ( number );

            strnum.clear ( );
        }
        in.close ( );

        //LOG_DEBUG(1, "Read in time mesh: num rows: " << data.size());
        //LOG_DEBUG(1, "Read in time mesh for level " << counter << ": " << string_from_range(data[0].begin(), data[0].end()));

        for ( int l = 0; l < data.size ( ); ++l )
        {
            if ( data[l].size ( ) > 0 )
            {
                //std::cout << data[l][0] << std::endl;
                this->add_time ( data[l][0], counter );
            }
        }
    }

    template<class DataType>
    void TimeMesh<DataType>::save_all ( const std::string& prefix ) const
    {
        int num_mesh = this->times_.size ( );
        for ( int l = 0; l < num_mesh; ++l )
        {
            this->save ( prefix, l );
        }
    }

    template<class DataType>
    void TimeMesh<DataType>::save ( const std::string& prefix, int counter ) const
    {
        assert ( counter >= 0 );
        assert ( counter < this->times_.size ( ) );

        std::stringstream pre;
        pre << prefix << "." << counter << ".csv";
        std::string filename = pre.str ( );

        std::ofstream myfile;
        myfile.open ( filename.c_str ( ) );

        for ( int t = 0; t<this->num_intervals ( counter ) + 1; ++t )
        {
            myfile << this->time ( t, counter ) << "\n";
        }
        myfile << "\n";
        myfile.close ( );
    }

    template class TimeMesh <double>;
}
