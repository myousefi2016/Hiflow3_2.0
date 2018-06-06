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

/// @author Staffan Ronnas, Thomas Gengenbach

#ifndef _TIMER_H_
#    define _TIMER_H_

#    include <sys/time.h>
#    include <iostream>

namespace hiflow
{

    /// Stopwatch class useful for measuring periods of time in code execution.
    /// Based on gettimeofday.
    /// In general, this class is not thread-safe.

    class Timer
    {
      public:
        /// Construct timer object and starts the timer.
        inline Timer ( );

        /// Restart the timer.
        inline void start ( );

        /// Stop the timer and add the time since the last call to start() (or the constructor)
        /// to its duration.
        inline void stop ( );

        /// Reset the duration to zero.
        inline void reset ( );

        /// Get duration in seconds as a double value.
        inline double get_duration ( ) const;

        /// Get integer part of duration (in seconds)
        inline long get_seconds ( ) const;

        /// Get fractional part of duration (in nanoseconds)
        inline long get_micro_seconds ( ) const;

      private:
        static const long us_per_s = 1000000L;

        inline void update_duration ( timeval start, timeval stop );
        inline timeval add_ts ( timeval ts1, timeval ts2 );
        inline timeval subtract_ts ( timeval begin, timeval end );

        timeval _start;
        timeval _duration;
        bool _is_running;

    };

    /// Output operator
    std::ostream& operator<< ( std::ostream& os, const Timer& timer );

    /// Inline implementation

    Timer::Timer ( )
    {
        reset ( );
        start ( );
    }

    void Timer::start ( )
    {
        _is_running = true;

        // take time as late as possible
        gettimeofday ( &_start, 0 );
    }

    void Timer::stop ( )
    {
        // take the time as soon as possible
        timeval stop;
        gettimeofday ( &stop, 0 );

        if ( _is_running )
        {
            _is_running = false;
            update_duration ( _start, stop );
        }
    }

    void Timer::reset ( )
    {
        _is_running = false;
        _duration.tv_sec = 0;
        _duration.tv_usec = 0;
    }

    double Timer::get_duration ( ) const
    {

        return get_seconds ( ) +
                static_cast < double > ( get_micro_seconds ( ) ) / static_cast < double > ( us_per_s );
    }

    long Timer::get_seconds ( ) const
    {
        return _duration.tv_sec;
    }

    long Timer::get_micro_seconds ( ) const
    {
        return _duration.tv_usec;
    }

    void Timer::update_duration ( timeval start, timeval stop )
    {
        timeval diff_ts = subtract_ts ( start, stop );
        _duration = add_ts ( _duration, diff_ts );
    }

    timeval Timer::add_ts ( timeval ts1, timeval ts2 )
    {
        timeval res;
        res.tv_sec = ts1.tv_sec + ts2.tv_sec;
        res.tv_usec = ts1.tv_usec + ts2.tv_usec;

        // carry over extra seconds
        while ( res.tv_usec > us_per_s )
        {
            ++res.tv_sec;
            res.tv_usec -= us_per_s;
        }

        return res;

    }

    timeval Timer::subtract_ts ( timeval begin, timeval end )
    {
        timeval res;

        res.tv_sec = end.tv_sec - begin.tv_sec;
        res.tv_usec = end.tv_usec - begin.tv_usec;

        // correct for negative micro-second difference     
        if ( res.tv_usec < 0 )
        {
            --res.tv_sec;
            res.tv_usec += us_per_s;
        }

        return res;
    }
}

#endif /* _TIMER_H_ */
