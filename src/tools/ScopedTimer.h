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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-06-15

#ifndef SRC_TOOLS_TIMER_H_
#    define SRC_TOOLS_TIMER_H_

#    include <string>
#    include <sys/time.h>

namespace hiflow
{

    /**
     * Measure the time between construction and destruction of an ScopedTimer object.
     * Time is printed automatically at destruction if print was not called before.
     */
    class ScopedTimer
    {
      public:
        /// Constructor starts the timer

        ScopedTimer ( std::string const& label, bool startImmediately = true )
        : startTime ( 0.0 ), duration ( 0.0 ), label ( label ), nbRepetitions ( 0 ), running ( false ), printed ( false )
        {
            if ( startImmediately ) start ( );
        }

        /// Destructor stops the timer and print the time

        ~ScopedTimer ( )
        {
            if ( !printed ) print ( );
        }

        /// Manually start the timer. Reseting the timer if already running.

        void start ( )
        {
            gettimeofday ( &now, NULL );
            startTime = now.tv_sec * 1000000.0 + now.tv_usec;
            running = true;
            if ( !nbRepetitions ) repetition ( );
        }

        /// Signal that a new repetition is started

        void repetition ( )
        {
            ++nbRepetitions;
        }

        /// Manually stop the timer. Do nothing if already not running.

        void stop ( )
        {
            if ( !running ) return;
            gettimeofday ( &now, NULL );
            double endTime = now.tv_sec * 1000000.0 + now.tv_usec;
            duration += ( endTime - startTime ) / 1000000;
            running = false;
        }

        /// Manually print the time.

        void print ( )
        {
            stop ( );
            std::cout << "Time for " << label << ": " << duration / nbRepetitions << " s";
            if ( nbRepetitions > 1 ) std::cout << " (number of repetitions = " << nbRepetitions << ")";
            std::cout << std::endl;
            printed = true;
        }

      private:

        /// Current time
        struct timeval now;

        /// Start time in microseconds
        double startTime;

        /// Time in seconds
        double duration;

        /// Label for printing
        std::string label;

        /// Number of repetitions for more accurate benchmarks
        /// Duration will be divided by nbRepetitions
        int nbRepetitions;

        /// Is timer running?
        bool running;

        /// Was time already printed?
        bool printed;
    };

    template <bool Timer>
    void start ( ScopedTimer *timer );

    template <>
    inline void start<false>( ScopedTimer *timer )
    {
    }

    template <>
    inline void start<true>( ScopedTimer *timer )
    {
#    ifdef WITH_CUDA
        cuda_sync_threads ( );
#    endif
        timer->start ( );
    }

    template <bool Timer>
    void stop ( ScopedTimer *timer );

    template <>
    inline void stop<false>( ScopedTimer *timer )
    {
    }

    template <>
    inline void stop<true>( ScopedTimer *timer )
    {
#    ifdef WITH_CUDA
        cuda_sync_threads ( );
#    endif
        timer->stop ( );
    }

    template <bool Timer>
    void reset ( ScopedTimer *timer );

    template <>
    inline void reset<false>( ScopedTimer *timer )
    {
    }

    template <>
    inline void reset<true>( ScopedTimer *timer )
    {
        timer->repetition ( );
    }

} // namespace hiflow

#endif /* SRC_TOOLS_TIMER_H_ */
