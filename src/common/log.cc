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

#include "common/log.h"

#include <cassert>
#include <cstring>
#include <ctime>
#include <stdexcept>

// INFO log is actived by default
bool INFO = true;

namespace hiflow
{

    //////////////// Log ////////////////

    Log::Log ( )
    : target_ ( &std::cerr )
    {
    }

    Log::Log ( std::ostream* target )
    : target_ ( target )
    {
    }

    Log::~Log ( )
    {
        flush ( );
    }

    void Log::flush ( )
    {
        if ( target_ )
        {
            for ( std::vector<LogEntry*>::const_iterator it = entries_.begin ( );
                  it != entries_.end ( ); ++it )
            {
                if ( *it )
                {
                    *target_ << "[" << ( *it )->file_name << ":" << ( *it )->line_number << "] "
                            << ( *it )->message.str ( ) << "\n";
                }
            }
        }
        clear ( );
    }

    void Log::clear ( )
    {
        for ( std::vector<LogEntry*>::iterator it = entries_.begin ( );
              it != entries_.end ( ); ++it )
        {
            delete *it;
            *it = 0;
        }
        entries_.clear ( );
    }

    std::ostream* Log::get_target ( ) const
    {
        return target_;
    }

    void Log::set_target ( std::ostream* target )
    {
        target_ = target;
    }

    std::ostream& Log::create_entry ( const char* file, const char* function,
                                      int line_number, LogCategory category )
    {
        entries_.push_back ( new LogEntry ( ) );
        LogEntry* entry = entries_.back ( );

        // strip path from filename (only on UNIX-like systems)
        const char* last_slash = strrchr ( file, '/' );
        last_slash = last_slash ? last_slash + 1 : file;
        entry->file_name = std::string ( last_slash );
        entry->function_name = std::string ( function );
        entry->line_number = line_number;
        entry->category = category;
        return entry->message;
    }

    //////////////// LogKeeper ////////////////
    LogKeeper LogKeeper::instance_;

    Log& LogKeeper::get_log ( const char* log_name )
    {
        const std::string name ( log_name );
        if ( instance ( )->logs_.count ( name ) == 0 )
        {
            assert ( 0 );
            std::cerr << "Log " << log_name << " not found!\n";
            throw std::runtime_error ( "Log not found!" );
        }

        return *( instance ( )->logs_[name] );
    }

    Log& LogKeeper::create_log ( const char* log_name, std::ostream* target )
    {
        const std::string name ( log_name );
        if ( instance ( )->logs_.count ( name ) == 0 )
        {
            instance ( )->logs_[name] = new Log ( target );
#ifdef INFO_LOG_ACTIVE
            if ( name == "info" )
            {
                time_t tp;
                time ( &tp );
                instance ( )->logs_[name]->create_entry ( "INFO", "", 0, LOG_INFO_TAG ) << "Time: " << std::asctime ( localtime ( &tp ) );
            }
#endif
        }
        return *( instance ( )->logs_[name] );
    }

    LogKeeper::LogKeeper ( )
    {
        create_log ( "debug", 0 );
        create_log ( "info", 0 );
        create_log ( "warning", 0 );
        create_log ( "error", 0 );
    }

    LogKeeper::~LogKeeper ( )
    {
        for ( std::map<std::string, Log*>::iterator it = logs_.begin ( ); it != logs_.end ( ); ++it )
        {
            Log* log = it->second;
            log->flush ( );
            delete log;
            log = 0;
        }
    }

    LogKeeper* LogKeeper::instance ( )
    {
        return &instance_;
    }

} // namespace hiflow
