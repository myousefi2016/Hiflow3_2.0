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

#ifndef HIFLOW_LOG_H
#    define HIFLOW_LOG_H

#    include "config.h"

#    include <iostream>
#    include <iterator>
#    include <map>
#    include <sstream>
#    include <string>
#    include <vector>

/// @brief This file contains facilities for logging output.
/// @author Staffan Ronnas, Thomas Gengenbach

/// @brief Turn log messages of category INFO on and off by setting
/// this variable.
extern bool INFO;

namespace hiflow
{

    /////////////////////////////////////////////////////////////////////////////
    //////////////// Log Machinery //////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    /// \brief Available categories for log messages.

    enum LogCategory
    {
        LOG_DEBUG_TAG, LOG_INFO_TAG, LOG_WARNING_TAG, LOG_ERROR_TAG
    };

    /// \brief Class representing a single entry in a log. For
    /// internal use by the class Log.

    class LogEntry
    {
        friend class Log;

        std::string file_name; //< source file where message originated
        std::string function_name; //< name of function where message originated
        int line_number; //< number of line where message originated
        LogCategory category; //< category of message
        std::ostringstream message; //< message data
    };

    /// \brief Class managing a collection of log entries.

    class Log
    {
      public:

        /// \brief Output the available entries to the output stream.
        void flush ( );

        /// \brief Delete all available entries.
        void clear ( );

        /// \return Pointer to current output stream.
        std::ostream* get_target ( ) const;

        /// \brief Set output stream
        /// \param[in] target   new output stream
        void set_target ( std::ostream* target );

        /// \brief Create a new entry in the log.
        /// \param[in] file          name of origin source file
        /// \param[in] function      name of origin function
        /// \param[in] line_number   number of origin line
        /// \param[in] category      category of log message
        /// \return output stream into which log message can be inserted.
        std::ostream& create_entry ( const char* file, const char* function,
                                     int line_number, LogCategory category );
      private:
        // All logs are managed by the singleton LogKeeper class.
        friend class LogKeeper;

        Log ( );
        explicit Log ( std::ostream* target );
        ~Log ( );

        // Copy and assignment is not allowed.
        Log ( const Log& );
        Log& operator= ( const Log& );

        std::vector<LogEntry*> entries_; //< Set of currently stored log entries
        std::ostream* target_; //< Target output stream
    };

    /// \brief Singleton class for managing Log objects.

    class LogKeeper
    {
      public:
        /// \brief Access a Log by its name.
        /// \param[in] log_name   name of log
        /// \return reference to the log
        static Log& get_log ( const char* log_name );

        /// \brief Create a new log.
        ///
        /// \param[in] log_name   name of new log
        /// \param[in] target     pointer to output target
        /// \return reference to newly created log
        static Log& create_log ( const char* log_name, std::ostream* target );
      private:
        static LogKeeper instance_; //< Singleton instance 

        LogKeeper ( );
        ~LogKeeper ( );

        static LogKeeper* instance ( );

        std::map<std::string, Log*> logs_; //< Map log_name -> log
    };

    /////////////////////////////////////////////////////////////////////////////
    //////////////// Convenience functions for output formatting ////////////////
    /////////////////////////////////////////////////////////////////////////////

    /// \brief Create a string representation of a range.
    ///
    /// \details It must be possible to output the elements of the
    /// range with operator<<.
    ///
    /// \param[in] begin   iterator to beginning of range
    /// \param[in] end     iterator to end of range
    /// \return std::string with space-delimited representation of elements in the range.

    template<typename Iterator>
    std::string string_from_range ( const Iterator& begin, const Iterator& end )
    {
        std::stringstream stream;
        std::copy ( begin, end, std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type > ( stream, " " ) );
        return stream.str ( );
    }

    /// \brief Create a high-precision string representation of a
    /// range of floating-point numbers.
    ///
    /// \details It must be possible to output the elements of the
    /// range with operator<<.
    ///
    /// \param[in] begin   iterator to beginning of range
    /// \param[in] end     iterator to end of range
    ///
    /// \return std::string with space-delimited representation of
    /// elements in the range, with 16-digit precision.

    template<typename Iterator>
    std::string precise_string_from_range ( const Iterator& begin, const Iterator& end )
    {
        std::stringstream stream;
        stream.precision ( 16 );
        std::copy ( begin, end, std::ostream_iterator<typename std::iterator_traits<Iterator>::value_type > ( stream, " " ) );
        return stream.str ( );
    }

    /// \brief Create a string representation of a pointer range.
    ///
    /// \details It must be possible to output the elements of the
    /// range with operator<<.
    ///
    /// \todo This function is a copy of string_from_range(), and it
    /// should be possible to replace its use with that function.
    ///
    /// \param[in] begin   pointer to beginning of range
    /// \param[in] end     pointer to end of range
    /// \return std::string with space-delimited representation of elements in the range.

    template<typename ValueType>
    std::string string_from_pointer_range ( const ValueType* begin, const ValueType* end )
    {
        std::stringstream stream;
        std::copy ( begin, end, std::ostream_iterator<ValueType>( stream, " " ) );
        return stream.str ( );
    }
} // namespace hiflow

/////////////////////////////////////////////////////////////////////////////
//////////////// Logging Macros /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/// These macros should be used to log messages in the different
/// categories. The use of macros makes it possible to capture
/// information about the position in the source code (file, function,
/// line number), and furthermore to completely remove e.g. log
/// messages for debugging when compiling in Release mode. Finally,
/// the "DEBUG_LEVEL" mechanism is implemented at this level.
///
/// Each macro takes its input as a stream, which makes it possible to
/// chain several values together easily. An example of how it can be used:
///
/// LOG_DEBUG(1, "Hello, " << your_name << ", the number of the day is " << number);
///
/// Here, "your_name" and "number" are variables in the current code context.

/// Use of the LOG_DEBUG macros is only activated in Debug compilation
/// mode. In Release mode, the macro expands to an empty line.

#    ifndef NDEBUG  

/// Log a message with a given level "lvl". The message is only logged if lvl is not greater than the
/// integer variable DEBUG_LEVEL, which must be defined whereever the
/// macro is used. Debug messages are logged to the "debug" log.
#        define LOG_DEBUG(lvl, stream) { if (lvl <= DEBUG_LEVEL) {            \
            hiflow::LogKeeper::get_log("debug").create_entry(          \
            __FILE__, __FUNCTION__, __LINE__, hiflow::LOG_DEBUG_TAG)   \
            << "DBG: " << stream ; hiflow::LogKeeper::get_log("debug").flush(); }}

#    else 

#        define LOG_DEBUG(lvl, stream) 

#    endif   // NDEBUG

/// Log an error message to the log "error".
#    define LOG_ERROR(stream) { hiflow::LogKeeper::get_log("error").create_entry(   \
            __FILE__, __FUNCTION__, __LINE__, hiflow::LOG_ERROR_TAG)             \
            << "ERR: " << stream ; hiflow::LogKeeper::get_log("error").flush(); }

/// Log a warning message to the log "warning".
#    define LOG_WARNING(stream) { hiflow::LogKeeper::get_log("warning").create_entry( \
            __FILE__, __FUNCTION__, __LINE__, hiflow::LOG_WARNING_TAG)             \
            << "WRN: " << stream ; hiflow::LogKeeper::get_log("warning").flush(); }

/// Log an informational message to the log "info". Info messages are
/// tagged with an additional arbitrary "name", which makes it
/// possible to filter them more easily.
#    ifdef INFO_LOG_ACTIVE
#        define LOG_INFO(name, stream) { if (INFO) { \
            hiflow::LogKeeper::get_log("info").create_entry(            \
            __FILE__, __FUNCTION__, __LINE__, hiflow::LOG_INFO_TAG)     \
            << name << " : " << stream ; hiflow::LogKeeper::get_log("info").flush(); }}
#    else
#        define LOG_INFO(name, stream) {}
#    endif
#endif // HIFLOW_LOG_H
