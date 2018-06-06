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

/// @author Philipp Gerstner

#ifndef HIFLOW_EIGEN_VALUE_SLEPC_ENVIRONMENT_H_
#    define HIFLOW_EIGEN_VALUE_SLEPC_ENVIRONMENT_H_

namespace hiflow
{
    namespace la
    {

        /// Control SLEPc environment initialization
        /// Before any SLEPc function will be called the initialize function must called.
        /// The initialize and finalize functions are idempotent.

        class SLEPcEnvironment
        {
          public:
            /// Destructor
            ~SLEPcEnvironment ( );

            /// Prepare SLEPc environment with arguments
            static void initialize ( int argc, char **argv );

            /// Prepare SLEPc environment without arguments
            static void initialize ( );

            /// Destroy SLEPc environment
            static void finalize ( );

          private:
            /// No construction
            SLEPcEnvironment ( );

            /// No copy
            SLEPcEnvironment ( const SLEPcEnvironment& );
            SLEPcEnvironment& operator= ( const SLEPcEnvironment& );

            static bool initialized_;
        };

    } // namespace la
} // namespace hiflow

#endif  
