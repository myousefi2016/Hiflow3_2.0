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
/// @date 2015-12-04

#ifndef HIFLOW_LINEARALGEBRA_PETSC_ENVIRONMENT_H_
#    define HIFLOW_LINEARALGEBRA_PETSC_ENVIRONMENT_H_

namespace hiflow
{
    namespace la
    {

        /// Control PETSc environment initialization
        /// Before any PETSc function will be called the initialize function must called.
        /// This is typically done within the PETSc wrapper classes for Vec, Mat, or KSP.
        /// The initialize and finalize functions are idempotent.

        class PETScEnvironment
        {
          public:
            /// Destructor
            ~PETScEnvironment ( );

            /// Prepare PETSc environment with arguments
            static void initialize ( int argc, char **argv );

            /// Prepare PETSc environment without arguments
            static void initialize ( );

            /// Destroy PETSc environment
            static void finalize ( );

          private:
            /// No construction
            PETScEnvironment ( );

            /// No copy
            PETScEnvironment ( const PETScEnvironment& );
            PETScEnvironment& operator= ( const PETScEnvironment& );

            static bool initialized_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_ENVIRONMENT_H_
