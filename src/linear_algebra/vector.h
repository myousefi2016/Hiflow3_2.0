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

/// @author Simon Gawlok, Martin Wlotzka

#ifndef HIFLOW_LINEAR_ALGEBRA_VECTOR_H
#    define HIFLOW_LINEAR_ALGEBRA_VECTOR_H

namespace hiflow
{
    namespace la
    {

        /// \brief Abstract base class for distributed vector implementations.

        template<class DataType>
        class Vector
        {
          public:
            /// Standard constructor

            Vector ( )
            {
            }
            /// Destructor

            virtual ~Vector ( )
            {
            }

            /// Local size of Vector
            virtual int size_local ( ) const = 0;
            /// Global size of Vector
            virtual int size_global ( ) const = 0;

            /// Clone this Vector
            virtual Vector<DataType>* Clone ( ) const = 0;
            /// Update operator, i.e. exchange values for distributed vectors
            virtual void Update ( ) = 0;
            /// Initiate update
            virtual void begin_update ( ) = 0;
            /// Finalize update
            virtual void end_update ( ) = 0;

            /// Set Vector to zero
            virtual void Zeros ( ) = 0;
            /// Get value at a known index
            virtual DataType GetValue ( const int index ) const = 0;
            /// Get values at given indices
            virtual void GetValues ( const int* indices, const int length, DataType* values ) const = 0;
            /// Euclidean length of vector
            virtual DataType Norm2 ( ) const = 0;
            /// l1 norm of vector
            virtual DataType Norm1 ( ) const = 0;
            /// maximum absolute value of entries
            virtual DataType NormMax ( ) const = 0;

            /// Scalar product
            virtual DataType Dot ( const Vector<DataType>& vec ) const = 0;

            /// Add value to a given index
            virtual void Add ( int index, DataType scalar ) = 0;
            /// Add values to given indices
            virtual void Add ( const int* indices, int length, const DataType* values ) = 0;
            /// Set given global index to value
            virtual void SetValue ( const int index, const DataType value ) = 0;
            /// Set given global indices to given values
            virtual void SetValues ( const int* indices, const int length, const DataType* values ) = 0;

            /// this <- this + alpha * vecx
            virtual void Axpy ( const Vector<DataType>& vecx, const DataType alpha ) = 0;

            /// this <- alpha * this + vecx
            virtual void ScaleAdd ( const Vector<DataType>& vecx, const DataType alpha ) = 0;

            /// this <- alpha * this
            virtual void Scale ( const DataType alpha ) = 0;
        };

    }
}

#endif
