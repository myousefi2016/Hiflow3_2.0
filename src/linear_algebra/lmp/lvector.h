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

/// @author Dimitar Lukarski

#ifndef __LVECTOR_H
#    define __LVECTOR_H

#    include <iostream>
#    include <cstring>
#    include <cmath>
#    include "la_global.h"
#    include "lmp_log.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Provides the base class to all local multi-platform vectors
        /// @author Dimitar Lukarski
        ///
        /// lVector (local Vector) - lmpLAtoolbox (Local Multi-Platform Linear Algebra Tool Box)
        /// Provides the base class to all local vectors for different platforms
        /// (like CPU,GPU,etc) with different implementations (sequential,parallel,mkl,cblas,...)

        template <typename ValueType>
        class lVector
        {
          public:

            lVector ( );
            /// The virtual destructor free the buffers of the vector
            virtual ~lVector ( );

            /// Return the size of the vector
            /// @return The size of the vector
            inline int get_size ( void ) const;

            /// Return the name of the vector
            /// @return The name of the vector
            std::string get_name ( void ) const;

            /// Provides human-readable information for the vector - vector name,
            /// size, precision, platform, implementation
            /// @param out - output stream
            void print ( std::ostream &out = std::cout ) const;

            /// Return the implementation method of the Blas Level 1 routines
            /// @return The implementation method of the Blas Level 1 routines
            enum IMPLEMENTATION get_implementation ( void ) const;

            /// Return the implementation name
            /// @return The implementation name
            std::string get_implementation_name ( void ) const;

            /// Return the platform of the vector
            /// @return The platform of the vector
            enum PLATFORM get_platform ( void ) const;

            /// Return the platform name
            /// @return The platform name
            std::string get_platform_name ( void ) const;

            /// Add a value to a specific index
            /// @param i - vector index
            /// @param val - the value
            virtual void add_value ( const int i, ValueType val ) = 0;

            /// Add a values to a specific indeces
            /// @param indices - vector indices
            /// @param length - number of indices
            /// @param values - the values
            virtual void add_values ( const int* indices, int length, const ValueType* values ) = 0;

            /// Set values of the vector
            /// @param index - index set (the memory should be
            ///                           allocated on the CPU)
            /// @param size - size of the index and values set
            /// @param values - values set (the memory should be
            ///                            allocated on the CPU)
            virtual void SetValues ( const int *index,
                                     const int size,
                                     const ValueType *values ) = 0;

            /// Get values of the vector
            /// @param index - index set (the memory should be
            ///                           allocated on the CPU)
            /// @param size - size of the index and values set
            /// @return values - return values set (the memory should be
            ///                           allocated on the CPU)
            virtual void GetValues ( const int *index,
                                     const int size,
                                     ValueType *values ) const = 0;

            /// Set block values of the vector
            /// @param start_i - beginning of the block
            /// @param end_i - ending of the block
            /// @param values - values set (the memory should be
            ///                             allocated on the CPU)
            virtual void SetBlockValues ( const int start_i,
                                          const int end_i,
                                          const ValueType *values ) = 0;

            /// Get block values of the vector
            /// @param start_i - beginning of the block
            /// @param end_i - ending of the block
            /// @return values - return values set (the memory should be
            ///                           allocated on the CPU)
            virtual void GetBlockValues ( const int start_i,
                                          const int end_i,
                                          ValueType *values ) const = 0;

            /// Initialize the indexset of the vector
            /// @param indexset - the index set
            /// @param size - the size of the index set
            /// the memory should be allocated on the CPU
            virtual void set_indexset ( const int *indexset,
                                        const int size ) = 0;

            /// Get the indexset of the vector;
            /// the buffer should be allocated on the CPU;
            /// the size of the buffer should be = get_indexset_size();
            /// @return indexset - the index set
            virtual void get_indexset ( int *indexset ) const = 0;

            /// Return the indexset size
            /// @return the indexset size
            virtual int get_indexset_size ( void ) const = 0;

            /// Copy the Indexset from a vector
            /// @param vec - the input vector
            virtual void CopyFromIndexset ( const lVector<ValueType>& vec );

            /// Copy the indexed elements to the *values buffer
            /// @param values - the buffer/ptr where the elements
            /// associate with the indexes will be copied;
            /// the size of the values buffer should be the same as the
            /// size specified in the set_indexset();
            /// the memory should be allocated on the CPU
            virtual void GetIndexedValues ( ValueType *values ) const = 0;

            /// Copy the *values buffer to the indexed elements of the vector
            /// @param values - the buffer/ptr where the elements
            /// associate with the indexes will be copied;
            /// the size of the values buffer should be the same as the
            /// size specified in the set_indexset();
            /// the memory should be allocated on the CPU
            virtual void SetIndexedValues ( const ValueType *values ) = 0;

            /// Extract a sub vector. The returned vector
            /// has the same platform/implementation as the
            /// original one.
            /// @param start_i - beginning of the block
            /// @param end_i - ending of the block
            /// @return lVector<ValueType>
            virtual lVector<ValueType> *extract_subvector ( const int start_i,
                                                            const int end_i ) const = 0;

            /// Partial replace of a vector to the vector;
            /// this->buffer[start_i+i] = sub_vec->buffer[start_sub_vec+i];
            /// @param start_i - the starting index of the vector
            /// @param start_sub_vec - the starting index of the sub vector
            /// @param size - the size of the data which has to be copied
            /// @param sub_vec - the sub vector
            virtual void partial_replace_subvector ( const int start_i,
                                                     const int start_sub_vec,
                                                     const int size,
                                                     const lVector<ValueType> &sub_vec ) = 0;

            /// Partial add of a weighted values of the sub_vector
            /// to the vector;
            /// this->buffer[start_i+i] += weight*sub_vec->buffer[start_sub_vec+i];
            /// @param start_i - the starting index of the vector
            /// @param start_sub_vec - the starting index of the sub vector
            /// @param size - the size of the data which has to be copied
            /// @param weight - the weight
            /// @param sub_vec - the sub vector
            virtual void partial_add_subvector ( const int start_i,
                                                 const int start_sub_vec,
                                                 const int size,
                                                 const ValueType weight,
                                                 const lVector<ValueType> &sub_vec ) = 0;

            /// Initialize (i.e. allocate) the vector
            /// @param size - size of the vector
            /// @param name - name of the vector
            virtual void Init ( const int size,
                                const std::string name ) = 0;

            /// Clear (i.e. free) the vector
            virtual void Clear ( void ) = 0;

            /// Set the vector to zero
            virtual void Zeros ( void ) = 0;

            /// Read data from ASCII file
            virtual void ReadFile ( const char* filename );

            /// Write data to ASCII file
            virtual void WriteFile ( const char* filename );

            /// Reordering the vector with respect to a dest permutation index set
            /// @param index - the dest index permutation vector
            virtual void Reorder ( const int *index ) = 0;

            /// See CopyTo() and CopyFrom()
            virtual lVector<ValueType> &operator= ( const lVector<ValueType> &vec2 ) = 0;

            /// Copy from and Copy to operators - can be use for different platforms
            virtual void CopyFrom ( const lVector<ValueType>& vec ) = 0;

            virtual void CopyTo ( lVector<ValueType>& vec ) const = 0;

            /// copy the structure; at the moment there is nothing
            virtual void CopyStructureFrom ( const lVector<ValueType> &vec2 );

            /// Copy and casting from double vector
            /// @param vec - the input vector<double>
            virtual void CastFrom ( const lVector<double>& vec ) = 0;
            virtual void CastFrom ( const lVector<float>& vec ) = 0;

            virtual void CastTo ( lVector<double>& vec ) const = 0;
            virtual void CastTo ( lVector<float>& vec ) const = 0;

            // Clone a vector - clean(), init(), copystructurefrom(), copyfrom()
            virtual void CloneFrom ( const lVector<ValueType>& vec );

            /// Clone without content of the class - only initialization of the vector is done
            /// @return the same class as this
            virtual lVector<ValueType> *CloneWithoutContent ( void ) const;

            /// Synchronize (sync threads in CUDA, OpenCL, etc)
            virtual void Sync ( void ) const;

            /// Component wise vector-vector multiplication
            virtual void ElementWiseMult ( const lVector<ValueType> &vec ) const = 0;

            /// Return the index of the minimum absolute value element.
            /// The index is 0-based.
            /// @return The index of the minimum absolute value element
            virtual int ArgMin ( void ) const = 0;

            /// Return the index of the maximum absolute value element of a vector
            /// The index is 0-based.
            /// @return The index of the maximum absolute value element of a vector
            virtual int ArgMax ( void ) const = 0;

            /// Return the L1 norm
            /// @return The L1 norm
            virtual ValueType Norm1 ( void ) const = 0;

            /// Return the Euclidean norm (L2 norm)
            /// @return The Euclidean norm (L2 norm)
            virtual ValueType Norm2 ( void ) const = 0;

            /// Return the Maximum Norm
            /// @return The Maximum Norm
            virtual ValueType NormMax ( void ) const = 0;
            /// Perform dot product (scalar product)
            /// @param vec - the second vector for the dot product
            /// @return the result of the dot product
            virtual ValueType Dot ( const lVector<ValueType> &vec ) const = 0;

            /// Scale the vector with scalar factor
            /// @param scalar - scalar factor
            virtual void Scale ( const ValueType scalar ) = 0;

            /// Update the vector (this = this + scalar*vec)
            /// @param vec - second vector
            /// @param scalar - scalar factor
            virtual void Axpy ( const lVector<ValueType> &vec,
                                const ValueType scalar ) = 0;

            /// Update the vector (this = scalar*this + vec)
            /// @param scalar - scalar factor
            /// @param vec - second vector
            virtual void ScaleAdd ( const ValueType scalar,
                                    const lVector<ValueType> &vec ) = 0;

            /// Return the sum of the absolute values of the elements in the vector
            /// @return sum of the absolute values of the elements
            //  virtual ValueType Sum (void) const = 0 ;

            /// Plane rotation
            virtual void Rot ( lVector<ValueType> *vec,
                               const ValueType &sc,
                               const ValueType &ss ) = 0;

            /// Returns the parameters for a Givens rotation
            virtual void Rotg ( ValueType *sa,
                                ValueType *sb,
                                ValueType *sc,
                                ValueType *ss ) const = 0;

            /// Returns rotation of points in the modified plane
            virtual void Rotm ( lVector<ValueType> *vec,
                                const ValueType &sparam ) = 0;

            /// Returns parameters for a modified Givens rotation
            virtual void Rotmg ( ValueType *sd1,
                                 ValueType *sd2,
                                 ValueType *x1,
                                 const ValueType &x2,
                                 ValueType *sparam ) const = 0;

            /// Returns entries of lvector

            virtual ValueType* GetBuffer ( ) const
            {
                return 0;
            }

          protected:

            int size_; // size of the vector

            std::string name_; // name of the vector
            std::string platform_name_; // platform name
            std::string implementation_name_; // implementation specification
            enum IMPLEMENTATION implementation_id_; //  implementation ID
            enum PLATFORM platform_id_; // platform ID

        };

        template <typename ValueType>
        inline int lVector<ValueType>::get_size ( void ) const
        {
            return this->size_;
        }

    } // namespace la
} // namespace hiflow

#endif
