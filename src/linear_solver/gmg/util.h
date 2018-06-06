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

/// \author Aksel Alpay

#ifndef UTIL_H
#    define UTIL_H

#    include <vector>
#    include <numeric>
#    include "boost/type_traits.hpp"

#    include "common/pointers.h"

#    if __GNUC__
#        if __GNUC__ < 4 || \
              (__GNUC__ == 4 && (__GNUC_MINOR__ < 7 ))
#            define PRE_GCC47_COMPAT
#        endif
#    endif

#    ifdef PRE_GCC47_COMPAT
#        define IMPORT_FROM_BASECLASS(baseclass, type) \
 typedef typename baseclass::type type
#    else
#        define IMPORT_FROM_BASECLASS(baseclass, type) \
 using typename baseclass::type
#    endif

#    define TYPE_FROM_CLASS(class, type) \
 typedef typename class::type type

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {
            namespace util
            {
                /// Gets access to the raw internal array of a std::vector
                /// @return The interal array of a std::vector, or NULL if the
                /// vector is empty.
                /// @param v The std::vector

                template<typename T>
                inline T* raw_array ( std::vector<T>& v )
                {
#    if __cplusplus >= 201103L
                    return v.data ( );
#    else
                    return hiflow::vec2ptr ( v );
#    endif
                }

                /// Gets access to the raw internal array of a std::vector
                /// @return The interal array of a std::vector, or NULL if the
                /// vector is empty.
                /// @param v The std::vector

                template<typename T>
                inline const T* raw_array ( const std::vector<T>& v )
                {
#    if __cplusplus >= 201103L
                    return v.data ( );
#    else
                    return hiflow::vec2ptr ( v );
#    endif
                }

                /// This class implements an ostream that will only write to the stream
                /// on the master level. Calls from other processes will be ignored.
                /// E.g. the following code will print "Hello world" only the process
                /// of rank 0 in MPI_COMM_WORLD:
                ///
                /// \code{.cpp}
                /// master_ostream stream(std::cout, 0);
                /// stream << "Hello World!\n";
                /// \endcode
                /*
                class master_ostream
                {
                public:
                  /// Initializes the object. Collective on MPI_COMM_WORLD.
                  /// @param ostr The output stream that shall be used. The reference
                  /// must remain valid throughout the existence of the master_ostream
                  /// object.
                  /// @param master_rank The rank of the process (in MPI_COMM_WORLD)
                  /// on which data shall be written.

                  master_ostream(std::ostream& ostr, int master_rank)
                  : ostr_(ostr), master_rank_(master_rank)
                  {
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
                  }

                  /// Conversion to std::ostream&

                  operator std::ostream& ()
                  {
                    return ostr_;
                  }

                  operator const std::ostream& () const
                  {
                    return ostr_;
                  }

                  /// Return the master rank.

                  int get_master_rank() const
                  {
                    return master_rank_;
                  }

                  /// Get the process rank

                  int get_rank() const
                  {
                    return rank_;
                  }

                  /// @return whether the calling process is the master process

                  inline bool is_master_process() const
                  {
                    return rank_ == master_rank_;
                  }

                  typedef std::ostream& (*io_manip_type)(std::ostream&);
                private:
                  std::ostream& ostr_;
                  int master_rank_;
                  int rank_;
                };

                /// Implements \c operator<<.
                /// Only calls from the process specified as \c master_process during
                /// the construction of the \c master_ostream object will lead to output
                /// being written to the output stream. Calls from other processes will
                /// be ignored. This call is not collective, but obviously at least
                /// the specified master process has to call the function if any
                /// output is to be written at all.

                template<class T>
                master_ostream& operator<<(master_ostream& ostr, const T& msg)
                {
                  if (ostr.get_rank() == ostr.get_master_rank())
                    (std::ostream&)ostr << msg;

                  return ostr;
                }

                /// This version enables the use of io manips

                master_ostream& operator<<(master_ostream& ostr,
                        master_ostream::io_manip_type io_manip)
                {
                  if (ostr.get_rank() == ostr.get_master_rank())
                    io_manip((std::ostream&)ostr);

                  return ostr;
                }*/

                /// Checks if a condition is not met, and if so, throws an exception
                /// of type \c Exception_Type. This type must have a constructor
                /// taking a const char* string to describe the exception.
                /// This "assert" function will not be optimized away in release builds,
                /// hence it can be used to easily check conditions in release configurations
                /// as well.
                /// @param statement The expression that shall be checked.
                /// @param failure_msg A message to describe the problem if \c statement
                /// evaluates to \c false.

                template<class Exception_Type>
                inline void exceptioned_assert ( bool statement, const char* failure_msg )
                {
                    if ( !statement )
                        throw Exception_Type ( failure_msg );
                }

                /// @return the size of a dimension of a simple, C-style multi dimensional
                /// array. Works up to three dimensional arrays.
                /// @param x the array
                /// @param dim the dimension of which the size is to be returned

                template <typename T, std::size_t N1>
                inline std::size_t c_array_size ( const T ( &x ) [N1], std::size_t dim )
                {
                    assert ( dim == 0 );
                    return N1;
                }

                template <typename T, std::size_t N1, std::size_t N2>
                inline std::size_t c_array_size ( const T ( &x ) [N1][N2], std::size_t dim )
                {
                    assert ( dim <= 1 );
                    if ( dim == 0 )
                        return N1;
                    else
                        return N2;
                }

                template <typename T, std::size_t N1, std::size_t N2, std::size_t N3>
                inline std::size_t c_array_size ( const T ( &x ) [N1][N2][N3], std::size_t dim )
                {
                    assert ( dim <= 2 );
                    if ( dim == 0 )
                        return N1;
                    else if ( dim == 1 )
                        return N2;
                    else if ( dim == 2 )
                        return N3;
                }

                /// A multi dimensional array

                template<typename T>
                class multi_array
                {
                  public:
                    typedef std::size_t size_type;
                    typedef std::size_t index_type;

                    typedef T* iterator;
                    typedef const T* const_iterator;

                    /// Construct empty array with no dimensions

                    multi_array ( )
                    : data_ ( NULL )
                    {
                    }

                    /// Construct multi dimensional array with the dimensions given
                    /// as a \c std::vector.
                    /// @param sizes Specifies the dimensions. \c sizes.size() is the number
                    /// of dimensions, and \c sizes[i] the extent in the i-th dimension.
                    /// E.g, to construct a 2x3 array, \c sizes has to contain the elements {2, 3}.

                    explicit multi_array ( const std::vector<size_type>& sizes )
                    : sizes_ ( sizes ), data_ ( NULL )
                    {
                        init ( );
                    }

                    /// Construct multi dimensional array with the dimensions given
                    /// as a simple stack-based C-style array.
                    /// @param sizes Specifies the dimensions. \c sizes.size() is the number
                    /// of dimensions, and \c sizes[i] the extent in the i-th dimension.
                    /// E.g, to construct a 2x3 array, \c sizes has to contain the elements {2, 3}.

                    template<size_type N>
                    explicit multi_array ( const size_type ( &sizes ) [N] )
                    : data_ ( NULL ), sizes_ ( sizes + 0, sizes + N )
                    {
                        init ( );
                    }

                    /// Construct multi dimensional array with the dimensions given
                    /// as a simple C-style array.
                    /// @param sizes Specifies the dimensions. \c sizes.size() is the number
                    /// of dimensions, and \c sizes[i] the extent in the i-th dimension.
                    /// E.g, to construct a 2x3 array, \c sizes has to contain the elements {2, 3}.
                    /// @param num_dimensions The number of elements of \c sizes and thus
                    /// the number of dimensions of the \c multi_array

                    multi_array ( const size_type* sizes, size_type num_dimensions )
                    : data_ ( NULL ), sizes_ ( sizes, sizes + num_dimensions )
                    {
                        init ( );
                    }

                    /// Construct two dimensional array
                    /// @param size_x The extent of the array in dimension 0
                    /// @param size_y The extent of the array in dimension 1

                    multi_array ( size_type size_x, size_type size_y )
                    : data_ ( NULL )
                    {
                        sizes_.reserve ( 2 );
                        sizes_.push_back ( size_x );
                        sizes_.push_back ( size_y );
                        init ( );
                    }

                    /// Construct three dimensional array
                    /// @param size_x The extent of the array in dimension 0
                    /// @param size_y The extent of the array in dimension 1
                    /// @param size_z The extent of the array in dimension 2

                    multi_array ( size_type size_x, size_type size_y, size_type size_z )
                    : data_ ( NULL )
                    {
                        sizes_.reserve ( 3 );
                        sizes_.push_back ( size_x );
                        sizes_.push_back ( size_y );
                        sizes_.push_back ( size_z );
                        init ( );
                    }

                    /// Copy Constructor. May Throw.

                    multi_array ( const multi_array<T>& other )
                    : sizes_ ( other.sizes_ ), buffer_size_ ( other.buffer_size_ ), data_ ( NULL ),
                    position_increments_ ( other.position_increments_ )
                    {
                        if ( sizes_.size ( ) != 0 )
                        {
                            init ( );
                            std::copy ( other.begin ( ), other.end ( ), data_ );
                        }
                    }

                    /// Assignment operator. Provides strong exception guarantee.

                    multi_array<T>& operator= ( multi_array<T> other )
                    {
                        // using the copy and swap idiom (note the by-value function
                        // parameter!) grants us a strong exception guarantee
                        swap ( *this, other );
                        return *this;
                    }

                    ~multi_array ( )
                    {
                        if ( data_ )
                            delete [] data_;
                    }

                    /// Swap two multi arrays. Their sizes do not have to equal.
                    /// @param a The first array
                    /// @param b The second array

                    friend void swap ( multi_array<T>& a, multi_array<T>& b )
                    {
                        using std::swap;

                        swap ( a.sizes_, b.sizes_ );
                        swap ( a.buffer_size_, b.buffer_size_ );
                        swap ( a.data_, b.data_ );
                        swap ( a.position_increments_, b.position_increments_ );
                    }

                    /// @return The extent of a dimension
                    /// @param dim The index of the dimension

                    size_type get_extent_of_dimension ( std::size_t dim ) const
                    {
                        assert ( dim < sizes_.size ( ) );
                        return sizes_[dim];
                    }

                    /// @return The dimension of the array

                    size_type get_dimension ( ) const
                    {
                        return sizes_.size ( );
                    }

                    /// @return An iterator type to the beginning of the array

                    iterator begin ( )
                    {
                        return data_;
                    }

                    /// @return An iterator type to the beginning of the array

                    const_iterator begin ( ) const
                    {
                        return data_;
                    }

                    /// @return An iterator type pointing to one element beyond the
                    /// last element of the array

                    iterator end ( )
                    {
                        return data_ + buffer_size_;
                    }

                    /// @return An iterator type pointing to one element beyond the
                    /// last element of the array

                    const_iterator end ( ) const
                    {
                        return data_ + buffer_size_;
                    }

                    /// @return The total number of elements in the array

                    size_type get_num_elements ( ) const
                    {
                        return buffer_size_;
                    }

                    /// @return The total number of elements in the array

                    size_type size ( ) const
                    {
                        return get_num_elements ( );
                    }

                    /// Access an element of the array
                    /// @param position Contains the indices of the element to look up
                    /// @return A reference to the specified element

                    T& operator[] ( const std::vector<index_type>& position )
                    {
                        assert ( position.size ( ) == get_dimension ( ) );
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position Contains the indices of the element to look up
                    /// @return A reference to the specified element

                    const T& operator[] ( const std::vector<index_type>& position ) const
                    {
                        assert ( position.size ( ) == get_dimension ( ) );
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position Contains the indices of the element to look up
                    /// @return A reference to the specified element

                    T& operator[] ( const std::vector<int>& position )
                    {
                        assert ( position.size ( ) == get_dimension ( ) );
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position Contains the indices of the element to look up
                    /// @return A reference to the specified element

                    const T& operator[] ( const std::vector<int>& position ) const
                    {
                        assert ( position.size ( ) == get_dimension ( ) );
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position A simple C-array containing the indices of the
                    /// element to look up
                    /// @return A reference to the specified element

                    template<size_type N>
                    T& operator[] ( const index_type ( &position ) [N] )
                    {
                        assert ( N == get_dimension ( ) );
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position A simple C-array containing the indices of the
                    /// element to look up
                    /// @return A reference to the specified element

                    template<size_type N>
                    const T& operator[] ( const index_type ( &position ) [N] ) const
                    {
                        assert ( N == get_dimension ( ) );
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position A pointer to asimple C-array containing the
                    /// indices of the element to look up
                    /// @return A reference to the specified element

                    T& operator[] ( const index_type* position )
                    {
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                    /// Access an element of the array
                    /// @param position A pointer to asimple C-array containing the
                    /// indices of the element to look up
                    /// @return A reference to the specified element

                    const T& operator[] ( const index_type* position ) const
                    {
                        assert ( data_ != NULL );

                        size_type pos = calculate_position ( position );
                        return data_[pos];
                    }

                  private:
                    /// Based on the indices of an element, calculates the position
                    /// of the element in the flat, one-dimensional data array.
                    /// @param position An object of a type offering operator[], that
                    /// stores the indices of the element
                    /// @return the position in the one dimensional data array.

                    template<typename Container>
                    inline size_type calculate_position ( const Container& position ) const
                    {
                        size_type pos = 0;
                        for ( size_type i = 0; i < get_dimension ( ); ++i )
                            pos += position[i] * position_increments_[i];

                        return pos;
                    }

                    /// Initializes the data array and the position increments of each dimension

                    void init ( )
                    {
                        assert ( sizes_.size ( ) != 0 );
                        for ( std::size_t i = 0; i < sizes_.size ( ); ++i )
                            assert ( sizes_[i] != 0 );

                        buffer_size_ = std::accumulate ( sizes_.begin ( ), sizes_.end ( ), 1,
                                                         std::multiplies<size_type>( ) );

                        data_ = new T [buffer_size_];

                        position_increments_.resize ( get_dimension ( ) );
                        position_increments_[0] = 1;

                        for ( std::size_t i = 1; i < sizes_.size ( ); ++i )
                        {
                            position_increments_[i] = position_increments_[i - 1] * sizes_[i - 1];
                        }
                    }

                    std::vector<size_type> sizes_;
                    std::vector<index_type> position_increments_;

                    size_type buffer_size_;

                    T* data_;
                };

                /// \c operator<< implementation for vectors, to simplify printing vectors.
                /// @param ostr The std::ostream to use
                /// @param data The vector to write to \c ostr

                template<typename T>
                std::ostream& operator<< ( std::ostream& ostr, const std::vector<T>& data )
                {
                    for ( typename std::vector<T>::const_iterator it = data.begin ( );
                          it != data.end ( ); ++it )
                        ostr << *it << " ";

                    return ostr;
                }

            } // namespace util
        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
