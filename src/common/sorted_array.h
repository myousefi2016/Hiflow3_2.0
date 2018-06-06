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

#ifndef HIFLOW_SORTED_ARRAY_H
#    define HIFLOW_SORTED_ARRAY_H

#    include <algorithm>
#    include <cassert>
#    include <functional>
#    include <vector>

/// @author Staffan Ronnas, Simon Gawlok

namespace hiflow
{

    /// \brief Array of sorted elements.
    ///
    /// \details Array of elements of type T, sorted according to some
    /// strict weak ordering. By default operator< for the type T is
    /// used. The requirement on T itself is that it has operator==
    /// defined. The class is a wrapper around std::vector<T>, and one
    /// of its uses is to provide a more memory-efficient alternative
    /// to a std::map or std::set. 
    /// 

    template < class T, class Cmp = typename std::less<T> >
    class SortedArray
    {
      public:
        typedef std::vector<T> array_type;

        // Typedefs exported from std::vector.
        typedef typename array_type::iterator iterator;
        typedef typename array_type::const_iterator const_iterator;
        typedef typename array_type::size_type size_type;
        typedef typename array_type::value_type value_type;
        typedef typename array_type::reference reference;
        typedef typename array_type::const_reference const_reference;

        /// \brief Create an empty array.
        /// \param[in] cmp   a comparison functor object.
        explicit SortedArray ( const Cmp& cmp = Cmp ( ) );

        /// \brief Create an array of a given size.
        /// \param[in] size    initial size of the array
        /// \param[in] value   default value for all elements of the array
        /// \param[in] cmp     comparison function object
        explicit SortedArray ( size_type size, const T& value = T ( ), const Cmp& cmp = Cmp ( ) );

        /// \brief Create an array from a range of input iterators.
        ///
        /// \details The array will be sorted after creation.
        ///
        /// \param[in] first   start of range
        /// \param[in] last    end of range
        /// \param[in] cmp     comparison functor object
        template <class InputIterator>
        SortedArray ( InputIterator first, InputIterator last, const Cmp& cmp = Cmp ( ) );

        /// \brief Create an array by copying a std::vector<T>.
        /// 
        /// \details The array will be sorted after creation.
        ///
        /// \param[in] x       input values
        /// \param[in] cmp     comparison functor object
        explicit SortedArray ( const array_type& x, const Cmp& cmp = Cmp ( ) );

        /// \brief Insert a value into the array. The array will be kept sorted.
        /// \details Time complexity = O(log size() ), where 
        ///
        /// \param[in] value   value to insert
        void insert ( const T& value );

        /// \brief Erase a value from the array.
        ///
        /// \details  This function works like std::vector<T>::erase().
        ///
        /// \param[in]  it   iterator pointing to element that should be erased
        /// \return an iterator pointing to the next element in the array, or array::end() if no such element exists.

        iterator erase ( iterator it )
        {
            // TODO(Staffan, Thomas): Create range version of this.
            assert ( it < data_.end ( ) );
            return data_.erase ( it ); // std::vector::erase will advance iterator
        }

        /// \brief Erase a value from the array.
        ///
        /// \param[in] pos  position of element to erase.
        /// \return an iterator pointing to the next element in the array, or array::end() if no such element exists.

        iterator erase ( int pos )
        {
            assert ( pos >= 0 );
            assert ( pos < size ( ) );
            return erase ( begin ( ) + pos );
        }

        /// \brief Search for a value in the array.
        ///
        /// \param[in]  value  value to search for
        /// \param[out] pos    position where value was found (optional output argument)
        /// \return true if @p value was found in the array, false otherwise
        bool find ( const T& value, int* pos ) const;

        /// \brief Search for a value in the array and insert it if it is not found.
        ///
        /// \param[in]  value  value to search and insert if necessary
        /// \return true if @p value already existed in the array; false if @p value was inserted in this function
        bool find_insert ( const T& value );

        /// \return iterator to beginning of array.

        iterator begin ( )
        {
            return data_.begin ( );
        }
        /// \return iterator to beginning of array.

        const_iterator begin ( ) const
        {
            return data_.begin ( );
        }
        /// \return iterator to end of array.

        iterator end ( )
        {
            return data_.end ( );
        }
        /// \return iterator to end of array.

        const_iterator end ( ) const
        {
            return data_.end ( );
        }

        /// \return size of the array.

        size_type size ( ) const
        {
            return data_.size ( );
        }
        /// \return whether the array is empty.

        bool empty ( ) const
        {
            return data_.empty ( );
        }

        /// \return Reference to data array

        const array_type& data ( ) const
        {
            return this->data_;
        }

        /// \return Reference to data array

        array_type& data ( )
        {
            return this->data_;
        }

        /// \brief Resize the array.
        ///
        /// \details The array will be sorted after resizing.
        ///
        /// \param[in] size  new size
        /// \param[in] val   value for newly created elements

        void resize ( size_type size, T val = T ( ) )
        {
            const int old_size = this->size ( );
            const T last_elem = back ( );
            data_.resize ( size, val );

            // re-sort if necessary
            if ( size > old_size && cmp_ ( val, last_elem ) )
            {
                std::stable_sort ( begin ( ), end ( ) );
            }
        }

        /// \brief Reserve memory for array.
        ///
        /// \details Memory for size elements will be reserved.
        ///
        /// \param[in] size  number of elements for which memory is reserved

        void reserve ( size_type size )
        {
            // check whether new size is really larger than current size
            if ( size > this->size ( ) )
            {
                data_.reserve ( size );
            }
        }

        /// \brief Access an element of the array.  
        /// 
        /// \details The array is not kept sorted when array element
        /// is modified. Subsequent calls to find and insert may fail
        /// or return incorrect results until the array is sorted
        /// again.
        ///
        /// \param[in] i     index of element.
        /// \return mutable reference to element i.

        reference operator[] ( size_type i )
        {
            /* TODO range_check */ return data_[i];
        }

        /// \brief Access an element of the array.  
        /// 
        /// \param[in] i     index of element.
        /// \return immutable reference to element i.

        const_reference operator[] ( size_type i ) const
        {
            /* TODO range_check */ return data_[i];
        }

        /// \return mutable reference to first element of the array

        reference front ( )
        {
            return data_.front ( );
        }

        /// \return immutable reference to first element of the array

        const_reference front ( ) const
        {
            return data_.front ( );
        }

        /// \return mutable reference to last element of the array

        reference back ( )
        {
            return data_.back ( );
        }

        /// \return immutable reference to last element of the array

        const_reference back ( ) const
        {
            return data_.back ( );
        }

        /// \brief Swap data with a std::vector<T>
        ///
        /// \details The array will be sorted after the swap.
        ///
        /// \param[in] vec    the vector to swap data with.

        void swap ( array_type& vec )
        {
            data_.swap ( vec );
            std::stable_sort ( begin ( ), end ( ) );
        }

        /// \brief Swap data with another SortedArray.
        /// \param[in] array   the array to swap data with

        void swap ( SortedArray<T, Cmp>& array )
        {
            data_.swap ( array.data_ );
        }

        /// \brief Clear the array.

        void clear ( )
        {
            data_.clear ( );
        }

        /// \return immutable reference to underlying std::vector<T>

        operator const array_type& ( ) const
        {
            return data_;
        }

      private:
        array_type data_;
        Cmp cmp_;
    };

    template <typename T, typename Cmp>
    SortedArray<T, Cmp>::SortedArray ( const Cmp& cmp )
    : cmp_ ( cmp )
    {
    }

    template <typename T, typename Cmp>
    SortedArray<T, Cmp>::SortedArray ( typename SortedArray<T, Cmp>::size_type size,
                                       const T& value, const Cmp& cmp )
    : data_ ( size, value ), cmp_ ( cmp )
    {
    }

    template <typename T, typename Cmp>
    template <class InputIterator>
    SortedArray<T, Cmp>::SortedArray ( InputIterator first, InputIterator last, const Cmp& cmp )
    : data_ ( first, last ), cmp_ ( cmp )
    {
        std::stable_sort ( data_.begin ( ), data_.end ( ), cmp_ );
    }

    template <typename T, typename Cmp>
    SortedArray<T, Cmp>::SortedArray ( const typename SortedArray<T, Cmp>::array_type& vec, const Cmp& cmp )
    : data_ ( vec ), cmp_ ( cmp )
    {
        std::stable_sort ( data_.begin ( ), data_.end ( ), cmp_ );
    }

    // Insertion

    template <typename T, typename Cmp>
    void SortedArray<T, Cmp>::insert ( const T& value )
    {

        // reserve more memory if needed
        if ( data_.size ( ) == data_.capacity ( ) )
        {
            data_.reserve ( static_cast < int > ( data_.capacity ( ) * 2.0 ) );
        }

        if ( data_.empty ( ) )
        {
            data_.push_back ( value );
        }
        else
        {
            // special binary search with only one cmp per iteration
            // http://en.wikipedia.org/wiki/Binary_search_algorithm
            int lower = 0;
            int upper = data_.size ( );

            do
            {
                const int pos = lower + static_cast < int > ( ( upper - lower ) / 2 );
                if ( cmp_ ( data_[pos], value ) )
                {
                    lower = pos + 1;
                }
                else
                {
                    upper = pos;
                }
            }
            while ( lower < upper );

            assert ( lower == upper );

            // Here lower == upper
            // TODO(Staffan): check that the is really correct
            assert ( lower >= 0 );
            assert ( lower <= data_.size ( ) );
            assert ( lower == 0 || cmp_ ( data_[lower - 1], value ) );
            assert ( lower == data_.size ( ) || cmp_ ( value, data_[lower] ) );

            data_.insert ( data_.begin ( ) + lower, value );
        }
    }

    // Lookup

    template <typename T, typename Cmp>
    bool SortedArray<T, Cmp>::find ( const T& value, int* pos ) const
    {
        if ( data_.empty ( ) )
        {
            return false;
        }

        int lower = 0;
        int upper = data_.size ( );

        // binary search, like in insert()
        do
        {
            const int cur = lower + static_cast < int > ( ( upper - lower ) / 2 );
            if ( cmp_ ( data_[cur], value ) )
            {
                lower = cur + 1;
            }
            else
            {
                upper = cur;
            }
        }
        while ( lower < upper );

        assert ( lower == upper );
        assert ( lower >= 0 );
        assert ( lower <= data_.size ( ) );

        // if lower != data_.size(), check the data_ array to see if
        // the value matches
        const bool found = lower < data_.size ( ) && ( data_[lower] == value );

        if ( found && pos != 0 )
        {
            *pos = lower;
        }
        return found;
    }

    // Combine Lookup and Insertion

    template <typename T, typename Cmp>
    bool SortedArray<T, Cmp>::find_insert ( const T& value )
    {

        // reserve more memory if needed
        if ( data_.size ( ) == data_.capacity ( ) )
        {
            data_.reserve ( static_cast < int > ( data_.capacity ( ) * 2.0 ) );
        }

        if ( data_.empty ( ) )
        {
            data_.push_back ( value );
            return false;
        }
        else
        {
            // special binary search with only one cmp per iteration
            // http://en.wikipedia.org/wiki/Binary_search_algorithm
            int lower = 0;
            int upper = data_.size ( );

            do
            {
                const int pos = lower + static_cast < int > ( ( upper - lower ) / 2 );
                if ( cmp_ ( data_[pos], value ) )
                {
                    lower = pos + 1;
                }
                else
                {
                    upper = pos;
                }
            }
            while ( lower < upper );

            assert ( lower == upper );

            // Here lower == upper
            // TODO(Staffan): check that the is really correct
            assert ( lower >= 0 );
            assert ( lower <= data_.size ( ) );

            // if lower != data_.size(), check the data_ array to see if
            // the value matches
            const bool found = lower < data_.size ( ) && ( data_[lower] == value );

            // if value does not already exists, insert it
            if ( !found )
            {
                assert ( lower == 0 || cmp_ ( data_[lower - 1], value ) );
                assert ( lower == data_.size ( ) || cmp_ ( value, data_[lower] ) );

                data_.insert ( data_.begin ( ) + lower, value );
            }

            return found;
        }

        return false;
    }
}

#endif 
