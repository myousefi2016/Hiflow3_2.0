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

// 

#ifndef HIFLOW_COMMON_TENSOR_H
#    define HIFLOW_COMMON_TENSOR_H

#    include <cassert>

/// \author Staffan Ronnas
namespace hiflow
{

    template<int R>
    class TensorIndex
    {
      public:

        TensorIndex ( )
        {
        }

        TensorIndex ( const int* strides )
        {
            std::copy ( strides, strides + R + 1, &strides_[0] );
        }

        // Return stride for index r (0 <= r < R)

        int stride ( int r ) const
        {
            assert ( r >= 0 );
            assert ( r < R );
            return strides_[r + 1];
        }

        // Return dimension for index r (0 <= r < R)

        int dim ( int r ) const
        {
            return strides_[r] / strides_[r + 1];
        }

        int size ( ) const
        {
            return strides_[0];
        }

        int operator() ( int i ) const
        {
            return checked_index ( i * stride ( 0 ) );
        }

        int operator() ( int i, int j ) const
        {
            return checked_index ( i * stride ( 0 ) + j * stride ( 1 ) );
        }

        int operator() ( int i, int j, int k ) const
        {
            return checked_index ( i * stride ( 0 ) + j * stride ( 1 ) + k * stride ( 2 ) );
        }

        int operator() ( int i, int j, int k, int l ) const
        {
            return checked_index ( i * stride ( 0 ) + j * stride ( 1 ) + k * stride ( 2 ) + l * stride ( 3 ) );
        }

        int operator() ( int i, int j, int k, int l, int m ) const
        {
            return checked_index ( i * stride ( 0 ) + j * stride ( 1 ) + k * stride ( 2 ) + l * stride ( 3 ) + m * stride ( 4 ) );
        }

        int operator() ( int i, int j, int k, int l, int m, int n ) const
        {
            return checked_index ( i * stride ( 0 ) + j * stride ( 1 ) + k * stride ( 2 )
                                   + l * stride ( 3 ) + m * stride ( 4 ) + n * stride ( 5 ) );
        }

      private:

        int checked_index ( int i ) const
        {
            assert ( i >= 0 );
            assert ( i < strides_[0] );
            return i;
        }

        ///< Strides for ranks 0 ,..., R-1
        /// where strides_[0] is the total size of the array
        /// and strides_[R] = 1. By convention, the strides are sorted
        /// in decreasing order.
        int strides_[R + 1];
    };

    inline
    TensorIndex<0> make_tensor_index ( )
    {
        int s[1] = { 1 };
        return TensorIndex<0>( s );
    }

    inline
    TensorIndex<1> make_tensor_index ( int dim1 )
    {
        int s[2] = { dim1, 1 };
        return TensorIndex<1>( s );
    }

    inline
    TensorIndex<2> make_tensor_index ( int dim1, int dim2 )
    {
        int s[3] = { dim1 * dim2, dim2, 1 };
        return TensorIndex<2>( s );
    }

    inline
    TensorIndex<3> make_tensor_index ( int dim1, int dim2, int dim3 )
    {
        int s[4] = { dim1 * dim2 * dim3, dim2 * dim3, dim3, 1 };
        return TensorIndex<3>( s );
    }

    inline
    TensorIndex<4> make_tensor_index ( int dim1, int dim2, int dim3, int dim4 )
    {
        int s[5] = { dim1 * dim2 * dim3 * dim4, dim2 * dim3 * dim4, dim3 * dim4, dim4, 1 };
        return TensorIndex<4>( s );
    }

    inline
    TensorIndex<5> make_tensor_index ( int dim1, int dim2, int dim3, int dim4, int dim5 )
    {
        int s[6] = { dim1 * dim2 * dim3 * dim4 * dim5, dim2 * dim3 * dim4 * dim5, dim3 * dim4 * dim5, dim4 * dim5, dim5, 1 };
        return TensorIndex<5>( s );
    }

    inline
    TensorIndex<6> make_tensor_index ( int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 )
    {
        int s[7] = { dim1 * dim2 * dim3 * dim4 * dim5 * dim6,
                    dim2 * dim3 * dim4 * dim5 * dim6,
                    dim3 * dim4 * dim5 * dim6,
                    dim4 * dim5 * dim6,
                    dim5 * dim6,
                    dim6,
                    1 };
        return TensorIndex<6>( s );
    }

    template<int R>
    class Tensor
    {
      public:

        Tensor ( )
        : is_owner_ ( false ), values_ ( 0 )
        {
        }

        Tensor ( const TensorIndex<R>& index, double* values = 0, bool assume_ownership = false )
        : index_ ( index ), is_owner_ ( assume_ownership )
        {
            std::cout << "Index size = " << index_.size ( ) << "\n";

            if ( values )
            {
                if ( assume_ownership )
                {
                    values_ = values;
                }
                else
                {
                    values_ = new double[index_.size ( )];
                    std::copy ( values, values + index_.size ( ), values_ );
                }
            }
            else
            {
                values_ = new double[index_.size ( )];
                is_owner_ = true;
            }
        }

        ~Tensor ( )
        {
            clear ( );
        }

        Tensor ( const Tensor& t )
        {
            copy ( t );
        }

        const Tensor& operator= ( const Tensor& t )
        {
            if ( &t != this )
            {
                clear ( );
                copy ( t );
            }

            return *this;
        }

        // Access

        double& operator() ( int i )
        {
            return values_[index_ ( i )];
        }

        const double& operator() ( int i ) const
        {
            return values_[index_ ( i )];
        }

        double& operator() ( int i, int j )
        {
            return values_[index_ ( i, j )];
        }

        const double& operator() ( int i, int j ) const
        {
            return values_[index_ ( i, j )];
        }

        double& operator() ( int i, int j, int k )
        {
            return values_[index_ ( i, j, k )];
        }

        const double& operator() ( int i, int j, int k ) const
        {
            return values_[index_ ( i, j, k )];
        }

        double& operator() ( int i, int j, int k, int l )
        {
            return values_[index_ ( i, j, k, l )];
        }

        const double& operator() ( int i, int j, int k, int l ) const
        {
            return values_[index_ ( i, j, k, l )];
        }

        double& operator() ( int i, int j, int k, int l, int m )
        {
            return values_[index_ ( i, j, k, l, m )];
        }

        const double& operator() ( int i, int j, int k, int l, int m ) const
        {
            return values_[index_ ( i, j, k, l, m )];
        }

        double& operator() ( int i, int j, int k, int l, int m, int n )
        {
            return values_[index_ ( i, j, k, l, m, n )];
        }

        const double& operator() ( int i, int j, int k, int l, int m, int n ) const
        {
            return values_[index_ ( i, j, k, l, m, n )];
        }

        double* const as_array ( )
        {
            return values_;
        }

      private:

        void copy ( const Tensor& t )
        {
            index_ = t.index_;
            is_owner_ = t.is_owner_;

            if ( is_owner_ )
            {
                // deep copy
                values_ = new double[index_.size ( )];
                std::copy ( t.values_, t.values_ + index_.size ( ), values_ );
            }
            else
            {
                // shallow copy
                values_ = t.values_;
            }
        }

        void clear ( )
        {
            if ( is_owner_ && values_ )
            {
                delete[] values_;
                values_ = 0;
            }

            is_owner_ = false;
        }

        double* values_;
        TensorIndex<R> index_;
        bool is_owner_;
    };

}

#endif
