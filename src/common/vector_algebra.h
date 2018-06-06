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

#ifndef HIFLOW_VECTOR_ALGEBRA_H
#    define HIFLOW_VECTOR_ALGEBRA_H

#    include <cassert>
#    include <cmath>
#    include <vector>
#    include <iostream>

/// @brief This file contains template classes for representing small
/// floating-point vectors and matrices with sizes fixed at
/// compile-time; as well as common mathematical operations for these
/// objects.

/// @author Staffan Ronnas, Simon Gawlok

namespace
{
    // Tolerance for comparing elements of the vector.
    const double COMPARISON_TOL = 1.e-14;
}

namespace hiflow
{

    template<size_t M, size_t N, class DataType> class Mat;

    /// \brief Class representing a floating-point vector of size N.
    ///
    /// \details The class also supports common mathematical operations.

    template<size_t N, class DataType>
    class Vec
    {
      public:
        typedef DataType value_type;

        // Constructors

        Vec ( )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] = value_type ( 0. );
            }
        }

        Vec ( const Vec<N, DataType>& v )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] = v.v_[i];
            }
        }

        explicit Vec ( const value_type* values )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] = values[i];
            }
        }

        explicit Vec ( const std::vector<value_type> values )
        {
            assert ( values.size ( ) <= N );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                if ( values.size ( ) > i )
                    this->v_[i] = values[i];
                else
                    this->v_[i] = value_type ( 0. );
            }
        }

        // Assignment

        Vec<N, DataType>& operator= ( const Vec<N, DataType>& v )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] = v.v_[i];
            }
            return *this;
        }

        // Access operator

        inline value_type operator[] ( size_t i ) const
        {
            assert ( i >= 0 && i < N );
            return v_[i];
        }

        inline value_type& operator[] ( size_t i )
        {
            assert ( i >= 0 && i < N );
            return v_[i];
        }

        // Multiplication by scalar

        inline Vec<N, DataType>& operator*= ( const value_type s )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] = this->v_[i] * s;
            }
            return *this;
        }

        // Division by scalar

        inline Vec<N, DataType>& operator/= ( const value_type s )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] = this->v_[i] / s;
            }
            return *this;
        }

        // Addition

        inline Vec<N, DataType>& operator+= ( Vec<N, DataType> v )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] += v.v_[i];
            }
            return *this;
        }

        // Subtraction

        inline Vec<N, DataType>& operator-= ( Vec<N, DataType> v )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] -= v.v_[i];
            }
            return *this;
        }

        // Comparison

        inline bool operator== ( Vec<N, DataType> v )
        {
            for ( size_t i = 0; i < N; ++i )
            {
                if ( std::abs ( v_[i] - v[i] ) > static_cast < DataType > ( COMPARISON_TOL ) )
                {
                    return false;
                }
            }
            return true;
        }

        // Size

        inline size_t size ( ) const
        {
            return N;
        }

        // Add multiple of second vector

        inline void Axpy ( const Vec<N, DataType>& vec, const DataType alpha )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < N; ++i )
            {
                this->v_[i] += vec.v_[i] * alpha;
            }
        }

      private:
        value_type v_[N];
    };

    // Binary vector operations (do not modify the arguments).

    template<size_t N, class DataType>
    inline Vec<N, DataType> operator* ( const Vec<N, DataType>& v, const typename Vec<N, DataType>::value_type s )
    {
        Vec<N, DataType> tmp ( v );
        tmp *= s;
        return tmp;
    }

    template<size_t N, class DataType>
    inline Vec<N, DataType> operator* ( const typename Vec<N, DataType>::value_type s, const Vec<N, DataType>& v )
    {
        return v * s;
    }

    template<size_t N, class DataType>
    inline Vec<N, DataType> operator+ ( const Vec<N, DataType>& v1, const Vec<N, DataType>& v2 )
    {
        Vec<N, DataType> tmp ( v1 );
        tmp += v2;
        return tmp;
    }

    template<size_t N, class DataType>
    inline Vec<N, DataType> operator- ( const Vec<N, DataType>& v1, const Vec<N, DataType>& v2 )
    {
        Vec<N, DataType> tmp ( v1 );
        tmp -= v2;
        return tmp;
    }

    /// \brief Computes the scalar product of two vectors.

    template<size_t N, class DataType>
    inline typename Vec<N, DataType>::value_type dot ( const Vec<N, DataType>& v1, const Vec<N, DataType>& v2 )
    {
        typename Vec<N, DataType>::value_type res = typename Vec<N, DataType>::value_type ( );
        for ( size_t i = 0; i < N; ++i )
        {
            res += v1[i] * v2[i];
        }

        return res;
    }

    /// \brief Computes the cross product of two 3d vectors.

    template<class DataType>
    inline Vec<3, DataType> cross ( const Vec<3, DataType>& v1, const Vec<3, DataType>& v2 )
    {
        Vec<3, DataType> v3;
        v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
        return v3;
    }

    /// \brief Computes the sum of the vector components.

    template<size_t N, class DataType>
    inline typename Vec<N, DataType>::value_type sum ( const Vec<N, DataType>& v )
    {
        typename Vec<N, DataType>::value_type res = typename Vec<N, DataType>::value_type ( );
        for ( size_t i = 0; i < N; ++i )
        {
            res += v[i];
        }
        return res;
    }

    /// \brief Computes the Euclidean norm of a vector.

    template<size_t N, class DataType>
    inline typename Vec<N, DataType>::value_type norm ( const Vec<N, DataType>& v1 )
    {
        typename Vec<N, DataType>::value_type res = typename Vec<N, DataType>::value_type ( );
        res = std::sqrt ( dot ( v1, v1 ) );
        return res;
    }

    /// \brief Computes the normed normal of a 2d vector.

    template<class DataType>
    inline Vec<2, DataType> normal ( const Vec<2, DataType>& v )
    {
        Vec<2, DataType> v2;
        v2[0] = -v[1];
        v2[1] = v[0];
        DataType v2norm = norm ( v2 );
        v2 /= v2norm;
        return v2;
    }

    /// \brief Computes the normed normal of two 3d vectors.

    template<class DataType>
    inline Vec<3, DataType> normal ( const Vec<3, DataType>& v1, const Vec<3, DataType>& v2 )
    {
        Vec<3, DataType> v3 = cross ( v1, v2 );
        DataType v3norm = norm ( v3 );
        v3 /= v3norm;
        return v3;
    }

    /// \brief Computes the normed normal of two 3d vectors.

    template<class DataType>
    inline Vec<3, DataType> normal ( const Vec<6, DataType>& v1v2 )
    {
        DataType v1_array[3] = { v1v2[0], v1v2[1], v1v2[2] };
        DataType v2_array[3] = { v1v2[3], v1v2[4], v1v2[5] };
        Vec<3, DataType> v1 ( v1_array );
        Vec<3, DataType> v2 ( v2_array );
        Vec<3, DataType> v3 = cross ( v1, v2 );
        DataType v3norm = norm ( v3 );
        v3 /= v3norm;
        return v3;
    }

    /// \brief Output operator for vectors.

    template<size_t N, class DataType>
    std::ostream& operator<< ( std::ostream& os, const Vec<N, DataType>& v )
    {
        os << "[ ";
        for ( size_t i = 0; i < N; ++i )
        {
            os << v[i] << " ";
        }
        os << "]";
        return os;
    }

    /// \brief Class representing a floating-point matrix with M rows
    /// and N columns.

    template<size_t M, size_t N, class DataType>
    class Mat
    {
      public:

        typedef DataType value_type;

        // Constructors

        Mat ( )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] = value_type ( 0. );
            }
        }

        Mat ( const Mat<M, N, DataType>& m )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] = m.m_[i];
            }
        }

        explicit Mat ( const value_type* values )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] = values[i];
            }
        };

        inline Mat<M, N, DataType>& operator= ( const Mat<M, N, DataType>& m )
        {
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] = m.m_[i];
            }
            return *this;
        }

        // Element access

        inline value_type operator() ( size_t i, size_t j ) const
        {
            assert ( i >= 0 && i < M );
            assert ( j >= 0 && j < N );

            return m_[i * N + j];
        }

        inline value_type& operator() ( size_t i, size_t j )
        {
            assert ( i >= 0 && i < M );
            assert ( j >= 0 && j < N );
            return m_[i * N + j];
        }

        // Multiplication by scalar

        inline Mat<M, N, DataType>& operator*= ( const value_type s )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                m_[i] *= s;
            }
            return *this;
        }

        // Division by scalar

        inline Mat<M, N, DataType>& operator/= ( const value_type s )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] /= s;
            }
            return *this;
        }

        // Matrix addition

        inline Mat<M, N, DataType>& operator+= ( const Mat<M, N, DataType>& m )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] += m.m_[i];
            }
            return *this;
        }

        // Matrix subtraction

        inline Mat<M, N, DataType>& operator-= ( const Mat<M, N, DataType>& m )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] -= m.m_[i];
            }
            return *this;
        }

        // Matrix multiplication with square matrix

        inline Mat<M, N, DataType>& operator*= ( const Mat<N, N, DataType>& m )
        {
            // copy values of this to new array
            value_type cpy[M * N];
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                cpy[i] = this->m_[i];
                this->m_[i] = 0;
            }

            // perform matrix-matrix multiplication
            for ( size_t i = 0; i < M; ++i )
            { // loop over rows
                const size_t i_offset = i*N;
                for ( size_t k = 0; k < N; ++k )
                { // inner loop
                    const size_t k_offset = k*N;
#    pragma clang loop vectorize(enable)
                    for ( size_t j = 0; j < N; ++j )
                    { // loop over columns
                        this->m_[i_offset + j] += cpy[i_offset + k] * m.m_[k_offset + j];
                    }
                }
            }
            return *this;
        }

        // Vector multiplication

        inline void VectorMult ( const Vec<N, DataType>& in, Vec<M, DataType>& out ) const
        {
            for ( size_t i = 0, e_i = M; i < e_i; ++i )
            {
                out[i] = static_cast < DataType > ( 0 );
                const size_t i_ind = i*N;
#    pragma clang loop vectorize(enable)
                for ( size_t j = 0, e_j = N; j < e_j; ++j )
                {
                    out[i] += this->m_[i_ind + j] * in[j];
                }
            }
        }

        // Add multiple of second matrix

        inline void Axpy ( const Mat<M, N, DataType>& mat, const DataType alpha )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0, e_i = M * N; i < e_i; ++i )
            {
                this->m_[i] += mat.m_[i] * alpha;
            }
        }

      private:
        value_type m_[M*N];
    };

    // Specializations for (unused) cases M = 0 or N = 0.

    template<size_t M, class DataType>
    class Mat<M, 0, DataType>
    {
    };

    template<size_t N, class DataType>
    class Mat<0, N, DataType>
    {
    };

    // Matrix-Matrix multiplication: out = mat1 * mat2

    template<size_t M, size_t N, size_t P, class DataType>
    inline void MatrixMatrixMult ( Mat<M, P, DataType>& out, const Mat<M, N, DataType>& mat1, const Mat<N, P, DataType>& mat2 )
    {
        for ( size_t i = 0; i < M; ++i )
        { // loop over rows
            for ( size_t k = 0; k < N; ++k )
            { // inner loop
#    pragma clang loop vectorize(enable)
                for ( size_t j = 0; j < P; ++j )
                { // loop over columns
                    out ( i, j ) += mat1 ( i, k ) * mat2 ( k, j );
                }
            }
        }
    }

    // Binary matrix operations (do not modify the arguments).

    template<size_t M, size_t N, class DataType>
    inline Mat<M, N, DataType> operator* ( const Mat<M, N, DataType>& m, const typename Mat<M, N, DataType>::value_type s )
    {
        Mat<M, N, DataType> tmp ( m );
        tmp *= s;
        return tmp;
    }

    template<size_t M, size_t N, class DataType>
    inline Mat<M, N, DataType> operator* ( const typename Mat<M, N, DataType>::value_type s, const Mat<M, N, DataType>& m )
    {
        Mat<M, N, DataType> tmp ( m );
        tmp *= s;
        return tmp;
    }

    template<size_t M, size_t N, class DataType>
    inline Mat<M, N, DataType> operator/ ( const Mat<M, N, DataType>& m, typename Mat<M, N, DataType>::value_type s )
    {
        Mat<M, N, DataType> tmp ( m );
        tmp /= s;
        return tmp;
    }

    template<size_t M, size_t N, class DataType>
    inline Mat<M, N, DataType> operator+ ( const Mat<M, N, DataType>& m1, const Mat<M, N, DataType>& m2 )
    {
        Mat<M, N, DataType> tmp ( m1 );
        tmp += m2;
        return tmp;
    }

    template<size_t M, size_t N, class DataType>
    inline Mat<M, N, DataType> operator- ( const Mat<M, N, DataType>& m1, const Mat<M, N, DataType>& m2 )
    {
        Mat<M, N, DataType> tmp ( m1 );
        tmp -= m2;
        return tmp;
    }

    // Matrix-matrix multiplication with square matrix.

    template<size_t M, size_t N, class DataType>
    inline Mat<M, N, DataType> operator* ( const Mat<M, N, DataType>& m1, const Mat<N, N, DataType>& m2 )
    {
        Mat<M, N, DataType> tmp ( m1 );
        tmp *= m2;
        return tmp;
    }

    // General matrix-matrix multiplication.

    template<size_t M, size_t N, size_t P, class DataType>
    inline Mat<M, P, DataType> operator* ( const Mat<M, N, DataType>& A, const Mat<N, P, DataType>& B )
    {
        Mat<M, P, DataType> C;
        for ( size_t i = 0; i < M; ++i )
        {
            for ( size_t j = 0; j < N; ++j )
            {
                for ( size_t k = 0; k < P; ++k )
                {
                    C ( i, k ) += A ( i, j ) * B ( j, k );
                }
            }
        }
        return C;
    }

    // Matrix-vector multiplication

    template<size_t M, size_t N, class DataType>
    inline Vec<M, DataType> operator* ( const Mat<M, N, DataType>& m, const Vec<N, DataType>& v )
    {
        Vec<M, DataType> mv;
        for ( size_t i = 0; i < M; ++i )
        {
            for ( size_t j = 0; j < N; ++j )
            {
                mv[i] += m ( i, j ) * v[j];
            }
        }
        return mv;
    }

    // Vector-matrix multiplication

    template<size_t M, size_t N, class DataType>
    inline Vec<N, DataType> operator* ( const Vec<M, DataType>& v, const Mat<M, N, DataType>& m )
    {
        Vec<N, DataType> mv;
        for ( size_t i = 0; i < M; ++i )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t j = 0; j < N; ++j )
            {
                mv[j] += m ( i, j ) * v[i];
            }
        }
        return mv;
    }

    // Matrix minor computation (helper function for det and inv of 3x3 matrix)

    template<class DataType>
    inline typename Mat<3, 3, DataType>::value_type matrix_minor ( const Mat<3, 3, DataType>& m, size_t i, size_t j )
    {
        const int indI0 = ( i + 1 ) % 3;
        const int indI1 = ( i + 2 ) % 3;
        const int indJ0 = ( j + 1 ) % 3;
        const int indJ1 = ( j + 2 ) % 3;
        return m ( indI0, indJ0 ) * m ( indI1, indJ1 ) - m ( indI0, indJ1 ) * m ( indI1, indJ0 );
    }

    // Sign of sum of two integers.

    template<class DataType>
    inline typename Mat<3, 3, DataType>::value_type sign ( size_t i, size_t j )
    {
        return (i + j ) % 2 == 0 ? 1. : -1.;
    }

    // Determinant for 1x1 matrix.

    template<class DataType>
    inline DataType det ( const Mat<1, 1, DataType>& m )
    {
        return m ( 0, 0 );
    }

    ///\brief Determinant for 2x2 matrix.

    template<class DataType>
    inline DataType det ( const Mat<2, 2, DataType>& m )
    {
        return m ( 0, 0 ) * m ( 1, 1 ) - m ( 0, 1 ) * m ( 1, 0 );
    }

    ///\brief Determinant for 3x3 matrix.

    template<class DataType>
    inline DataType det ( const Mat<3, 3, DataType>& m )
    {
        return m ( 0, 0 ) * ( m ( 1, 1 ) * m ( 2, 2 ) - m ( 1, 2 ) * m ( 2, 1 ) )
                - m ( 0, 1 ) * ( m ( 1, 0 ) * m ( 2, 2 ) - m ( 1, 2 ) * m ( 2, 0 ) )
                + m ( 0, 2 ) * ( m ( 1, 0 ) * m ( 2, 1 ) - m ( 1, 1 ) * m ( 2, 0 ) );
    }

    /// \brief Inverse of 1x1 matrix.

    template<class DataType>
    inline void inv ( const Mat<1, 1, DataType>& m, Mat<1, 1, DataType>& m_inv )
    {
        m_inv ( 0, 0 ) = 1. / m ( 0, 0 );
    }

    /// \brief Inverse of 2x2 matrix.

    template<class DataType>
    inline void inv ( const Mat<2, 2, DataType>& m, Mat<2, 2, DataType>& m_inv )
    {
        typename Mat<2, 2, DataType>::value_type d = det ( m );
        assert ( d != 0. );
        m_inv ( 0, 0 ) = m ( 1, 1 ) / d;
        m_inv ( 0, 1 ) = -m ( 0, 1 ) / d;
        m_inv ( 1, 0 ) = -m ( 1, 0 ) / d;
        m_inv ( 1, 1 ) = m ( 0, 0 ) / d;
    }

    /// \brief Inverse of 3x3 matrix.

    template<class DataType>
    inline void inv ( const Mat<3, 3, DataType>& m, Mat<3, 3, DataType>& m_inv )
    {
        // compute determinant
        const DataType d = det ( m );
        assert ( d != static_cast < DataType > ( 0 ) );

        m_inv ( 0, 0 ) = ( m ( 1, 1 ) * m ( 2, 2 ) - m ( 1, 2 ) * m ( 2, 1 ) ) / d;
        m_inv ( 0, 1 ) = -( m ( 0, 1 ) * m ( 2, 2 ) - m ( 0, 2 ) * m ( 2, 1 ) ) / d;
        m_inv ( 0, 2 ) = ( m ( 0, 1 ) * m ( 1, 2 ) - m ( 0, 2 ) * m ( 1, 1 ) ) / d;

        m_inv ( 1, 0 ) = -( m ( 1, 0 ) * m ( 2, 2 ) - m ( 1, 2 ) * m ( 2, 0 ) ) / d;
        m_inv ( 1, 1 ) = ( m ( 0, 0 ) * m ( 2, 2 ) - m ( 0, 2 ) * m ( 2, 0 ) ) / d;
        m_inv ( 1, 2 ) = -( m ( 0, 0 ) * m ( 1, 2 ) - m ( 0, 2 ) * m ( 1, 0 ) ) / d;

        m_inv ( 2, 0 ) = ( m ( 1, 0 ) * m ( 2, 1 ) - m ( 1, 1 ) * m ( 2, 0 ) ) / d;
        m_inv ( 2, 1 ) = -( m ( 0, 0 ) * m ( 2, 1 ) - m ( 0, 1 ) * m ( 2, 0 ) ) / d;
        m_inv ( 2, 2 ) = ( m ( 0, 0 ) * m ( 1, 1 ) - m ( 0, 1 ) * m ( 1, 0 ) ) / d;
    }

    /// \brief Inverse-Transpose of 1x1 matrix.

    template<class DataType>
    inline void invTransp ( const Mat<1, 1, DataType>& m, Mat<1, 1, DataType>& m_inv )
    {
        m_inv ( 0, 0 ) = 1. / m ( 0, 0 );
    }

    /// \brief Inverse-Transpose of 2x2 matrix.

    template<class DataType>
    inline void invTransp ( const Mat<2, 2, DataType>& m, Mat<2, 2, DataType>& m_inv )
    {
        typename Mat<2, 2, DataType>::value_type d = det ( m );
        assert ( d != 0. );
        m_inv ( 0, 0 ) = m ( 1, 1 ) / d;
        m_inv ( 0, 1 ) = -m ( 1, 0 ) / d;
        m_inv ( 1, 0 ) = -m ( 0, 1 ) / d;
        m_inv ( 1, 1 ) = m ( 0, 0 ) / d;
    }

    /// \brief Inverse-Transpose of 3x3 matrix.

    template<class DataType>
    inline void invTransp ( const Mat<3, 3, DataType>& m, Mat<3, 3, DataType>& m_inv )
    {
        // copy into  inverse matrix
        // compute determinant
        const DataType d = det ( m );
        assert ( d != static_cast < DataType > ( 0 ) );

        m_inv ( 0, 0 ) = ( m ( 1, 1 ) * m ( 2, 2 ) - m ( 1, 2 ) * m ( 2, 1 ) ) / d;
        m_inv ( 0, 1 ) = -( m ( 1, 0 ) * m ( 2, 2 ) - m ( 1, 2 ) * m ( 2, 0 ) ) / d;
        m_inv ( 0, 2 ) = ( m ( 1, 0 ) * m ( 2, 1 ) - m ( 1, 1 ) * m ( 2, 0 ) ) / d;

        m_inv ( 1, 0 ) = -( m ( 0, 1 ) * m ( 2, 2 ) - m ( 0, 2 ) * m ( 2, 1 ) ) / d;
        m_inv ( 1, 1 ) = ( m ( 0, 0 ) * m ( 2, 2 ) - m ( 0, 2 ) * m ( 2, 0 ) ) / d;
        m_inv ( 1, 2 ) = -( m ( 0, 0 ) * m ( 2, 1 ) - m ( 0, 1 ) * m ( 2, 0 ) ) / d;

        m_inv ( 2, 0 ) = ( m ( 0, 1 ) * m ( 1, 2 ) - m ( 0, 2 ) * m ( 1, 1 ) ) / d;
        m_inv ( 2, 1 ) = -( m ( 0, 0 ) * m ( 1, 2 ) - m ( 0, 2 ) * m ( 1, 0 ) ) / d;
        m_inv ( 2, 2 ) = ( m ( 0, 0 ) * m ( 1, 1 ) - m ( 0, 1 ) * m ( 1, 0 ) ) / d;
    }

    /// \brief Gaussian elimination of a linear system.

    template<size_t M, size_t N, class DataType>
    inline bool gauss ( Mat<M, N, DataType>& mat, Vec<M, DataType>& vec )
    {
        // Gaussian elimination with pivoting of a linear system of equations
        // transforms the given Matrix and vector
        // if the submatrix m[1:M,1:M] is regular, the solution of m[1:M,1:M] x = v
        // ist stored into vec
        // if the matrix consists of m = [m[1:M,1:M] Id_M], the invers of m
        // is stored in the second half of the matrix

        // the current version needs to get equal or more columns than rows
        assert ( N >= M );

        for ( size_t i = 0; i < M; ++i )
        {
            // find pivot row
            int pivot = i;
            for ( size_t p = i; p < M; ++p )
            {
                if ( std::abs ( mat ( pivot, i ) ) < std::abs ( mat ( p, i ) ) )
                {
                    pivot = p;
                }
            }
            // check if system is solvable
            if ( std::abs ( mat ( pivot, i ) ) < COMPARISON_TOL )
            {
                return false;
            }
            // swap rows
            if ( pivot != i )
            {
                for ( size_t n = i; n < N; ++n )
                {
                    std::swap ( mat ( pivot, n ), mat ( i, n ) );
                }
                std::swap ( vec[pivot], vec[i] );
            }
            // reduce
            vec[i] /= mat ( i, i );
            for ( size_t n = N - 1; n >= i; --n )
            {
                mat ( i, n ) /= mat ( i, i );
            }
            // elimination forwards
            for ( size_t m = i + 1; m < M; ++m )
            {
                vec[m] -= vec[i] * mat ( m, i );
                for ( size_t n = N - 1; n >= i; --n )
                {
                    mat ( m, n ) -= mat ( i, n ) * mat ( m, i );
                }
            }
        }

        // elimination backwards
        for ( int i = M - 1; i > 0; --i )
        {
            for ( int m = i - 1; m >= 0; --m )
            {
                vec[m] -= vec[i] * mat ( m, i );
                for ( size_t n = N - 1; n >= i; --n )
                {
                    mat ( m, n ) -= mat ( i, n ) * mat ( m, i );
                }
            }
        }
        return true;
    }

    /// \brief Transpose of a general matrix.

    template<size_t M, size_t N, class DataType>
    inline void trans ( const Mat<M, N, DataType>& m, Mat<N, M, DataType>& m_trans )
    {
#    pragma clang loop vectorize(enable)
        for ( size_t i = 0; i < N; ++i )
        {
            for ( size_t j = 0; j < M; ++j )
            {
                m_trans ( i, j ) = m ( j, i );
            }
        }
    }

    /// \brief Computes the Frobenius product of two quadratic matrices.

    template<size_t N, class DataType>
    inline DataType frob ( const Mat<N, N, DataType>& m1, const Mat<N, N, DataType>& m2 )
    {

        DataType res = 0.0;
        assert ( N > 0 );
        for ( size_t i = 0; i < N; ++i )
        {
#    pragma clang loop vectorize(enable)
            for ( size_t k = 0; k < N; ++k )
            {
                res += m1 ( i, k ) * m2 ( i, k );
            }
        }
        return res;
    }

    /// \brief Trace of a quadratic matrix.

    template<size_t M, class DataType>
    inline DataType trace ( const Mat<M, M, DataType>& m )
    {
        DataType trace = 0;
#    pragma clang loop vectorize(enable)
        for ( size_t i = 0; i < M; ++i )
        {
            trace += m ( i, i );
        }
        return trace;
    }

    /// \brief Output operator for a general matrix.

    template<size_t M, size_t N, class DataType>
    std::ostream& operator<< ( std::ostream& os, const Mat<M, N, DataType>& m )
    {
        os << "[";
        for ( size_t i = 0; i < M; ++i )
        {
            for ( size_t j = 0; j < N; ++j )
            {
                os << " " << m ( i, j );
            }
            if ( i < M - 1 )
            {
                os << ";";
            }
        }
        os << " ]";
        return os;
    }

}

#endif /* _VECTOR_ALGEBRA_H_ */
