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

#ifndef __LMATRIX_FORMATS_H
#    define __LMATRIX_FORMATS_H

namespace hiflow
{
    namespace la
    {

        /// @brief COO Matrix structure
        /// @author Dimitar Lukarski

        template <typename ValueType>
        struct COO_lMatrixType
        {
            int *row; // row index
            int *col; // col index
            ValueType *val; // values
        };

        /// @brief CSR Matrix structure
        /// @author Dimitar Lukarski

        template <typename ValueType>
        struct CSR_lMatrixType
        {
            int *row; // row pointer
            int *col; // col index
            ValueType *val; // values
        };

        /// @brief MCSR Matrix structure
        /// @author Dimitar Lukarski

        template <typename ValueType>
        struct MCSR_lMatrixType
        {
            int *row; // row pointer
            int *col; // col index
            ValueType *val; // values
            ValueType *diag; // diagonal values
        };

        /// @brief DIA Matrix structure
        /// @author Dimitar Lukarski

        template <typename ValueType>
        struct DIA_lMatrixType
        {
            // extra matrix info
            int num_diagonals; // number of diagonals

            // data
            int *offsets; // the offsets wrt the diagonal
            ValueType *data; // values
        };

        /// @brief ELL Matrix structure
        /// @author Dimitar Lukarski

        template <typename ValueType>
        struct ELL_lMatrixType
        {
            // extra matrix info
            int num_cols; // max number of column per row

            // data
            int *index; //  column index set
            ValueType *data; // values
        };

        /// @brief DENSE Matrix structure
        /// @author Niels Wegh

        template <typename ValueType>
        struct DENSE_lMatrixType
        {
            ValueType *val; // values
        };

    } // namespace la
} // namespace hiflow

#endif
