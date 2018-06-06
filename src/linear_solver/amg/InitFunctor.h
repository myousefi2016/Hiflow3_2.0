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

/// @author Bernd Doser
/// @date 2015-10-12

#ifndef SRC_LINEAR_SOLVER_AMG_INITFUNCTOR_H_
#    define SRC_LINEAR_SOLVER_AMG_INITFUNCTOR_H_

namespace hiflow
{
    namespace AMG
    {

        /// For local matrices

        template <class MatrixType>
        struct InitFunctor
        {

            void operator() ( MatrixType &matrix, int nnz, int num_row, int num_col,
                    std::string const& name, MatrixType const& ref_matrix ) const
            {
                matrix.Init ( nnz, num_row, num_col, name );
            }
        };

        /// For coupled matrices

        template <class T>
        struct InitFunctor<hiflow::la::CoupledMatrix<T> >
        {

            void operator() ( hiflow::la::CoupledMatrix<T> &matrix, int nnz, int num_row,
                    int num_col, std::string const& namee, hiflow::la::CoupledMatrix<T> const& ref_matrix ) const
            {
                matrix.Init (
                              ref_matrix.comm ( ),
                              ref_matrix.row_couplings ( ),
                              ref_matrix.col_couplings ( ),
                              ref_matrix.diagonal ( ).get_platform ( ),
                              ref_matrix.diagonal ( ).get_implementation ( ),
                              ref_matrix.diagonal ( ).get_matrix_format ( )
                              );
            }
        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_INITFUNCTOR_H_ */
