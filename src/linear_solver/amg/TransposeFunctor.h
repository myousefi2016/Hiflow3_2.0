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

#ifndef SRC_LINEAR_SOLVER_AMG_TRANSPOSEFUNCTOR_H_
#    define SRC_LINEAR_SOLVER_AMG_TRANSPOSEFUNCTOR_H_

namespace hiflow
{
    namespace AMG
    {

        /// For local matrices

        template <class MatrixType>
        struct TransposeFunctor
        {

            void operator() ( MatrixType &matrix ) const
            {
                matrix.transpose_me ( );
            }
        };

        /// For coupled matrices

        template <class T>
        struct TransposeFunctor<hiflow::la::CoupledMatrix<T> >
        {

            void operator() ( hiflow::la::CoupledMatrix<T> &matrix ) const
            {
                // TODO
            }
        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_TRANSPOSEFUNCTOR_H_ */
