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

/// \author Martin Baumann

#ifndef __POLYNOMIALS_POLYCOEFF_H_
#    define __POLYNOMIALS_POLYCOEFF_H_

#    include <vector>

namespace hiflow
{

    template<class T>
    class PolyCoeff : public std::vector<T>
    {
      public:

        PolyCoeff ( );
        PolyCoeff ( int );

        void reinit ( int );

        using std::vector<T>::begin;
        using std::vector<T>::clear;
        using std::vector<T>::end;
        using std::vector<T>::empty;
        using std::vector<T>::reserve;
        using std::vector<T>::size;

    };

    /* --------------------------------------------------------------
    /                                                               /
    /  INLINE FUNCTIONS for PolyCoeff                               /
    /                                                               /
    /  ----------------------------------------------------------- */

} // namespace hiflow

#endif
