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

#include <cassert>

#include "polycoeff.h"

using namespace std;

namespace hiflow
{

    /// Constructor

    template<class T>
    PolyCoeff<T>::PolyCoeff ( ) : vector<T>( )
    {
    }

    template<class T>
    PolyCoeff<T>::PolyCoeff ( int k ) : vector<T>( k )
    {
        reinit ( k );
    }

    /// reinit

    template<class T>
    void PolyCoeff<T>::reinit ( int deg )
    {
        if ( !empty ( ) ) this->erase ( begin ( ), end ( ) );

        this->reserve ( deg + 1 );
        this->insert ( begin ( ), deg + 1, T ( ) );

        assert ( size ( ) == deg + 1 );
    }

    /// template instanciation

    template class PolyCoeff<double>;

} // namespace hiflow
