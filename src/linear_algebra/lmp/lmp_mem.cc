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

#include "lmp_mem.h"

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

template<typename ValueType>
void memcpyhost ( ValueType *dest, const ValueType *src, const size_t size )
{
    memcpy ( dest, src, sizeof (ValueType ) * size );
}

void memsethost ( void *p, int value, size_t size, const size_t size_type )
{
    memset ( p, value, size * size_type );
}

template void memcpyhost<double >( double *dest, const double *src, const size_t size );
template void memcpyhost<float >( float *dest, const float *src, const size_t size );
template void memcpyhost<int >( int *dest, const int *src, const size_t size );
template void memcpyhost<short int>( short int *dest, const short int *src, const size_t size );
