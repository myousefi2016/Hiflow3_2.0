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

#ifndef HIFLOW_MPITOOLS__
#    define HIFLOW_MPITOOLS__

///
/// \author Martin Wlotzka
///

#    include <mpi.h>

template <typename DataType>
struct mpi_data_type
{
};

template<> struct mpi_data_type<long double>
{

    static MPI_Datatype get_type ( )
    {
        return MPI_LONG_DOUBLE;
    }
};

template<> struct mpi_data_type<double>
{

    static MPI_Datatype get_type ( )
    {
        return MPI_DOUBLE;
    }
};

template<> struct mpi_data_type<float>
{

    static MPI_Datatype get_type ( )
    {
        return MPI_FLOAT;
    }
};

template<> struct mpi_data_type<int>
{

    static MPI_Datatype get_type ( )
    {
        return MPI_INT;
    }
};

template<> struct mpi_data_type<unsigned int>
{

    static MPI_Datatype get_type ( )
    {
        return MPI_UNSIGNED;
    }
};

template<> struct mpi_data_type<char>
{

    static MPI_Datatype get_type ( )
    {
        return MPI_CHAR;
    }
};

#endif
