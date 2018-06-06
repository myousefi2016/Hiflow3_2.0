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

/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski

#ifndef __OPENCL_UTILS_H
#    define __OPENCL_UTILS_H

#    include <iostream>
#    include <sys/stat.h>
#    include <fstream>
#    include <stdlib.h>

#    include "config.h"

#    include "../lmp_log.h"

#    ifdef WITH_OPENCL

#        ifdef __APPLE__
#            include <cl.h>
#        else
#            include <CL/cl.h>
#        endif

static void checkOpenCLErr ( cl_int err, const char * name )
{

    if ( err != CL_SUCCESS )
    {
        LOG_ERROR ( "OpenCL - " << name << "Error:" << err );
        exit ( -1 );
    }
}

static char *load_program_source ( const char *filename )
{

    struct stat statbuf;
    FILE *fh = NULL;
    char *source;

    fh = fopen ( filename, "r" );
    if ( fh == 0 ) return 0;

    stat ( filename, &statbuf );
    source = ( char * ) malloc ( statbuf.st_size + 1 );
    if ( fread ( source, statbuf.st_size, 1, fh ) != 1 )
    {
        LOG_ERROR ( "OpenCL - Error in reading kernel file" );
        exit ( -1 );
    }
    source[statbuf.st_size] = '\0';

    return source;
}
#    endif

#endif
