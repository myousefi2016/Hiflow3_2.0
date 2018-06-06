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

///
/// \file write_csv.h
/// \brief Function to write data to .csv files.
///
/// \author Simon Gawlok
///

#ifndef _WRITE_CSV_H_
#    define _WRITE_CSV_H_

#    include <fstream>
#    include <iostream>
#    include <utility>
#    include <string>
#    include <vector>
#    include <assert.h>
using namespace std;

void write_csv ( std::vector<std::string> names, std::vector<std::vector<double> > data, std::string filename );

#endif
