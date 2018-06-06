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

/// \author Thomas Gengenbach, Staffan Ronnas

#ifndef _TEST_H_
#    define _TEST_H_

#    include <cstdlib>

#    define TEST(x) {                                                       \
        if (!(x)) {                                                     \
            std::clog << "Test failed: [" << __FILE__ << ":" <<         \
                __LINE__ << "\n\n\t" << #x << "\n";                     \
            exit(-1); }}

#    define TEST_EQUAL(x, y) {                                              \
        if ((x) != (y)) {                                               \
                     std::clog << "Test failed: [" << __FILE__ << ":"   \
                     << __LINE__ << "\n\n\t" << #x << " == "            \
                     <<  #y << "\n" << #x << " = " << x << "; " \
                               << #y << " = " << y << "\n";     \
                     exit(-1); }}

#    define TEST_NOT_EQUAL(x, y) {                                              \
        if ((x) == (y)) {                                               \
                     std::clog << "Test failed: [" << __FILE__ << ":"   \
                     << __LINE__ << "\n\n\t" << #x << " != "            \
                     <<  #y << "\n" << #x << " = " << x << "; " \
                               << #y << " = " << y << "\n";     \
                     exit(-1); }}

#    define TEST_EQUAL_EPS(x, y, eps) {                                     \
        if (((x) - (y)) >= eps || ((x) - (y)) <= -eps) {                \
                     std::clog << "Test failed: [" << __FILE__ << ":"   \
                     << __LINE__ << "\n\n\t" << #x << " == "            \
                     <<  #y << "\n" << #x << " = " << x << "; " \
                               << #y << " = " << y << "\n";     \
                     exit(-1); }}

#    define TEST_LESS(x, y) {                                               \
        if ((x) >= (y)) {                                               \
            std::clog << "Test failed: [" << __FILE__ << ":"            \
                      << __LINE__ << "\n\n\t" << #x << " < "            \
                      <<  #y << "\n" << #x << " = " << x << "; "        \
                      << #y << " = " << y << "\n"                       \
                      << "ERROR\n";     \
            exit(-1); }}

#endif /* _TEST_H_ */
