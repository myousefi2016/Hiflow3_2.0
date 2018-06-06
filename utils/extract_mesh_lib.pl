#!/usr/bin/perl -w

# Copyright (C) 2011-2017 Vincent Heuveline
#
# HiFlow3 is free software: you can redistribute it and/or modify it under the
# terms of the European Union Public Licence (EUPL) v1.2 as published by the
#/ European Union or (at your option) any later version.
#
# HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
# details.
#
# You should have received a copy of the European Union Public Licence (EUPL) v1.2
# along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

# Perl script to extract mesh module from HiFlow directory structure.
# author: Thomas Gengenbach

use strict;

my $my_actual_path = `pwd`;
chomp($my_actual_path);
print "$my_actual_path\n";
my @full_hiflow_path = split(/\//, $my_actual_path);

my $utils_path = pop(@full_hiflow_path);
my $hiflow_path = pop(@full_hiflow_path);

# check if we are in the correct path, otherwise DIE.
unless ($utils_path eq "utils" && $hiflow_path eq "hiflow") {
    die "Run this script only in the utils directory of hiflow...\n";
}

# change to directory on the same level as hiflow
chdir("../../");

my $libmesh_dir = "elemesh";

# dies if libmesh already exists
system("mkdir $libmesh_dir");
system("mkdir $libmesh_dir/src/");

print "\nCreated new directory $my_actual_path/../../$libmesh_dir\n\n";

# copy mesh directory to newly created directory in home
print ("cp -r $my_actual_path/../src/mesh/ $libmesh_dir/src/\n");
system("cp -r $my_actual_path/../src/mesh/ $libmesh_dir/src/");

# copy common dir
print ("cp -r $my_actual_path/../src/common/ $libmesh_dir/src/common/\n");
system("cp -r $my_actual_path/../src/common/ $libmesh_dir/src/common/");

# copy all tests
print ("cp -r $my_actual_path/../test/ $libmesh_dir/\n");
system("cp -r $my_actual_path/../test/ $libmesh_dir/");

# copy contrib directory (tinyxml)
print ("cp -r $my_actual_path/../contrib/ $libmesh_dir/\n");
system("cp -r $my_actual_path/../contrib/ $libmesh_dir/");

# copy documentation
print ("cp -r $my_actual_path/../doc/ $libmesh_dir/\n");
system("cp -r $my_actual_path/../doc/ $libmesh_dir/");

# copy CMakeLists.txt of src
print ("cp -r $my_actual_path/../src/CMakeLists.txt $libmesh_dir/src/\n");
system("cp -r $my_actual_path/../src/CMakeLists.txt $libmesh_dir/src/");

# go to newly created libmesh dir
chdir($libmesh_dir);

# change CMakeLists.txt in test
print "\nChange original test CMakeLists.txt to new one that only contains mesh specific tests.\n";
my $cmake_test_orig = "test/CMakeLists.txt";
my $cmake_test_new = "test/CMakeLists";
open(CMAKE_TEST_ORIG, $cmake_test_orig);
open(CMAKE_TEST_NEW, ">$cmake_test_new");
while(<CMAKE_TEST_ORIG>) {
    if ((
        (/include_directories/ ||
         /link_directories/) && (/mesh/                ||
                                 /Boost/               ||
                                 /common/              ||
                                 /src\/utils/          ||
                                 /MPI_INCLUDE_PATH/ )) ||
        /Common defines/ ||
        /add_definitions/ ||
        /DEMOS/ ||
        ((/add_executable/        ||
          /target_link_libraries/ ||
          /add_test/) && (/partitioned_lung_demo/             ||
                          /db_synchronize_demo/               ||
                          /iteration_demo/                    ||
                          /advanced_iteration_demo/           ||
                          /refinement_demo/                   ||
                          /vertex_sync_demo/                  ||
                          /cell_type_demo/                    ||
                          /uniform_refinement_perf_test/      ||
                          /visualization_of_random_data_demo/ ||
                          /partitioning_demo/                 ||
                          /global_partitioning_demo/          ||
                          /meshdb_test/                       ||
                          /reader_test/                       ||
                          /vtk_reader_test/                   ||
                          /parallel_vtk_reader_test/          ||
                          /big_mesh_test/                     ||
                          /partitioned_lung_test/             ||
                          /writer_test/                       ||
                          /incidence_test/                    ||
                          /vertex_lookup_test/                ||
                          /split_mesh_test/                   ||
                          /refinement_test/                   ||
                          /tetra_refinement_test/             ||
#                         /boundary_extraction_test/          ||
                          /material_number_test/              ||
                          /arbitrary_data_test/
         ))
        ||
        /WITH_METIS/
        ||
        /TESTS/
        ||
        /mesh module tests/
        )
    {
        print CMAKE_TEST_NEW $_;
    }
}
close(CMAKE_TEST_NEW);
close(CMAKE_TEST_ORIG);
rename($cmake_test_new, "$cmake_test_new.txt");

# change CMakeLists.txt in src
print "\nChange original src CMakeLists.txt to new one that only contains mesh specific srcs.\n";
my $cmake_src_orig = "src/CMakeLists.txt";
my $cmake_src_new = "src/CMakeLists";
open(CMAKE_SRC_ORIG, $cmake_src_orig);
open(CMAKE_SRC_NEW, ">$cmake_src_new");
while(<CMAKE_SRC_ORIG>) {
    if (/add_subdirectory/ && (/mesh/ ||
                               /common/ ))
    {
        print CMAKE_SRC_NEW $_;
    }
}

print CMAKE_SRC_NEW "set(DOXYGEN_DIRS common mesh)

foreach(dir in \${DOXYGEN_DIRS})
  list(APPEND DOXYGEN_SRC_DIRS_SRC \${PROJECT_SOURCE_DIR}/src/\${dir})
endforeach(dir in DOXYGEN_DIRS)

set(DOXYGEN_SRC_DIRS_SRC \${DOXYGEN_SRC_DIRS_SRC} PARENT_SCOPE)";

close(CMAKE_SRC_NEW);
close(CMAKE_SRC_ORIG);
rename($cmake_src_new, "$cmake_src_new.txt");

# create new CMakeLists.txt file and write it.
my $cmake_list_filename = "CMakeLists.txt";
open(CMAKE_FILE_SRC, ">$cmake_list_filename");
print "Create $cmake_list_filename in $libmesh_dir\n";
print CMAKE_FILE_SRC "cmake_minimum_required(VERSION 2.6.2)

project(libmesh)

include(CTest)
enable_testing()

message(\"Build type = \${CMAKE_BUILD_TYPE}\")

find_package(Boost)
find_package(MPI)

\#\#\# OPTIONAL EXTERNAL COMPONENTS \#\#\#

\# check for METIS if it is desired

set(WITH_METIS 0)
option(WANT_METIS \"Compile mesh module with METIS partitioning.\")
if (WANT_METIS)
  find_library(METIS_LIB NAMES metis)
  if (METIS_LIB)
    set(WITH_METIS 1)
  endif (METIS_LIB)
endif (WANT_METIS)
list(APPEND OPTIONS_FLAGS \"-DWITH_METIS=\${WITH_METIS}\")

option(RUN_TEDIOUS_TESTS \"Run all tests, including those that take a long time.\")

# set up variables

# installation directories
set(HIFLOW_LIB_DIR lib)
set(HIFLOW_BIN_DIR bin)

# DEBUG / RELEASE build flags
set(HIFLOW_COMMON_FLAGS \"-g -pg\")
set(HIFLOW_PEDANTIC_WARNINGS \" -Wall -Woverloaded-virtual -Wreorder -Wnon-virtual-dtor\")
if(\"\${CMAKE_BUILD_TYPE}\" STREQUAL \"Debug\")
  set(HIFLOW_CUSTOM_FLAGS \"-O0 -fno-inline \${HIFLOW_PEDANTIC_WARNINGS}\")
else(\"\${CMAKE_BUILD_TYPE}\" STREQUAL \"Debug\")
  set(HIFLOW_CUSTOM_FLAGS \"\")
endif(\"\${CMAKE_BUILD_TYPE}\" STREQUAL \"Debug\")

set(CMAKE_CXX_FLAGS \"\${CMAKE_CXX_FLAGS} \${HIFLOW_COMMON_FLAGS} \${HIFLOW_CUSTOM_FLAGS} \${OPTIONS_FLAGS} \${MPI_LINK_FLAGS}\")

message(\"Compiler flags = \${CMAKE_CXX_FLAGS}\")

\# process subdirectories
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(contrib)

set(DOXYGEN_SRC_DIRS \${DOXYGEN_SRC_DIRS_SRC})
add_subdirectory(doc)";

close(CMAKE_FILE_SRC);

my $header =
"=======================================================================

                      EleMESH - Elegant MESH

=======================================================================";
# Create README
open(README, ">README");
print README
"
$header

This is the mesh library created by Staffan Ronnas and Thomas
Gengenbach.

Some concepts and ideas are based on the paper \"Efficient
Representation of Computational Meshes\" by Anders Logg.

http://home.simula.no/~logg/pub/papers/Logg2009a.pdf

TO GET STARTED
~~~~~~~~~~~~~~

There are many classes in the library, but only a few are important
for using the existing functionality. For reading and writing meshes,
look at the Reader and Writer classes and their sub-classes. The Mesh
interface is used to interact with the mesh as a whole, while the
Entity class provides a view of just one entity (a cell, face, edge or
vertex). For looping over entities, the EntityIterator and
IncidentEntityIterator classes are used.

The documentation, which can be generated in HTML format with \"make
doxygen\" is still incomplete, but there are some descriptions for the
classes mentioned above. Another way to understand how the library
works is to look at the demo code in the test/ directory. Especially
iteration_demo.cc, refinement_demo.cc, reader_test.cc and
writer_test.cc are simple examples of the use of the library.

DEPENDENCIES
~~~~~~~~~~~~

This module depends heavily on MPI and METIS. MPI is required, at
least to make a good use of this library, METIS is optional but
recommended.

DISCLAIMER
~~~~~~~~~~

This library is still under heavy development and may (almost
certainly) contain bugs. Use at your own risk.

CONTACT
~~~~~~~

Please give us feedback on the usability of the library and report on
bugs.

Mailto: staffan.ronnas\@student.kit.edu and thomas.gengenbach\@kit.edu.

";

close(README);

# Create README
open(INSTALL, ">INSTALL");
print INSTALL
"
$header

INSTALLATION
~~~~~~~~~~~~

To install this module, you need to have cmake in version >2.6.2
installed. Then you follow the usual steps building a cmake-based
project: - create a build directory and go in there

         - type \"cmake path/to/source/\" or \"ccmake path/to/source\"
           if you prefer it graphically.

         - type \"make\"

         - type \"make test\" to run all the provided tests (more
           tests if you activate the RUN_TEDIOUS_TESTS switch in
           cmake)

         - type \"make doxygen\" to create the doxygen documentation.

";

close(INSTALL);

# create tarball
chdir("../");
print "\nArchiving $libmesh_dir...\n";
my $tarball = $libmesh_dir . ".tar.bz2";
my $tarcmd="tar -jpscf " . $tarball . " " . $libmesh_dir;
print "$tarcmd\n";
system($tarcmd);

#clean up
system("rm -r $libmesh_dir");
