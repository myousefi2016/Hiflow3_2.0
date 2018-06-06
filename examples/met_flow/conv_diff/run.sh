#!/bin/bash

WORKSPACE="/home/phil/Programme/Hiflow/hiflow_dev/hiflow/application/met_flow_wavetank_full"

#mpirun -np 2 valgrind --tool=memcheck ./tehd_exe ${WORKSPACE}
mpirun -np 2  ./tehd_exe ${WORKSPACE}
