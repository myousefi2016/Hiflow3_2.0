set(general_SOURCES
  vorticity.cc
  write_csv.cc
  met_flow_main.cc
  goal_functional.cc   
)

set(general_PUBLIC_HEADERS
  write_csv.h
  vorticity.h
  my_newton.h
  met_flow_main.h
  goal_functional.h
  source_term.h
  convection_term.h
)

include_directories(${METFLOW_BASE_DIR}/general/cyl)
include(${METFLOW_BASE_DIR}/general/cyl/HiFlowFiles.cmake)

foreach (i ${cyl_SOURCES})
  list(APPEND general_SOURCES "cyl/${i}")
endforeach()

foreach (i ${cyl_PUBLIC_HEADERS})
  list(APPEND general_PUBLIC_HEADERS "cyl/${i}")
endforeach()

