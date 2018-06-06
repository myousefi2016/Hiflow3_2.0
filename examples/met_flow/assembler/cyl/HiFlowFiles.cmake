set(cyl_SOURCES
  scalar_prod_cyl_assembler.cc
  mass_matrix_cyl_assembler.cc
)

set(cyl_PUBLIC_HEADERS
  scalar_prod_cyl_assembler.h
  mass_matrix_cyl_assembler.h
)

include_directories(${METFLOW_BASE_DIR}/assembler/cyl/primal)
include(${METFLOW_BASE_DIR}/assembler/cyl/primal/HiFlowFiles.cmake)

foreach (i ${primal_SOURCES})
  list(APPEND cyl_SOURCES "primal/${i}")
endforeach()

foreach (i ${primal_PUBLIC_HEADERS})
  list(APPEND cyl_PUBLIC_HEADERS "primal/${i}")
endforeach()

include_directories(${METFLOW_BASE_DIR}/assembler/cyl/dual)
include(${METFLOW_BASE_DIR}/assembler/cyl/dual/HiFlowFiles.cmake)

foreach (i ${dual_SOURCES})
  list(APPEND cyl_SOURCES "dual/${i}")
endforeach()

foreach (i ${dual_PUBLIC_HEADERS})
  list(APPEND cyl_PUBLIC_HEADERS "dual/${i}")
endforeach()

