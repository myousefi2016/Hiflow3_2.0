set(assembler_SOURCES
  met_flow_assembler.cc
  met_flow_est_assembler.cc
  met_flow_convdiff_assembler.cc
  met_flow_incomp_assembler.cc
  met_flow_bous_assembler.cc
  scalar_prod_assembler.cc
  mass_matrix_assembler.cc
)

set(assembler_PUBLIC_HEADERS
  met_flow_assembler.h
  met_flow_est_assembler.h
  met_flow_convdiff_assembler.h
  met_flow_incomp_assembler.h
  met_flow_bous_assembler.h
  scalar_prod_assembler.h
  mass_matrix_assembler.h
)

include_directories(${METFLOW_BASE_DIR}/assembler/cyl)
include(${METFLOW_BASE_DIR}/assembler/cyl/HiFlowFiles.cmake)

foreach (i ${cyl_SOURCES})
  list(APPEND assembler_SOURCES "cyl/${i}")
endforeach()

foreach (i ${cyl_PUBLIC_HEADERS})
  list(APPEND assembler_PUBLIC_HEADERS "cyl/${i}")
endforeach()

include_directories(${METFLOW_BASE_DIR}/assembler/cart)
include(${METFLOW_BASE_DIR}/assembler/cart/HiFlowFiles.cmake)

foreach (i ${cart_SOURCES})
  list(APPEND assembler_SOURCES "cart/${i}")
endforeach()

foreach (i ${cart_PUBLIC_HEADERS})
  list(APPEND assembler_PUBLIC_HEADERS "cart/${i}")
endforeach()

include_directories(${METFLOW_BASE_DIR}/assembler/estimator)
include(${METFLOW_BASE_DIR}/assembler/estimator/HiFlowFiles.cmake)

foreach (i ${estimator_SOURCES})
  list(APPEND assembler_SOURCES "estimator/${i}")
endforeach()

foreach (i ${estimator_PUBLIC_HEADERS})
  list(APPEND assembler_PUBLIC_HEADERS "estimator/${i}")
endforeach()

