set(cart_SOURCES
  scalar_prod_cart_assembler.cc
  mass_matrix_cart_assembler.cc
)

set(cart_PUBLIC_HEADERS
  scalar_prod_cart_assembler.h
  mass_matrix_cart_assembler.h
)

include_directories(${METFLOW_BASE_DIR}/assembler/cart/primal)
include(${METFLOW_BASE_DIR}/assembler/cart/primal/HiFlowFiles.cmake)

foreach (i ${primal_SOURCES})
  list(APPEND cart_SOURCES "primal/${i}")
endforeach()

foreach (i ${primal_PUBLIC_HEADERS})
  list(APPEND cart_PUBLIC_HEADERS "primal/${i}")
endforeach()

include_directories(${METFLOW_BASE_DIR}/assembler/cart/dual)
include(${METFLOW_BASE_DIR}/assembler/cart/dual/HiFlowFiles.cmake)

foreach (i ${dual_SOURCES})
  list(APPEND cart_SOURCES "dual/${i}")
endforeach()

foreach (i ${dual_PUBLIC_HEADERS})
  list(APPEND cart_PUBLIC_HEADERS "dual/${i}")
endforeach()

include_directories(${METFLOW_BASE_DIR}/assembler/cart/estimator)
include(${METFLOW_BASE_DIR}/assembler/cart/estimator/HiFlowFiles.cmake)

foreach (i ${estimator_SOURCES})
  list(APPEND cart_SOURCES "estimator/${i}")
endforeach()

foreach (i ${estimator_PUBLIC_HEADERS})
  list(APPEND cart_PUBLIC_HEADERS "estimator/${i}")
endforeach()

