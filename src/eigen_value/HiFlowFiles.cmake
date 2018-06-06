set(eigen_value_SOURCES

)

set(eigen_value_PUBLIC_HEADERS
  eigen_solver.h 
)

if (WITH_SLEPC)
  list(APPEND eigen_value_SOURCES
    slepc_environment.cc
    slepc_eigen_solver.cc
    slepc_general_eps.cc	
  )
  list(APPEND eigen_value_PUBLIC_HEADERS
	slepc_environment.h    
    slepc_eigen_solver.h
	slepc_general_eps.h
  )
endif()

