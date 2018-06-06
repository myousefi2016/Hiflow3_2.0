set(nonlinear_SOURCES
  newton.cc
  damping_armijo.cc
  forcing_eisenstat_walker.cc)

set(nonlinear_PUBLIC_HEADERS
  nonlinear_problem.h
  nonlinear_solver.h
  nonlinear_solver_creator.h
  nonlinear_solver_factory.h
  newton.h
  damping_strategy.h
  damping_armijo.h
  forcing_strategy.h
  forcing_eisenstat_walker.h
  )
