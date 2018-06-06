set(dof_SOURCES
  degree_of_freedom.cc 
  dof_interpolation.cc 
  dof_interpolation_pattern.cc 
  fe_interface_pattern.cc 
  dof_partition.cc
  numbering_strategy.cc
  numbering_lagrange.cc)

set(dof_PUBLIC_HEADERS
  dof_fem_types.h
  degree_of_freedom.h 
  dof_interpolation.h 
  dof_interpolation_pattern.h 
  fe_interface_pattern.h 
  dof_partition.h
  numbering_strategy.h
  numbering_lagrange.h)
