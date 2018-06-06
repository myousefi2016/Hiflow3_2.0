# Use the MPI compiler wrapper to link MPI libraries.
CXX=mpicxx

# Base directories of HiFlow installation, Hiflow contrib directory and directory for external packages
HIFLOW_DIR  = /home/phil/Programme/Hiflow/hiflow_server/install
CONTRIB_DIR = /home/phil/Programme/Hiflow/hiflow_server/hiflow/contrib
HIFLOW_EXT  = /home/phil/Programme/Hiflow_external

PETSC_DIR   = $(HIFLOW_EXT)/petsc-3.7.3/install-debug-complex
SLEPC_DIR   = $(HIFLOW_EXT)/slepc-3.7.2-complex/install-debug
HYPRE_DIR   = $(HIFLOW_EXT)/hypre-2.11.1/src/hypre
METIS_DIR   = $(HIFLOW_EXT)/metis-5.1.0
ILUPP_DIR   = $(HIFLOW_EXT)/ILU++_1.1.1
BOOST_DIR   = $(CONTRIB_DIR)/boost_libraries

# don't edit below this line ---------------------------------------------------
# Compiler Flags
CXXFLAGS    = 	-I$(HIFLOW_DIR)/include/hiflow3 \
		-I$(BOOST_DIR) -I$(BOOST_DIR)/boost/tr1 \
		-I$(PETSC_DIR)/include \
		-I$(SLEPC_DIR)/include \
		-I$(ILUPP_DIR)/lib \
		-I$(METIS_DIR)/libmetis \
		-I$(HYPRE_DIR)/include \
		-O3 -g -fopenmp

# Linker Flags 		
SLEPC_LINK   = 	-Wl,-rpath,$(SLEPC_DIR)/lib -L$(SLEPC_DIR)/lib -lslepc -Wl,-rpath,$(SLEPC_DIR)/lib -L$(SLEPC_DIR)/lib 
#-lparpack -larpack
PETSC_LINK   = 	-L$(PETSC_DIR)/lib  -lpetsc -Wl,-rpath,$(PETSC_DIR)/lib -L$(PETSC_DIR)/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -L/usr/local/lib -lparmetis -lmetis -lscalapack -lumfpack -lklu -lcholmod -lbtf -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -lflapack -lfblas
GENERAL_LINK = 	-L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lgfortran -lm -lgfortran -lm -lquadmath -lm -lmpi_cxx -lstdc++ -L/usr/lib/openmpi/lib -lmpifort -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lhwloc -lgcc_s -lpthread -ldl 

LINK_FLAGS   = 	-L$(HIFLOW_DIR)/lib -lhiflow \
		-L$(HYPRE_DIR)/lib -lHYPRE \
		-L$(ILUPP_DIR) $(ILUPP_DIR)/libiluplusplus-1.1.1.a \
		-L$(METS_DIR) -L$(HIFLOW_DIR)/lib -lmetis $(SLEPC_LINK) $(PETSC_LINK) $(GENERAL_LINK)


include $(SLEPC_DIR)/lib/slepc/conf/slepc_common


# poisson_precond target (use if hiflow3 was installed without metis)
poisson_eigenvalue: poisson_eigenvalue.o
	$(CXX) -o poisson_eigenvalue poisson_eigenvalue.o $(CXXFLAGS) $(LINK_FLAGS)

# poisson_precond target (use if hiflow3 was installed with metis)
#poisson_precond: poisson_precond.o
#	$(CXX) -o poisson_precond poisson_precond.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow -lmetis


poisson_eigenvalue.o: poisson_eigenvalue.h




