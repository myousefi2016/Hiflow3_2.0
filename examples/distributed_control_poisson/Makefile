# Use the MPI compiler wrapper to link MPI libraries.
CXX = mpicxx

# Base directory of HiFlow installation and ILU++ installation
HIFLOW_DIR = /usr/local
ILU_DIR = 

# Uncomment the line you need, whether you installed with ILU++ or not, and comment the other.
# NO ILU++:
CXXFLAGS = -I$(HIFLOW_DIR)/include/hiflow3 -O3 -g -fopenmp 
# WITH ILU++:
#CXXFLAGS = -I$(HIFLOW_DIR)/include/hiflow3 -O3 -g -fopenmp -I$(ILU_DIR)/lib

#Depending on wether you've installed HiFlow with metis and/or with ILU++, uncomment the lines you need for your installation
# setup in the following and comment the others. 

# NO metis, NO ILU++:
distributed_control_poisson_tutorial: distributed_control_poisson_tutorial.o
	$(CXX) -o distributed_control_poisson_tutorial distributed_control_poisson_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow 

# ONLY metis, NO ILU++:
#distributed_control_poisson_tutorial: distributed_control_poisson_tutorial.o
#	$(CXX) -o distributed_control_poisson_tutorial distributed_control_poisson_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow -lmetis

# NO metis, ONLY ILU++:
#distributed_control_poisson_tutorial: distributed_control_poisson_tutorial.o
#	$(CXX) -o distributed_control_poisson_tutorial distributed_control_poisson_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow -L$(ILU_DIR)   -liluplusplus-1.1

#  metis AND ILU++:
#distributed_control_poisson_tutorial: distributed_control_poisson_tutorial.o
#	$(CXX) -o distributed_control_poisson_tutorial distributed_control_poisson_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow -lmetis  -L$(ILU_DIR)   -liluplusplus-1.1

distributed_control_poisson_tutorial.o: distributed_control_poisson_tutorial.h

clean:
	rm *.o distributed_control_poisson_tutorial
