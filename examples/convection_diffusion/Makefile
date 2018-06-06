CXX=mpicxx
HIFLOW_DIR=/usr/local

CXXFLAGS=-I$(HIFLOW_DIR)/include/hiflow3 -O3 -g -fopenmp -DBASEDIR=\"$(CURDIR)\" 

# stabilized_convdiff_tutorial target (linked to METIS)
# stabilized_convdiff_tutorial: stabilized_convdiff_tutorial.o
#	$(CXX) -o stabilized_convdiff_tutorial stabilized_convdiff_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow -lmetis

# stabilized_convdiff_tutorial target (not linked to METIS)
stabilized_convdiff_tutorial: stabilized_convdiff_tutorial.o
	$(CXX) -o stabilized_convdiff_tutorial stabilized_convdiff_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow

stabilized_convdiff_tutorial.o: stabilized_convdiff_tutorial.h

# instationary_convdiff_tutorial target (linked to METIS)
# instationary_convdiff_tutorial: instationary_convdiff_tutorial.o
#	$(CXX) -o instationary_convdiff_tutorial instationary_convdiff_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow -lmetis

# instationary_convdiff_tutorial target (not linked to METIS)
instationary_convdiff_tutorial: instationary_convdiff_tutorial.o
	$(CXX) -o instationary_convdiff_tutorial instationary_convdiff_tutorial.o $(CXXFLAGS) -L$(HIFLOW_DIR)/lib -lhiflow

instationary_convdiff_tutorial.o: instationary_convdiff_tutorial.h

all: stabilized_convdiff_tutorial instationary_convdiff_tutorial

clean:
	rm *.o stabilized_convdiff_tutorial instationary_convdiff_tutorial
