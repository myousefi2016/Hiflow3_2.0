// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
// details.
//
// You should have received a copy of the European Union Public Licence (EUPL) v1.2
// along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#ifndef LAMINAR_FLOW_BENCHMARK_H
#    define LAMINAR_FLOW_BENCHMARK_H

/// \author Thomas Gengenbach, Tobias Hahn, Staffan Ronnas and Michael Schick

#    include <cmath>
#    include <utility>
#    include <string>
#    include <fstream>
#    include <vector>

#    include <mpi.h>
#    include "hiflow.h"

using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

typedef LADescriptorCoupledD LAD;

typedef std::vector<double> Coord;

int main ( int argc, char** argv );

/// Application class used to control simulation state and keep track
/// of main variables.

class LaminarFlowApp : public NonlinearProblem<LAD>
{
  public:
    // Application constants. These could perhaps be changed, but do it at your own risk!
    static const int DIM = 3; // problem dimension
    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;

    LaminarFlowApp ( );
    LaminarFlowApp ( PropertyTree &config );
    ~LaminarFlowApp ( );

    void run ( );

    virtual void EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F );
    virtual void EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF );
  private:
    friend class LaminarFlowNL;

    void prepare_space ( );
    void prepare_linear_solver ( );
    void prepare_lin_alg_structures ( );
    void prepare_nls ( );
    void prepare_bc ( );
    void solve_nls ( );
    void visualize_solution ( );

    void compute_residual ( const LAD::VectorType& u, LAD::VectorType* F );
    void compute_matrix ( const LAD::VectorType& u, LAD::MatrixType* F );

    MPI_Comm comm_; // must use concrete object here -> hence C interface
    int num_partitions_;
    int rank_;
    MeshPtr mesh_, master_mesh_;
    VectorSpace<double> space_;
    std::vector<Quadrature<double>*> quadratures_;
    SYSTEM la_sys_;
    Couplings<double> couplings_;
    double rho_, nu_;
    double H_, Um_;

    double time_measurement_res_, time_measurement_mat_;
    int res_count_, mat_count_;

    LAD::MatrixType matrix_;
    LAD::VectorType sol_, rhs_, res_;
    std::vector<int> dirichlet_dofs_;
    std::vector<LAD::DataType> dirichlet_values_;

    // would be better to use pointer to base class here, but not
    // possible since these objects do not have sensible constructors...
    NonlinearSolver<LAD>* nls_;
    DampingStrategy<LAD>* damping_strategy_;
    ForcingStrategy<LAD>* forcing_strategy_;
    bool damping_;
    bool forcing_;
#    ifdef WITH_MUMPS
    ScopedPtr< MumpsSolver<LAD, MumpsStructureD> >::Type linear_solver_;
#    else
    LinearSolver<LAD>* linear_solver_;
#    endif

    int refinement_level_;

};

void convert_to_std_vector ( const LAD::VectorType& in, std::vector<double>& out );

#endif /* _LAMINAR_FLOW_BENCHMARK_H_ */
