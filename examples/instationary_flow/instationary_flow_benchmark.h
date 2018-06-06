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

#ifndef INSTATIONARY_FLOW_BENCHMARK_H
#    define INSTATIONARY_FLOW_BENCHMARK_H

/// \author Thomas Gengenbach, Michael Schick and Staffan Ronnas

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>

#    include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;
using namespace hiflow::la;

static const int DIM = 2;
//static const int DIM = 3;
typedef LADescriptorCoupledD LAD;

typedef std::vector<double> Coord;

int main ( int argc, char** argv );

// InstationaryFlowAssembler ////////////////

class InstationaryFlowAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    InstationaryFlowAssembler ( double nu, double rho )
    : nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void set_newton_solution ( const LAD::VectorType& newton_sol );
    void set_prev_solution ( const LAD::VectorType& prev_sol );

    void set_timestep_parameters ( const std::vector<double>& alphas );

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature );
    void initialize_for_facet ( const Element<double>& element, const Quadrature<double>& quadrature, int facet_number );

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

  private:
    LAD::VectorType const* newton_sol_;
    LAD::VectorType const* prev_sol_;

    double nu_, inv_rho_;
    double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
    FunctionValues<double> prev_ns_vel_[DIM];
    FunctionValues<double> prev_ts_vel_[DIM];
    FunctionValues<double> pressure_k_;
    FunctionValues< Vec<DIM, double> > grad_prev_ns_vel_[DIM];
    FunctionValues< Vec<DIM, double> > grad_prev_ts_vel_[DIM];
};

class InstationaryFlowApp : public NonlinearProblem<LAD>
{
  public:

    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE; //MKL;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;

    InstationaryFlowApp ( );
    ~InstationaryFlowApp ( );

    void run ( );

    virtual void EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F );
    virtual void EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF );

    void compute_residual ( const LAD::VectorType& u, LAD::VectorType* F );
    void compute_matrix ( const LAD::VectorType& u, LAD::MatrixType* F );

  private:

    void read_and_distribute_mesh ( const std::string& filename );

    void prepare_space ( );
    void prepare_linear_solver ( );
    void prepare_lin_alg_structures ( );
    void prepare_bc ( );
    void prepare_nls ( );

    void solve_nlp ( );
    void visualize_solution ( int time_step );

    void compute_alphas ( int sub_step, double delta_t, std::vector<double>* alphas );

    MPI_Comm comm_; // must use concrete object here -> hence C interface
    int rank_;
    int num_partitions_;
    const int master_rank_;

    bool visu_p2s_;

    MeshPtr mesh_;
    VectorSpace<double> space_;

    InstationaryFlowAssembler* local_asm_;

    std::string method_;

    SYSTEM la_sys_;
    Couplings<double> couplings_;
    double rho_, nu_;
    double H_, Um_;

    LAD::MatrixType matrix_;
    LAD::VectorType sol_, sol_prev_, rhs_, res_;
    std::vector<int> dirichlet_dofs_;
    std::vector<LAD::DataType> dirichlet_values_;

    // would be better to use pointer to base class here, but not
    // possible since these objects do not have sensible constructors...
    ScopedPtr< Newton<LAD> >::Type nls_;
#    ifdef WITH_MUMPS
    ScopedPtr< MumpsSolver<LAD, MumpsStructureD> >::Type linear_solver_;
#    else
    ScopedPtr< GMRES<LAD> >::Type linear_solver_;
#    endif

    int refinement_level_;

};

void convert_to_std_vector ( const LAD::VectorType& in, std::vector<double>& out );

#endif /* _INSTATIONARY_FLOW_BENCHMARK_H_ */
