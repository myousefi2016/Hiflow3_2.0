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

/// \author Michael Schick, Chen Song

#ifndef HIFLOW_EXAMPLES_POISSON_UNCERTAINTY_H_
#    define HIFLOW_EXAMPLES_POISSON_UNCERTAINTY_H_

#    include <mpi.h>
#    include "hiflow.h"

using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::polynomialchaos;

int main ( int argc, char** argv );

typedef LADescriptorCoupledD LAD;
/// Vector for coordinates in spatial variable
typedef std::vector<double> Coord;
/// Deterministic space dimension
static const int DIM = 2;

class PoissonMPIModeAssembler : private AssemblyAssistant<DIM, double>
{
  public:
    typedef la::LADescriptorCoupledD LAD;

    /// Constructor

    PoissonMPIModeAssembler ( int mode, double nu )
    : mode_ ( mode ), nu_ ( nu )
    {
    }

    /// Compute local finite-element matrix
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
    /// Compute local finite-element residual
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

  protected:
    /// Index of Polynomial Chaos mode
    int mode_;
    /// Diffusivity parameter
    double nu_;
};

class PoissonMPI
{
  public:
    /// Constructor
    PoissonMPI ( );
    /// Constructor using property tree
    explicit PoissonMPI ( PropertyTree &config, const MPI_Comm& comm_space, const MPI_Comm& comm_uq );
    /// Destructor
    virtual ~PoissonMPI ( );
    /// Run application
    void run ( );

  protected:
    /// Setup Finite Element mesh
    void setup_mesh ( );
    /// Setup application specific parameters
    void setup_application ( );
    /// Setup Finite-Element space
    void setup_space ( );
    /// Setup CPU system
    void setup_system ( );
    /// Setup Linear algebra structures
    void setup_la_structures ( );
    /// Setup linear solver
    void setup_linear_solver ( );
    /// Setup boundary conditions
    void setup_bc ( );
    /// Solve the discretized linear system
    void solve_linear_system ( );

    /// Assemble the stiffness Galerkin matrix
    void compute_matrix ( );
    /// Compute the residual
    void compute_residual ( );

    /// Visualize the Polynomial Chaos modes
    void visualize_solution ( LAD::VectorType& u, std::string const& filename ) const;

    /// MPI Communicator for spatial variable
    const MPI_Comm comm_space_;
    /// MPI Communicator for Polynomial Chaos modes
    const MPI_Comm comm_uq_;
    /// Configuration data provided in a separate xml file
    PropertyTree* config_;

    /// Number of mesh partitions
    int num_partitions_;
    /// Refinement level of mesh
    int refinement_level_;
    /// MPI rank of spatial communicator
    int rank_;
    /// Finite element mesh
    MeshPtr mesh_, master_mesh_;

    /// Finite Element dofs for Dirichlet boundary conditions
    std::vector<int> dirichlet_dofs_;
    /// Values of Dirichlet boundary conditions
    std::vector<LAD::DataType> dirichlet_values_;

    /// System variable
    SYSTEM la_sys_;
    /// Implementation types (CPU)
    IMPLEMENTATION matrix_impl_, vector_impl_;
    /// Implementation types (CPU)
    MATRIX_FORMAT la_matrix_format_;
    /// Implementation types (CPU)
    MATRIX_FREE_PRECOND matrix_precond_;

    /// Stiffness matrix of dimensions of spatial variable
    LAD::MatrixType matrix_;
    /// Allocated modes of stiffness Galerkin matrix
    std::vector<LAD::MatrixType*> galerkin_matrix_;
    /// Couplings of Finite-Element mesh and ansatz functions
    Couplings<double> couplings_;

    /// Dummy solution of dimension of spatial variable
    LAD::VectorType sol_;
    /// Temporary Galerkin modes allocated in application
    std::vector<CoupledVector<LAD::DataType>* > galerkin_sol_tmp_;

    /// Abstract linear solver interface
    la::LinearSolver<LADescriptorPolynomialChaosD>* linear_solver_;

    /// CG solver
    la::CG<LADescriptorPolynomialChaosD> cg_solver_;

#    ifdef WITH_UMFPACK
    /// Mean based preconditioning using UMFPACK
    polynomialchaos::MeanbasedPreconditioner<LADescriptorPolynomialChaosD>* mean_precond_;
#    endif
    /// Multilevel linear solver
    polynomialchaos::PCMultilevelSolver<LADescriptorPolynomialChaosD>* ml_precond_;

    /// Finite-Element space
    VectorSpace<double>* space_;

    /// Parameters for Polynomial Chaos expansion
    int No_, N_, No_input_;
    /// Basis class for Polynomial Chaos
    PCBasis pcbasis_;
    /// Third order Galerkin tensor
    PCTensor pctensor_;
    /// Galerkin vectors for solution, residuals and corrections
    PCGalerkinVector<LAD::DataType> *galerkin_sol_, *galerkin_res_, *galerkin_cor_;
    /// Galerkin stiffness matrix
    PCGalerkinMatrix<LAD::DataType>* system_matrix_;

    /// Assemler for stiffness matrix
    PoissonMPIModeAssembler* local_mode_asm_;
    /// Uncertainty coefficients for input parameter
    std::vector<double> nu_;
};

#endif
