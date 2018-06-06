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

/// \author Thomas Gengenbach
/// \author Eva Ketelaer

#ifndef _CONVDIFF_NEU_H_
#    define _CONVDIFF_NEU_H_

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>

#    include "hiflow.h"
#    include "config.h"
#    include "common/log.h"

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

static const int DIM = 2;

typedef LADescriptorCoupledD LAD;
typedef LAD::MatrixType CMatrix;
typedef LAD::VectorType CVector;
typedef LAD::DataType Scalar;

typedef std::vector<double> Coord;

int main ( int argc, char** argv );

// Exception types

struct UnexpectedParameterValue : public std::runtime_error
{

    UnexpectedParameterValue ( const std::string& name, const std::string& value )
    : std::runtime_error (
                           "Unexpected value '" + value + "' for parameter " + name )
    {
    }
};

struct TimingData
{
    double time_elapsed;
};

class TimingScope
{
  public:

    explicit TimingScope ( const std::string& name )
    {
        if ( report_ )
        {
            report_->begin_section ( name );
        }
    }

    explicit TimingScope ( int iteration )
    {
        if ( report_ )
        {
            std::stringstream sstr;
            sstr << "Iteration " << iteration;
            report_->begin_section ( sstr.str ( ) );
            timer_.reset ( );
            timer_.start ( );
        }
    }

    ~TimingScope ( )
    {
        timer_.stop ( );
        if ( report_ )
        {
            TimingData* data = report_->end_section ( );
            data->time_elapsed = timer_.get_duration ( );
        }
    }

    static void set_report ( HierarchicalReport<TimingData>* report )
    {
        report_ = report;
    }

  private:
    static HierarchicalReport<TimingData>* report_;
    Timer timer_;
};

HierarchicalReport<TimingData>* TimingScope::report_ = 0;

class TimingReportOutputVisitor
{
  public:

    explicit TimingReportOutputVisitor ( std::ostream& os )
    : os_ ( os ), level_ ( 0 )
    {
    }

    void enter ( const std::string& name, TimingData* data )
    {
        if ( name == "root" )
        {
            os_ << "+++ Timing Report +++\n\n";
        }
        else
        {
            for ( int l = 0; l < level_; ++l )
            {
                os_ << "  ";
            }
            os_ << name << " took " << data->time_elapsed << " s.\n";
            ++level_;
        }
    }

    void exit ( const std::string& name, TimingData* data )
    {
        if ( name == "root" )
        {
            os_ << "\n+++ End Timing Report +++\n\n";
        }
        else
        {
            --level_;
        }
    }

  private:
    std::ostream& os_;
    int level_;
};

// Homogenous Dirichlet boundary conditions

struct Dirichlet_0_BC
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        // homogenous Dirichlet boundary
        return std::vector<double>( coords_on_face.size ( ), 0. );
    }
};

// ConvectionDiffusionAssembler /////////////

class ConvDiffAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    ConvDiffAssembler ( const std::vector<double>& beta, double nu )
    : beta_ ( beta ), nu_ ( nu )
    {
    }

    ConvDiffAssembler ( )
    : beta_ ( 0 ), nu_ ( 0.0 )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv );

    const std::vector<double> beta_;
    const double nu_;
};

// Convection Diffusion Simulation

class ConvDiff
{
  public:
    explicit ConvDiff ( const std::string& param_filename );
    ~ConvDiff ( );

    void run ( );

  private:
    void read_and_distribute_mesh ( );
    void setup_system ( );
    void initialize_platform ( );
    void prepare_space ( );
    void prepare_bc ( );
    void prepare_linear_solver ( );
    void set_up_preconditioner ( );
    void prepare_lin_alg_structures ( );
    void assemble_system ( );
    void solve_system ( );
    void visualize_solution ( );

    PropertyTree param_;

    // variables for MPI
    MPI_Comm comm_;
    int rank_;
    int num_partitions_;
    const int master_rank_;

    // variables for mesh
    MeshPtr mesh_;
    int refinement_level_;

    // Linear algebra stuff
    bool init_platform_;
    SYSTEM la_sys_;
    IMPLEMENTATION matrix_impl_, vector_impl_;
    MATRIX_FORMAT la_matrix_format_;
    MATRIX_FREE_PRECOND matrix_precond_;

    // variables for assembling and FEM space
    VectorSpace<double> space_;
    ConvDiffAssembler* local_asm_;

    // names
    std::string method_;

    Couplings<double> couplings_;
    double nu_;
    std::vector<double> beta_;

    CMatrix matrix_, *dev_matrix_;
    CVector sol_, rhs_, *dev_sol_, *dev_rhs_;

    // variables for dirichlet boundary
    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;

    LinearSolver<LAD>* linear_solver_;
    Preconditioner<LAD> *precond_;
    // would be better to use pointer to base class here, but not
    // possible since these objects do not have sensible constructors...
    // GMRES<LAD> solver_;
    // ilupp preconditioner
#    ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
#    endif

    bool use_ilupp_;
    HierarchicalReport<TimingData> time_report_;
};

#endif /* _CONVDIFF_NEU_H_ */
