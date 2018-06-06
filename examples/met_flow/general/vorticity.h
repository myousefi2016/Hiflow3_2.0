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

#ifndef Vorticity_H
#    define Vorticity_H

///
/// \file vorticity.h
/// \brief AssemblyAssistant for vorticity and derived quantities.
/// \details Containing assistants for the assembly of the vorticity derived from a given velocity field, 
/// the assembly of the velocity derived from given stream function and the assembly of the stream function 
/// derived from given vorticity.
/// \author Volker Lange, Teresa Beck
///
/// TODO inlcude HYPRE laternative for ilupp

#    include <cmath>
#    include <fstream>
#    include <iostream>
#    include <string>

#    include "hiflow.h"
#    include "../../src/common/log.h"
#    include "assembly/assembly_assistant.h"
#    include "../tmp_config/met_flow_vars.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;
using namespace hiflow::la;

typedef la::SeqDenseMatrix<double> LocalMatrix;
typedef std::vector<double> LocalVector;

/// Assembler used for velocity -> vorticity

template<int DIM>
class VortAssembler : private AssemblyAssistant<DIM, double>
{
  public:
    /// If true the velocity is given by an analytical function, else a vector containing the velocity must be provided.

    void set_from_function ( bool frm_fnct )
    {
        from_function = frm_fnct;
    }

    /// Set the function describing the velocity. 

    void set_function ( double (*ptr2func )( Vec<DIM, double> ) )
    {
        velo_f = ptr2func;
    }

    /// Set the vectors containing the velocity.

    void set_velocity_solution ( const LAD::VectorType* velocity_u, const LAD::VectorType* velocity_v )
    {
        velocity_u_ = velocity_u;
        velocity_v_ = velocity_v;
    }
    /// Assembles mass-matrix    

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        // compute local matrix
        const int num_q = this->num_quadrature_points ( );
        const int num_d = this->num_dofs ( 0 );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );

            for ( int i = 0; i < num_d; ++i )
            {
                for ( int j = 0; j < num_d; ++j )
                {
                    //A=int(phi_i dot phi_j)
                    lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 0 ) ) +=
                            wq * this->phi ( i, q, 0 ) * this->phi ( j, q, 0 ) * dJ;
                } // j
            } // i
        }// q
    }// assemble_local_matrix

    /// Assembles residual vector.

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        if ( !from_function )
        {
            grad_[0].clear ( );
            grad_[1].clear ( );
            this->evaluate_fe_function_gradients ( *velocity_u_, 0, grad_[0] );
            this->evaluate_fe_function_gradients ( *velocity_v_, 0, grad_[1] );
        }

        const int num_q = this->num_quadrature_points ( );
        double value;

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = this->detJ ( q );

            if ( from_function )
            {
                value = ( *velo_f )( this->x ( q ) );
            }
            else
            {
                value = grad_[1][q][0] - grad_[0][q][1];
            }

            for ( int i = 0; i < this->num_dofs ( 0 ); ++i )
            {
                lv[this->dof_index ( i, 0 )] += wq * value * this->phi ( i, q ) * dJ;
            } // i
        } // q
    }//assemble_local_vector

  private:
    const LAD::VectorType* velocity_u_;
    const LAD::VectorType* velocity_v_;
    bool from_function;
    double (*velo_f )( Vec<DIM, double> );
    FunctionValues< Vec<DIM, double> > grad_[DIM];
};

// Assemblers used for: vorticity -> velocity

/// Assembler to calculate streamfunction from vorticity.

template<int DIM>
class StreamAssembler : private AssemblyAssistant<DIM, double>
{
  public:
    /// If true the vorticity is given by an analytical function, else a vector containing the vorticity must be provided.

    void set_from_function ( bool frm_fnct )
    {
        from_function = frm_fnct;
    }

    /// Set the function describing the vorticity. 

    void set_function ( double (*ptr2func )( Vec<DIM, double> ) )
    {
        vort_f = ptr2func;
    }

    /// Set the vector containing the vorticity.

    void set_vorticity_solution ( const LAD::VectorType* vorticity )
    {
        vorticity_ = vorticity;
    }

    /// Assembles Jacobi-Matrix    

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        // compute local matrix
        const int num_q = this->num_quadrature_points ( );

        //add number of DoFs of this element
        const int num_local_dofs = this->num_dofs ( 0 );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );

            for ( int i = 0; i < num_local_dofs; ++i )
            {
                for ( int j = 0; j < num_local_dofs; ++j )
                {//A=(grad_phi_i dot grad_phi_j)
                    lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 0 ) ) +=
                            wq * dot ( this->grad_phi ( j, q ), this->grad_phi ( i, q ) ) * dJ;
                }
            }
        }
    }//assemble_local_matrix

    /// Assembles residual vector.  

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        if ( !from_function )
        {
            vort_.clear ( );
            this->evaluate_fe_function ( *vorticity_, 0, vort_ );
        }

        const int num_q = this->num_quadrature_points ( );
        //number of DoFs of this element
        const int num_local_dofs = this->num_dofs ( 0 );
        double vort;

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = this->detJ ( q );

            if ( from_function )
            {
                vort = ( *vort_f )( this->x ( q ) );
            }
            else
            {
                vort = vort_[q];
            }

            for ( int i = 0; i < num_local_dofs; ++i )
            {
                //b=(xi dot phi)
                lv[this->dof_index ( i, 0 )] += wq * vort * this->phi ( i, q ) * dJ;
            }
        }
    }//assemble_local_vector

  private:
    const LAD::VectorType* vorticity_;
    bool from_function;
    double (*vort_f )( Vec<DIM, double> );
    FunctionValues<double> vort_;
};

/// Assembler to calculate velocity from stream function

template<int DIM>
class VeloAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    /// set the stream function solution

    void set_stream_solution ( const LAD::VectorType* stream_sol )
    {
        sol_str_ = stream_sol;
    }

    /// set the previous newton solution / initial guess

    void set_newton_solution ( const LAD::VectorType* newton_sol )
    {
        prev_newton_sol_ = newton_sol;
    }

    /// assemble the Jacobi matrix  

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {

        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        const int num_q = this->num_quadrature_points ( );

        // loop q
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    {
                        //[BUG] seg fault is caused in phi(i,q,u_var)
                        lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) +=
                                wq * this->phi ( j, q, u_var ) * this->phi ( i, q, u_var ) * dJ;
                    }
                }
            }

            const int p_var = DIM;
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                // assemble b = - inv_rho: * (grad_phi_i dot phi_j) for p   
                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < this->num_dofs ( p_var ); ++j )
                    {
                        lm ( this->dof_index ( i, u_var ), this->dof_index ( j, p_var ) ) +=
                                wq * (
                                -this->grad_phi ( i, q, u_var )[u_var] * this->phi ( j, q, p_var )
                                ) * dJ;
                    }
                }
            }
            //*** incompress. terms ***************************
            const int q_var = DIM;
            // assemble bT = - inv_rho: * (phi_i dot grad_phi_j)
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i < this->num_dofs ( q_var ); ++i )
                {
                    for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                    {
                        lm ( this->dof_index ( i, q_var ), this->dof_index ( j, u_var ) ) +=
                                wq * (
                                -this->phi ( i, q, q_var ) * this->grad_phi ( j, q, u_var )[u_var]
                                ) * dJ;
                    }
                }
            }

        }//loop q
    }//assemble_local_matrix

    /// assemble the residual vector  

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            Vector<double>& lv )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

        // recompute previous velocity solution values
        for ( int d = 0; d < DIM; ++d )
        {
            vel_[d].clear ( );
            grad_vel_[d].clear ( );

            this->evaluate_fe_function ( *prev_newton_sol_, d, vel_[d] );
            this->evaluate_fe_function_gradients ( *prev_newton_sol_, d, grad_vel_[d] );
        }//d

        // recompute pressure solution values
        p_.clear ( );
        this->evaluate_fe_function ( *prev_newton_sol_, DIM, p_ );

        // recompute streamfunction solution values
        grad_str_.clear ( );
        this->evaluate_fe_function_gradients ( sol_str_, 0, grad_str_ );

        const int num_q = this->num_quadrature_points ( );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );

            // l(u) =-(d_y xi_x * phi)
            for ( int i = 0; i < this->num_dofs ( 0 ); ++i )
            {
                lv[this->dof_index ( i, 0 )] += wq * (
                        ( vel_[0][q] + grad_str_[q][1] ) * this->phi ( i, q, 0 )
                        - p_[q] * this->grad_phi ( i, q, 0 )[0]
                        ) * dJ;
            }

            // l(v) =(d_x xi_y * phi)  
            for ( int i = 0; i < this->num_dofs ( 1 ); ++i )
            {
                lv[this->dof_index ( i, 1 )] += wq * (
                        ( vel_[1][q] - grad_str_[q][0] ) * this->phi ( i, q, 1 )
                        - p_[q] * this->grad_phi ( i, q, 1 )[1]
                        ) * dJ;
            }

            // l(q)
            const int q_var = DIM;
            double div_u_k = 0.;
            for ( int d = 0; d < DIM; ++d )
            {
                div_u_k += grad_vel_[d][q][d];
            }

            for ( int i = 0; i < this->num_dofs ( q_var ); ++i )
            {
                lv[this->dof_index ( i, q_var )] += wq * (
                        -div_u_k * this->phi ( i, q, q_var )
                        ) * dJ;
            }

        }//loop over q

    }//assemble_local_vector

  private:
    const LAD::VectorType* sol_str_;
    FunctionValues< Vec<DIM, double> > grad_str_;

    const LAD::VectorType* prev_newton_sol_;
    FunctionValues<double> vel_[DIM];
    // pressure at previous newton step
    FunctionValues<double> p_;
    // gradient of velocity at previous newton step
    FunctionValues< Vec<DIM, double> > grad_vel_[DIM];
};

/// Vorticity base class

template<int DIM>
class Vorticity
{
  public:

    Vorticity ( );
    Vorticity ( mesh::MeshPtr mesh, int degrees_vars );
    ~Vorticity ( );

    void initialize ( );

    inline bool is_initialized ( )
    {
        return initialized;
    }
    void setup ( mesh::MeshPtr mesh, int velo_degree );

    inline bool is_setup ( )
    {
        return was_setup;
    }
    void periodify_space ( );
    void free_mem ( );

    /// Calculate vorticity from a FEM solution vector containing the velocity
    void calc_vort ( LAD::VectorType& velo, VectorSpace<double> const& space );
    /// Calculate vorticity from an analytical function (using a function pointer)
    void calc_vort ( double (*ptr2func )( Vec<DIM, double> ) );

    void calc_velo ( double (*f_ )( Vec<DIM, double> ) );

    void calc_stream ( LAD::VectorType& vort ); // calculate streamfunction from a FEM solution vector containing the vorticity
    void calc_stream ( double (*ptr2func )( Vec<DIM, double> ) ); // calculate streamfunction from an analytical function (using a functionpointer)

    inline bool vorticity_solved ( )
    {
        return vort_solved;
    }

    inline bool stream_solved ( )
    {
        return streamsolved;
    }

    LAD::VectorType* get_vort ( )
    {
        /*assert(vort_solved);*/ return &sol_vort_;
    }

    LAD::VectorType* get_stream ( )
    {
        assert ( streamsolved );
        return &sol_stream_;
    }

    /// VectorSpace owned by Vorticity, related to scalar problems (u,v,vort).

    VectorSpace<double> const& space ( )
    {
        return *space_;
    }

    /// Couplings related to the VectorSpace owned by Vorticity

    Couplings<double>* couplings ( )
    {
        return couplings_;
    }

    MeshPtr* mesh ( )
    {
        return &mesh_;
    }

  private:

    const MPI_Comm& communicator ( )
    {
        return comm_;
    }

    int rank ( )
    {
        return rank_;
    }

    int num_partitions ( )
    {
        return num_partitions_;
    }

    PLATFORM la_platform ( ) const
    {
        return la_sys_.Platform;
    }

    IMPLEMENTATION la_implementation ( ) const
    {
        return la_impl_;
    }

    MATRIX_FORMAT la_matrix_format ( ) const
    {
        return la_matrix_format_;
    }

#    ifdef USE_HYPRE  
    MPI_Comm comm_;
#    else
    const MPI_Comm comm_;
#    endif
    int rank_;
    int num_partitions_;

    bool initialized, was_setup, vort_solved, streamsolved;

    SYSTEM la_sys_;
    IMPLEMENTATION la_impl_;
    MATRIX_FORMAT la_matrix_format_;
    Couplings<double> * couplings_;
    SparsityStructure * sparsity_;

    MeshPtr mesh_;
    VectorSpace<double> * space_;
    const std::vector<int> degrees;

    // Linear solver
    GMRES<LAD> solver_;
#    ifndef USE_HYPRE  
    PreconditionerIlupp<LAD> ilupp_;
#    endif  
    // Assemblers
    //  StandardGlobalAssembler * global_asm_;
    HpFemAssembler<double> * global_asm_;
    VortAssembler<DIM> * vort_asm_;
    StreamAssembler<DIM> * stream_asm_;

    // Matrizes and vectors
    LAD::MatrixType matrix_vort_, matrix_stream_;
    LAD::VectorType sol_, sol_vort_, sol_stream_, rhs_vort_, rhs_stream_, cor_;

    // Member function
    void setup_linear_algebra ( const std::string platform_str, const std::string impl_str, const std::string matrix_str );
    void setup_linear_algebra ( ); // using default values

    void setup_solver ( const int max_iter, const double abs_tol, const double rel_tol, const double div_tol );
    void setup_solver ( ); // using default values
};

#endif
