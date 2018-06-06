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

#ifndef _NONLINEAR_BALANCE_ASSEMBLER_H_
#    define _NONLINEAR_BALANCE_ASSEMBLER_H_

///
/// \file nonlinear_balance_assembler.h
///
/// \author Martin Baumann
///

/// Assembler for system matrix and rhs vector of the nonlinear balance problem
/// that has to be solved for Boussinesq problems to balance the Coriolis
/// force, the pressure gradient force, and the centrifugal force.
/// Given, the initial velocity field, the corresponding initial pressure
/// field is calculated. In a second step, the corresponding initial temperature
/// field can be calculated that prevents winds due to this pressure field.

template<int DIM>
class NonlinearBalanceAssembler : public AssemblyAssistant<DIM>
{
  public:

    // Constructor

    NonlinearBalanceAssembler ( const VectorSpace<double>& NS_space,
                                const LAD::VectorType* NS_sol,
                                double beta,
                                double f_cor
                                )
    : NS_space_ ( NS_space ),
    NS_sol_ ( NS_sol ),
    beta_ ( beta ),
    f_cor_ ( f_cor ),
    kappa_ ( 0.286 ),
    theta_z_ ( 3.0 ),
    theta_0_ ( 300. ),
    R_ ( 287.e-6 ),
    p0_ ( 1.e5 ),
    g_ ( 9.81e-3 )
    {
        ;
    };

    /// assemble the Jacobi matrix        

    void operator() ( const Element& element, // corresponds to temp-space
            const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        this->initialize_for_element ( element, quadrature );

        const int num_q = this->num_quadrature_points ( );

        // loop q
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );
            const int n_dofs = this->num_dofs ( 0 );
            const double z = this->x ( q )[2];

            // 1/rho * {\grad(u) : \grad(v)}
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                    lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 0 ) ) +=
                        wq
                        * 1. / rho ( z )
                    * dot ( this->grad_phi ( j, q, 0 ),
                            this->grad_phi ( i, q, 0 ) )
                    * dJ;
            }
        } //loop q
    } //assemble_local_matrix

    /// assemble the right-hand side vector

    void operator() ( const Element& element,
            const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        this->initialize_for_element ( element, quadrature );

        // Assembler for Navier-Stokes
        AssemblyAssistant<DIM> ns_aa;
        Element ns_el ( NS_space_, element.get_cell_index ( ) );
        ns_aa.initialize_for_element ( ns_el, quadrature );

        // Local data containers for Navier-Stokes velocity
        FunctionValues<double> local_vel [DIM];
        FunctionValues<Vec<DIM> > local_grad_vel[DIM];
        //     FunctionValues<Mat<DIM,DIM> > local_hessian_u;
        //     FunctionValues<Mat<DIM,DIM> > local_hessian_v;
        std::vector<FunctionValues<Mat<DIM, DIM> > > local_hessian;
        local_hessian.resize ( 3 );

        // Prepare local data containers    
        for ( int v = 0; v < DIM; ++v )
        {
            local_vel [v].clear ( );
            local_grad_vel[v].clear ( );

            ns_aa.evaluate_fe_function ( *NS_sol_, v, local_vel [v] );
            ns_aa.evaluate_fe_function_gradients ( *NS_sol_, v, local_grad_vel[v] );

            ns_aa.evaluate_fe_function_hessian ( *NS_sol_, v, local_hessian.at ( v ) );
        }
        //     ns_aa.evaluate_fe_function_hessian(*NS_sol_, 0, local_hessian_u);
        //     ns_aa.evaluate_fe_function_hessian(*NS_sol_, 1, local_hessian_v);

        // Prepare assembly loop
        const int n_q = this->num_quadrature_points ( );
        const int n_dofs = this->num_dofs ( 0 );

        // Loop over quadrature points
        for ( int q = 0; q < n_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );
            const double z = this->x ( q )[2];
            const double laplace_p = p_zz ( z );

            // basic state of pressure contribution
            for ( int i = 0; i < n_dofs; ++i )
            {
                lv[this->dof_index ( i, 0 )] += wq
                        * 1. / rho ( z )
                        * laplace_p
                        * this->phi ( i, q )
                        * dJ;
            }

            // convection contribution (version after differential by parts)
            for ( int v_var = 0; v_var < DIM; ++v_var ) // scalar-product of variational formulation
            {
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int t_var = 0; t_var < DIM; ++t_var ) // sum of laplacian
                        lv[this->dof_index ( i, 0 )] -= wq
                            * local_vel[t_var][q]
                            * local_grad_vel[v_var][q][t_var]
                            * this->grad_phi ( i, q, 0 )[v_var]
                            * dJ;
                }
            }

            // convection contribution (version without differential by parts)
            //       for (int v_var = 0; v_var < DIM; ++v_var)  // sum over divergence
            //       {
            //         for (int t_var = 0; t_var < DIM; ++t_var)  // sum of nonlinear operator
            //         {
            //           for (int i = 0; i < n_dofs; ++i)
            //             lv[this->dof_index(i, 0)] += wq
            //                                        * (     local_grad_vel[v_var][q][t_var]
            //                                              * local_grad_vel[t_var][q][v_var]
            //                                            +   local_vel[t_var][q]
            //                                              * local_hessian.at(v_var)[q](v_var,t_var)
            //                                           )
            //                                        * this->phi(i,q, 0)
            //                                        * dJ * (1);
            //         }
            //       }

            //       // laplace contribution
            //       const double nu = 5.e-3;
            //       for (int i = 0; i < n_dofs; ++i)
            //         lv[this->dof_index(i, 0)] -= wq
            //                                    * nu
            //                                    * (
            //                                        (   local_hessian_u[q](0,0)
            //                                          + local_hessian_u[q](1,1)
            //                                          + local_hessian_u[q](2,2)
            //                                        ) * this->grad_phi(i, q, 0)[0]
            //                                      + (   local_hessian_v[q](0,0)
            //                                          + local_hessian_v[q](1,1)
            //                                          + local_hessian_v[q](2,2)
            //                                        ) * this->grad_phi(i, q, 0)[1]
            //                                      ) * dJ;

            //       // coriolis contribution
            //       interminable_assert(fabs(Coriolis(this->x(q)[1]))<1.0e-14);
            //       for (int i = 0; i < n_dofs; ++i)
            //       {
            //         lv[this->dof_index(i, 0)] -= wq * ( - Coriolis(this->x(q)[1])
            //                                             * local_vel[1][q]
            //                                             * this->grad_phi(i, q, 0)[0]
            //                                             
            //                                             + Coriolis(this->x(q)[1]) 
            //                                             * local_vel[0][q]
            //                                             * this->grad_phi(i, q, 0)[1]
            //                                             ) * dJ;
            //       }

        }//loop over q
    }//assemble_local_vector

  private:

    // given Navier-Stokes solution vector (global vector) and space
    const LAD::VectorType* NS_sol_;
    const VectorSpace<double>& NS_space_;

    double Coriolis ( double y ) const
    {
        return this->f_cor_ + this->beta_*y;
    }

    double rho ( double z ) const
    {
        return p0_ * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, 0.1e1 / kappa_ - 0.1e1 ) / R_ / ( double ) ( theta_0_ + theta_z_ * z );
    }

    double p_zz ( double z ) const
    {
        return p0_ * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, 0.1e1 / kappa_ ) * ( double ) ( int ) pow ( ( double ) theta_0_, ( double ) ( -2 ) ) * ( double ) ( int ) pow ( ( double ) ( 1 + theta_z_ * z / theta_0_ ), ( double ) ( -2 ) ) * g_ * g_ * pow ( R_, -0.2e1 ) * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, -0.2e1 ) + p0_ * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, 0.1e1 / kappa_ ) * ( double ) ( int ) pow ( ( double ) theta_0_, ( double ) ( -2 ) ) * ( double ) ( int ) pow ( ( double ) ( 1 + theta_z_ * z / theta_0_ ), ( double ) ( -2 ) ) * g_ / R_ / ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1 ) * ( double ) theta_z_ - p0_ * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, 0.1e1 / kappa_ ) * ( double ) ( int ) pow ( ( double ) theta_0_, ( double ) ( -2 ) ) * ( double ) ( int ) pow ( ( double ) ( 1 + theta_z_ * z / theta_0_ ), ( double ) ( -2 ) ) * g_ * g_ * pow ( R_, -0.2e1 ) * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, -0.2e1 ) * kappa_;
    }

    // parameters
    double beta_, f_cor_, kappa_, theta_z_, theta_0_, R_, p0_, g_;
};

/// Assembler for system matrix and rhs vector of the stationary temperature
/// problem to calculate initial temperature field for usage in Boussinesq 
/// applications.

template<int DIM>
class StationaryHeatAssembler : public AssemblyAssistant<DIM>
{
  public:

    // Constructor

    StationaryHeatAssembler ( const VectorSpace<double>& NS_space,
                              const LAD::VectorType* NS_sol,
                              double thermal_diffusivity
                              )
    : NS_space_ ( NS_space ),
    NS_sol_ ( NS_sol ),
    thermal_diffusivity_ ( thermal_diffusivity )
    {
        ;
    };

    /// assemble the Jacobi matrix        

    void operator() ( const Element& element, // corresponds to temp-space
            const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        this->initialize_for_element ( element, quadrature );

        // Assembler for Navier-Stokes
        AssemblyAssistant<DIM> ns_aa;
        Element ns_el ( NS_space_, element.get_cell_index ( ) );
        ns_aa.initialize_for_element ( ns_el, quadrature );

        // Local data containers for Navier-Stokes velocity
        FunctionValues<double> local_vel[DIM];

        // Prepare local data containers    
        for ( int v = 0; v < DIM; ++v )
        {
            local_vel [v].clear ( );
            ns_aa.evaluate_fe_function ( *NS_sol_, v, local_vel [v] );
        }

        // loop q
        for ( int q = 0; q<this->num_quadrature_points ( ); ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );
            const int n_dofs = this->num_dofs ( 0 );

            // a * (\grad(vel), \grad(phi)) + (vel*nabla(theta), phi)
            for ( int i = 0; i < n_dofs; ++i ) // -> test
            {
                for ( int j = 0; j < n_dofs; ++j ) // -> trial
                {
                    lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 0 ) ) +=
                            wq * ( local_vel[0][q] * this->grad_phi ( j, q )[0]
                            + local_vel[1][q] * this->grad_phi ( j, q )[1]
                            + local_vel[2][q] * this->grad_phi ( j, q )[2] )
                            * this->phi ( i, q ) * dJ;
                    lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 0 ) ) +=
                            wq * thermal_diffusivity_
                            * dot ( this->grad_phi ( j, q, 0 ), this->grad_phi ( i, q, 0 ) )
                            * dJ;
                }
            }
        } //loop q
    } //assemble_local_matrix

    /// assemble the right-hand side vector

    void operator() ( const Element& element,
            const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        // - vel_3 * d/dz Theta_0   -> should be zero since vel_3 = 0

        /*
        this->initialize_for_element(element, quadrature);

        // Assembler for Navier-Stokes
        AssemblyAssistant<DIM> ns_aa;
        Element ns_el(NS_space_, element.get_cell_index ());
        ns_aa.initialize_for_element(ns_el, quadrature);

        // Local data containers for Navier-Stokes velocity
        FunctionValues<double>        local_vel     [DIM];

        // Prepare local data containers    
        for (int v=0; v<DIM; ++v)
        {
          local_vel     [v].clear();
          ns_aa.evaluate_fe_function          (*NS_sol_, v, local_vel     [v]);
        }

        // Prepare assembly loop
        const int n_q    = this->num_quadrature_points();
        const int n_dofs = this->num_dofs(0);

        // Loop over quadrature points
        for (int q = 0; q < n_q; ++q)
        {
          const double wq = this->w(q);
          const double dJ = std::abs(this->detJ(q));

    //       // u-convection contribution
    //       for (int i = 0; i < n_dofs; ++i)
    //       {
    //         lv[this->dof_index(i, 0)] -= wq * (   this->x(q)[2] //local_vel[0][q]
    //                                             * this->phi(i, q)
    //                                           ) * dJ;
    //       }

          // convection contribution
          for (int v_var = 0; v_var < DIM; ++v_var)  // scalar-product of variational formulation
          {
            for (int i = 0; i < n_dofs; ++i)
            {
              for (int t_var = 0; t_var < DIM; ++t_var)  // scalar-product of nonlinearity
                lv[this->dof_index(i, 0)] -= wq * (   local_vel[t_var][q]
         * local_grad_vel[v_var][q][t_var]
         * this->grad_phi(i, q, 0)[v_var] 
                                                  ) * dJ;
            }
          }

          // coriolis contribution
          for (int i = 0; i < n_dofs; ++i)
          {
            lv[this->dof_index(i, 0)] -= wq * ( - Coriolis(this->x(q)[1])
         * local_vel[1][q]
         * this->grad_phi(i, q, 0)[0]

                                                + Coriolis(this->x(q)[1]) 
         * local_vel[0][q]
         * this->grad_phi(i, q, 0)[1]
                                                ) * dJ;
          }

        }//loop over q

         */
    }//assemble_local_vector

  private:

    // given Navier-Stokes solution vector (global vector) and space
    const LAD::VectorType* NS_sol_;
    const VectorSpace<double>& NS_space_;

    // parameters
    double thermal_diffusivity_;
};

/// Assembler for mass matrix and rhs vector of the L2-projection of the
/// temperature field, corresponding to the balance pressure calculated
/// using the NonlinearBalanceAssembler.

template<int DIM>
class BalancedTemperatureAssembler : public AssemblyAssistant<DIM>
{
  public:

    // Constructor

    BalancedTemperatureAssembler ( const LAD::VectorType* pressure,
                                   double theta_0,
                                   double theta_z,
                                   double gravity
                                   )
    : pressure_ ( pressure ),
    theta_0_ ( theta_0 ),
    theta_z_ ( theta_z ),
    kappa_ ( 0.286 ),
    R_ ( 287.e-6 ),
    p0_ ( 1.e5 ),
    g_ ( gravity )
    {
        ;
    };

    /// assemble the Jacobi matrix        

    void operator() ( const Element& element, // corresponds to temp-space
            const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        this->initialize_for_element ( element, quadrature );

        const int num_q = this->num_quadrature_points ( );

        // loop q
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );
            const int n_dofs = this->num_dofs ( 0 );

            // scalar product (u,v)
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                    lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 0 ) ) += wq * this->phi ( j, q, 0 ) * this->phi ( i, q, 0 ) * dJ;
            }
        } //loop q
    } //assemble_local_matrix

    /// assemble the right-hand side vector

    void operator() ( const Element& element,
            const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        this->initialize_for_element ( element, quadrature );

        // Prepare gradient for pressure
        FunctionValues<Vec<DIM> > local_grad;
        local_grad.clear ( );
        this->evaluate_fe_function_gradients ( *pressure_, 0, local_grad );

        // Prepare assembly loop
        const int n_q = this->num_quadrature_points ( );
        const int n_dofs = this->num_dofs ( 0 );

        // Loop over quadrature points
        for ( int q = 0; q < n_q; ++q )
        {
            const double wq = this->w ( q );
            const double dJ = std::abs ( this->detJ ( q ) );

            const double z = this->x ( q )[2];
            const double theta = theta_0_ + z*theta_z_;

            // Contribution
            for ( int i = 0; i < n_dofs; ++i )
                lv[this->dof_index ( i, 0 )] += wq * ( theta * 1. / rho ( z ) * 1. / g_ )
                * ( local_grad[q][2] * this->phi ( i, q, 0 )
                    //+ p_z(z) * this->phi(i, q, 0)
                    ) * 1 * dJ;
        }//loop over q
    }//assemble_local_vector

  private:

    // given Navier-Stokes solution vector (global vector) and space
    const LAD::VectorType* pressure_;

    double rho ( double z ) const
    {
        return p0_ * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, 0.1e1 / kappa_ - 0.1e1 ) / R_ / ( double ) ( theta_0_ + theta_z_ * z );
    }

    double p_z ( double z ) const
    {
        return -p0_ * pow ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1, 0.1e1 / kappa_ ) / ( double ) theta_0_ / ( double ) ( 1 + theta_z_ * z / theta_0_ ) * g_ / R_ / ( -log ( ( double ) ( 1 + theta_z_ * z / theta_0_ ) ) * g_ * kappa_ / R_ / ( double ) theta_z_ + 0.1e1 );
    }

    // parameters
    double theta_0_, theta_z_, g_, kappa_, R_, p0_;
};

#endif
