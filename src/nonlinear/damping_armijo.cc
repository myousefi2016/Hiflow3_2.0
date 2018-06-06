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

#include "damping_armijo.h"
#include "newton.h"

namespace hiflow
{

    template<class LAD>
    DampingState ArmijoDamping<LAD>::Init ( ArmijoParam param, int data )
    {
        if ( param == ArmijoMaxLoop ) this->MaxLoop_ = data;
        else return kDampingError;
        return kDampingSuccess;
    }

    template<class LAD>
    DampingState ArmijoDamping<LAD>::Init ( ArmijoParam param, DataType data )
    {
        switch ( param )
        {
            case ArmijoInitial:
                this->Initial_ = data;
                break;

            case ArmijoMinimal:
                this->Minimal_ = data;
                break;

            case ArmijoDecrease:
                this->Decrease_ = data;
                break;

            case ArmijoSuffDec:
                this->SuffDec_ = data;
                break;

            default:
                return kDampingError;
        }

        return kDampingSuccess;
    }

    template<class LAD>
    DampingState ArmijoDamping<LAD>::Update ( const VectorType& cor, const VectorType& rhs, VectorType* res, VectorType* sol, void* myNewton )
    {
        Newton<LAD>* NewtonSolver = ( Newton<LAD>* ) myNewton;
        // trial step

        DataType res_start = res->Norm2 ( );
        LOG_INFO ( "Initial residual", res_start );

        DataType res_cur = res_start;
        this->residual_ = res_start;
        int iter = 0;
        DataType theta = this->Initial_;
        DataType t = this->SuffDec_;
        DataType eta = 0.;
        if ( NewtonSolver->Forcing ( ) )
            NewtonSolver->GetForcingTerm ( eta );

        // New solution and residual

        VectorType sol_backup;
        sol_backup.CloneFrom ( *sol );
        //sol_backup.Update();

        sol->Axpy ( cor, -theta );
        /*sol->Update();
        this->Average(sol);*/
        sol->Update ( );
        NewtonSolver->GetOperator ( )->ApplyFilter ( *sol );
        sol->Update ( );
        if ( NewtonSolver->GetNonConstMode ( ) )
            NewtonSolver->ComputeResidualNonConst ( *sol, rhs, res );
        else
            NewtonSolver->ComputeResidual ( *sol, rhs, res );
        res_cur = res->Norm2 ( );
        LOG_INFO ( "Residual norm (trial step)", res_cur );
        this->residual_ = res_cur;

        // -> Iterate if needed

        const DataType bound = ( 1. - t * ( 1. - eta ) ) * res_start;
        LOG_INFO ( "Armijo damping acceptance bound", bound );

        while ( ( res_cur > bound )
                && ( iter <= this->MaxLoop_ )
                && ( theta > this->Minimal_ ) )
        {
            // restore old solution

            // sol->Axpy(cor,theta);
            sol->CloneFrom ( sol_backup );
            //sol->Update();
            //this->Average(sol);

            // change damping since actual residual is rejected

            ++iter;
            theta *= this->Decrease_;
            eta = 1. - theta * ( 1 - eta );

            // New solution and residual

            sol->Axpy ( cor, -theta );
            /*sol->Update();
            this->Average(sol);*/
            sol->Update ( );
            NewtonSolver->GetOperator ( )->ApplyFilter ( *sol );
            sol->Update ( );
            if ( NewtonSolver->GetNonConstMode ( ) )
                NewtonSolver->ComputeResidualNonConst ( *sol, rhs, res );
            else
                NewtonSolver->ComputeResidual ( *sol, rhs, res );
            res_cur = res->Norm2 ( );

            LOG_INFO ( "Residual norm (damped)", res_cur );
            this->residual_ = res_cur;

            if ( NewtonSolver->Forcing ( ) )
                NewtonSolver->SetForcingTerm ( eta );
        }
        LOG_INFO ( "Damping factor", theta );
        return kDampingSuccess;
    }

    template<class LAD>
    void ArmijoDamping<LAD>::Average ( VectorType* x )
    {
        /*const_iterator first = _ic->begin();
        const_iterator last  = _ic->end();

        while(first != last)
          {
            vector<pair<int,DataType> >::const_iterator first_v;
            vector<pair<int,DataType> >::const_iterator last_v ;

            first_v = ((*first).second).begin();
            last_v  = ((*first).second).end();

            x[(*first).first] = 0.;

            while(first_v != last_v)
              {
                x[(*first).first] += x[first_v->first] * first_v->second;
                ++first_v;
              }

            ++first;

          } // while(first !=...*/
    }

    template<class LAD>
    ArmijoDamping<LAD>::ArmijoDamping ( )
    {
        Initial_ = 1.;
        Minimal_ = 1.e-4;
        Decrease_ = 0.5;
        SuffDec_ = 1.e-4;
        MaxLoop_ = 10;
        name_ = "Armijo";
    }

    template<class LAD>
    ArmijoDamping<LAD>::ArmijoDamping ( DataType init, DataType mini, DataType dec, DataType suffdec, int maxloop )
    : Initial_ ( init ), Minimal_ ( mini ), Decrease_ ( dec ), SuffDec_ ( suffdec ), MaxLoop_ ( maxloop )
    {
        name_ = "Armijo";
    }

    template<class LAD>
    ArmijoDamping<LAD>::~ArmijoDamping ( )
    {
        ForcingTerm_ = NULL;
    }

    /// template instantiation
    template class ArmijoDamping<la::LADescriptorCoupledD>;
    template class ArmijoDamping<la::LADescriptorCoupledS>;
    template class ArmijoDamping<la::LADescriptorHypreD>;
    template class ArmijoDamping<la::LADescriptorPolynomialChaosD>;
    template class ArmijoDamping<la::LADescriptorPolynomialChaosExpansionD>;

} // namespace hiflow
