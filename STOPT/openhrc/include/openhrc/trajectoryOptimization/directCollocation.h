/*
 *
 * Copyright (C) 2018
 * Gustavo Arechavaleta <garechav@cinvestav.edu.mx>, Alvaro Paz
 * CINVESTAV - Saltillo Campus
 *
 * This file is part of OpenHRC
 * OpenHRC is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * OpenHRC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

/**
 *	\file include/openhrc/trajectoryOptimization/directCollocation.h
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2018
 *
 *	Class to implement the Trajectory Optimization -> Direct Collocation
 */

#ifndef HR_TRAJECTORY_OPTIMIZATION_DIRECT_COLLOCATION_H
#define HR_TRAJECTORY_OPTIMIZATION_DIRECT_COLLOCATION_H

#include <memory>
#include "openhrc/types.h"
#include "openhrc/core.h"
#include "openhrc/dynamics.h"

namespace hr{
namespace core{


//! The type of the differentiation
/*! \ingroup core_module  */
enum DifferentiationType
{
    wrt_time,            /*!< Differentiate wrt the time. */
    wrt_state,           /*!< Differentiate wrt the robot state {q,dq}. */
    wrt_controlPoints    /*!< Differentiate wrt the control points. */
};


struct robotSettingsTrajectoryOptimization {
    //! General Parameters
    short int n;                            // degrees of freedom
    short int numberControlPoints;          // number of control points
    int numberPartitions;                   // numbver of partitions
    real_t si, sf;                          // initial and final time
    VectorXr S;                             // vector time parameter
    int numberConstraints;                  // number of constraints

    //! Differentiate with respect to
    DifferentiationType DifferentiationWRT;

    //! Boundaries
    VectorXr initialConfiguration;
    VectorXr finalConfiguration;
    VectorXr initialGeneralizedVelocity;
    VectorXr finalGeneralizedVelocity;

    //! consider include constraints customization
};


class DirectCollocation
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor
    DirectCollocation(std::shared_ptr< MultiBody > robot, robotSettingsTrajectoryOptimization robotSettings) : robotMotion{ robot } {

        this->robot = robot;
        this->n = robot->getDoF();

        this->robotSettings = robotSettings;
        this->S = robotSettings.S;

        this->f = 0;


        this->q_initial  = robotSettings.initialConfiguration;
        this->q_final    = robotSettings.finalConfiguration;
        this->dq_initial = robotSettings.initialGeneralizedVelocity;
        this->dq_final   = robotSettings.finalGeneralizedVelocity;


        this->q = VectorXr::Zero(n,1);
        this->dq = VectorXr::Zero(n,1);
        this->ddq = VectorXr::Zero(n,1);
        this->Tau = VectorXr::Zero(n,1);

        this->constraintsValue = VectorXr::Zero(robotSettings.numberConstraints);

        switch (robotSettings.DifferentiationWRT) {

            case wrt_state: {
                this->g = MatrixXr::Zero(1,2*n);
                this->H = MatrixXr::Zero(2*n,2*n);
                this->JacobianConstraints = MatrixXr::Zero(robotSettings.numberConstraints,2*n);

                this->Diff_Tau = MatrixXr::Zero(n,2*n);
                this->DDiff_Tau = MatrixXr::Zero(n,2*n*2*n);

                //! Variables for differentiation
                this->D_X  = MatrixXr::Identity(2*n,2*n);
                this->D_q  = D_X.block(0,0,n,2*n);
                this->D_dq = D_X.block(n,0,n,2*n);

                break;
            }
            case wrt_controlPoints: {
                short int ncp = robotSettings.n*robotSettings.numberControlPoints;
                this->g = MatrixXr::Zero(1,ncp);
                this->H = MatrixXr::Zero(ncp,ncp);
                this->JacobianConstraints = MatrixXr::Zero(robotSettings.numberConstraints,ncp);

                this->Diff_Tau = MatrixXr::Zero(n,ncp);
                this->DDiff_Tau = MatrixXr::Zero(n,ncp*ncp);

                //! Variables for differentiation
                this->D_q  = MatrixXr::Zero(n,ncp);
                this->D_dq = MatrixXr::Zero(n,ncp);

                break;
            }
            default:

                cout<<"Unrecognized differentiation settings"<<endl;
                break;

        }


    }


    ~DirectCollocation(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    protected:

    //! Dof
    short int n;

    //! MultiBody
    std::shared_ptr< MultiBody > robot;

    //! Robot dynamics object
    InverseDynamics robotMotion;

    //! Settings structure for NLP
    robotSettingsTrajectoryOptimization robotSettings;

    //! Time vector
    MatrixXr S;

    //! Cost function
    real_t f;

    //! Gradient of cost function
    MatrixXr g;

    //! Hessian of cost function
    MatrixXr H;

    //! Constraints
    VectorXr constraintsValue;

    //! Jacobian of constraints
    MatrixXr JacobianConstraints;

    //! basis function
    MatrixXr B;

    //! first derivative of B wrt s
    MatrixXr dB;

    //! second derivative of B wrt s
    MatrixXr ddB;

    //! initial configuration
    VectorXr q_initial;

    //! final configuration
    VectorXr q_final;

    //! initial generalized velocity
    VectorXr dq_initial;

    //! final generalized velocity
    VectorXr dq_final;

    //! Generalized vectors
    VectorXr q, dq, ddq, Tau;

    //! Variables for differentiation
    MatrixXr D_X, D_q, D_dq, D_ddq, Diff_Tau, DDiff_Tau;


    public:



    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    /*! Description: Builds the basis functions and its first two derivatives, i.e. B, dB, ddB
         * \parameters: Vector time S and the numer of control points per DoF N
         * \return: int 0
         */
    int buildBasisFunctions(  );    

    /*! Description: Builds the basis functions and its first two derivatives, i.e. B, dB, ddB
         * \parameters: Vector time S and the numer of control points per DoF N
         * \return: int 0
         */
    int buildBasisFunctions( short int N, VectorXr S );

    //! Compute the objective function
         /*! \param control points c and boolean flags for computing partial derivatives
         * \return void
         */
    void computeObjectiveFunction(const MatrixXr &c, const bool &computeFirstDerivative, const bool &computeSecondDerivative);

    //! Compute constraints of NLP
         /*! \param control points c and boolean flag for computing partials
         * \return void
         */
    void computeConstraints(const MatrixXr &c, const bool &computePartialDerivatives);

    //! Compute constraints 2 of NLP
         /*! \param control points c and boolean flag for computing partials
         * \return void
         */
    void computeConstraintsII(const MatrixXr &c, const bool &computePartialDerivatives);

    //! Set initial configuration
         /*! \param initial configuration
         * \return void
         */
    void setInitialConfiguration( const VectorXr &q_initial ){ this->q_initial = q_initial; }

    //! Set final configuration
         /*! \param final configuration
         * \return void
         */
    void setFinalConfiguration( const VectorXr &q_final ){ this->q_final = q_final; }

    //! Set initial generalized velocity
         /*! \param initial generalized velocity
         * \return void
         */
    void setInitialGeneralizedVelocity( const VectorXr &dq_initial ){ this->dq_initial = dq_initial; }

    //! Set final generalized velocity
         /*! \param final generalized velocity
         * \return void
         */
    void setFinalGeneralizedVelocity( const VectorXr &dq_final ){ this->dq_final = dq_final; }

    //! Get the cost function value
    /*! \param none
        \return a real_t precision number
             */
    real_t getCost(){ return f; }

    //! Get the cost function gradient
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getCostGradient(){ return g; }

    //! Get the cost function Hessian
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getCostHessian(){ return H; }

    //! Get the constraints value
    /*! \param none
        \return a real vector
             */
    VectorXr getConstraints(){ return constraintsValue; }

    //! Get the constraints Jacobian
    /*! \param none
        \return a MatrixXr element
             */
    MatrixXr getConstraintsJacobian(){ return JacobianConstraints; }

    //! Get the basis function
    /*! \param none
        \return Matrix real_t
             */
    MatrixXr getBasis(){ return B; }

    //! Get the basis function first derivative
    /*! \param none
         \return Matrix real_t
              */
    MatrixXr getDBasis(){ return dB; }

    //! Get the basis function second derivative
    /*! \param none
          \return Matrix real_t
               */
    MatrixXr getDDBasis(){ return ddB; }


protected:


};
} // end of namespace core
} // end of namespace hr

#endif // HR_TRAJECTORY_OPTIMIZATION_DIRECT_COLLOCATION_H
