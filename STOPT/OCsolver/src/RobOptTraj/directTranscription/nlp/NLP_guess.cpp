/*
 *
 * Copyright (C) 2019
 * Daniel S. Cardona <ingdanielcardona@gmail.com>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx>
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
 *	\file src/problem.h
 *	\author Daniel S. Cardona
 *	\version 1.0
 *	\date 2019
 *
 * Problem and algorithm class implementation.
 */

#include "RobOptTraj/directTranscription.h"

namespace OCSolver
{

namespace NLP {


        void setInitialGuess(Problem &problem, Alg &algorithm){

            int nSegments=problem.nSegments;
            int nControls=problem.nControls;
            int nStates=problem.nStates;

            double t00, tf0;

            Eigen::MatrixXd state_guess(nStates,nSegments+1);
            Eigen::MatrixXd control_guess(nControls,nSegments+1);

            //Obtain initial guess of t0 and tf

            t00=problem.guess.time(0,0);
            tf0=problem.guess.time(0,nSegments);

            //Obtain the user-defined guess for controls and states

            state_guess=problem.guess.states;
            control_guess=problem.guess.controls;

            //Obtain the scaling factor variables using the initial guess vector supplied by the user

            determineScalingFactorsDecVariables(problem,algorithm);

            //determineScalingFactorsDecVariables(problem,algorithm);

            //Scaling the x0

            t00*=problem.scale.time;
            tf0*=problem.scale.time;

            state_guess.array().colwise()*=problem.scale.states.array();

            control_guess.array().colwise()*=problem.scale.controls.array();


            //Vectorice the control and state matrix

            Eigen::Map<Eigen::VectorXd> vecControl(control_guess.data(),control_guess.size());
            Eigen::Map<Eigen::VectorXd> vecState(state_guess.data(),state_guess.size());


            //Now set the guess vector in the following order

            //Trapezoidal method  ---> [t0, tF, x1, u1, x2, u2...., xM, UM]
            //HSC method          ---> [t0, tF, x1, u1, u1.5, x2, u2, u2.5,...., xM, UM]
            //HSS method          ---> [t0, tF, x1, u1, x1.5, u1.5, x2, u2, x2.5, u2.5,...., xM, UM]

            Eigen::VectorXd x0(problem.x0.size());

            //If a trapezoidal defect is defined

            if(problem.discretizationMethod=="Trapezoidal"){

                Eigen::MatrixXd decVar(nStates+nControls,nSegments+1);

                decVar<<state_guess,control_guess;

                Eigen::Map<Eigen::VectorXd> vecX0(decVar.data(),decVar.size());

                problem.x0<<t00,tf0,vecX0;

            }

            //If a Hermite-Simpson Compressed defect is defined

            if(problem.discretizationMethod=="Hermite-Simpson"){

                Eigen::MatrixXd decVar(nStates+nControls*2,nSegments+1);

                decVar<<state_guess,control_guess,control_guess;

                Eigen::Map<Eigen::VectorXd> vecX0(decVar.data(),decVar.size());

                problem.x0<<t00,tf0,vecX0.head(vecX0.size()-nControls);

            }

            //Now obtain the scaling factors of the cost function

            determineObjectiveScaling(problem,algorithm);

            //Now obtain the scaling factor of the constraints functions

            determineConstraintScaling(problem,algorithm);


        }




} //END NLP namespace

} //END OCSolver namespace


