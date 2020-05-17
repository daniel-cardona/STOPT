#include "RobOptTraj/directTranscription.h"


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

namespace OCSolver
{

namespace collocationMethod
{


namespace trapezoidal {

    double costFunction(Eigen::VectorXd &x,Problem &problem,Alg &algorithm){

        //This functions implements the NLP cost function for analytic or numerical differentiation

        double retval=0;

        int nOrder=problem.nSegments;

        double t0,tf;

        int k_1; //Index variable

        double endpoint_cost=0;

        double sum_phase=0.0;

        double sum_cost=0.0;

        Eigen::MatrixXd decVarMatrix(problem.nStates+problem.nControls,problem.nNodes);

        Eigen::VectorXd controls(problem.nControls);

        Eigen::VectorXd states(problem.nStates);

        Eigen::VectorXd states_next(problem.nStates);

        utils::getDecVarMatrix(problem,algorithm,x,decVarMatrix,t0,tf);


        //Integrate the Lagrange term using trapezoidal rule

        for(int k=0;k<nOrder;k++){

            double interval_cost=0;

            double integrand;

            //***************

            k_1=k+1;

            //Get the unscaled states and controls in k segment

            utils::getVariables(problem,decVarMatrix,states,controls,k);

            //Get the unscaled states and controls in k+1 and save in the problem.state_k_1 variables

            utils::getVariables(problem,decVarMatrix,problem.states_k1,controls,k_1);

            //***************


            //Obtain the integration step

            double tk=utils::convert_to_original_time(problem.snodes(k),t0,tf);

            double tk1=utils::convert_to_original_time(problem.snodes(k+1),t0,tf);

            double h=tk1-tk;

            //Evaluate the integral part of the cost function in k

            interval_cost=problem.integrand_cost(states,controls,tk,problem);

            //Evaluate the integral part of the cost function in k+1

            //***************

            k_1=k+2;

            utils::getVariables(problem,decVarMatrix,states_next,controls,k+1);

            if(k_1>=problem.nNodes){ k_1=k+1;}

            utils::getVariables(problem,decVarMatrix,problem.states_k1,controls,k_1);

            //***************

            integrand=problem.integrand_cost(states_next,controls,tk1,problem);

            interval_cost+=integrand;

            interval_cost*=h/2.0;

            sum_phase+=interval_cost;

            }

        //Mayer term (end point cost evaluation)

        Eigen::VectorXd initial_states(problem.nStates);
        Eigen::VectorXd final_states(problem.nStates);

        sum_cost+=sum_phase;

        utils::getVariables(problem,decVarMatrix,initial_states,controls,0);              //Obtain the states in t0

        utils::getVariables(problem,decVarMatrix,final_states,controls,problem.nNodes-1); //Obtain the states in tF

        endpoint_cost=problem.endpoint_cost(initial_states,final_states,t0,tf,problem);

        sum_cost+=endpoint_cost;

        if(problem.scale.cost_fcn!=-1.0){

            retval=sum_cost*problem.scale.cost_fcn;

        }

        else
            retval=sum_cost;

        return retval;



    }

} //END trapezoidalMethod namespace

} //END NLP namespace

} //END OCSolver namespace

