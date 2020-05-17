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

namespace trapezoidal
{


void cnsFunction(Eigen::VectorXd &x,Eigen::VectorXd &cns,Problem &problem,Alg &algorithm){

    //This function implements the NLP inequality constraints for numerical or analytic differentiaton

       int nOrder=problem.nSegments;
       int nStates=problem.nStates;
       int nControls=problem.nControls;

       int nEvents=problem.nEvents;
       int nPath= problem.nPath;

       double resid=0.0;
       int offset=0;
       int path_offset=nStates*(nOrder)+nEvents; //The constraints are located -> g=[Defects Events Path]'

       double t0,tf;

       Eigen::MatrixXd decVarMatrix(problem.nStates+problem.nControls,problem.nNodes);

       Eigen::VectorXd controls(nControls);
       Eigen::VectorXd states(nStates);

       Eigen::VectorXd dx(nStates);
       Eigen::VectorXd path(nPath);

       Eigen::VectorXd controls_next(nControls);
       Eigen::VectorXd states_next(nStates);

       Eigen::VectorXd dx_next(nStates);
       Eigen::VectorXd path_next(nPath);

       Eigen::VectorXd g(problem.nCns);

       g.setZero(problem.nCns);
       cns.setZero(problem.nCns);

       //Obtain the decision variable vector in matrix form

       utils::getDecVarMatrix(problem,algorithm,x,decVarMatrix,t0,tf);


       //-----------Differential defects

       for(int k=0;k<nOrder+1;k++){  //Iterate along the knot points

           //Obtain the controls,states and time in k

           utils::getVariables(problem,decVarMatrix,states,controls,k);

           double tk=utils::convert_to_original_time(problem.snodes(k),t0,tf);

           problem.dae(states,controls,tk,dx,path,problem);

           //-----------Trapezoidal method---------------

               if(k!=(nOrder)){

                        //Obtain the controls,states and time in k+1

                        double tk1=utils::convert_to_original_time(problem.snodes(k+1),t0,tf);

                        double hk=(tf-t0)/problem.nSegments;

                        utils::getVariables(problem,decVarMatrix,states_next,controls_next,k+1);

                        //Evaluate the daes in k+1

                        problem.dae(states_next,controls_next,tk1,dx_next,path_next,problem);

                        //Obtain the value of the diferential constraint

                        for(int j=0;j<nStates;j++){

                            resid=states_next(j)-states(j)-hk*(dx(j)+dx_next(j))/2.0;

                            offset=(k*nStates)+j;

                            g(offset)=resid;//*(tf-t0)/(2.0*hk);
                           }

                }


               //--------- Path constraints

               for(int j=0;j<nPath;j++){

                   offset=path_offset+(k*nPath)+j;
                   g(offset)=path(j);

                   }

       }

        //-------------------End for k-----------------

       //---------------Event constraints-----------------

       Eigen::VectorXd initial_states(problem.nStates);
       Eigen::VectorXd final_states(problem.nStates);

       Eigen::VectorXd e(nEvents);

       utils::getVariables(problem,decVarMatrix,initial_states,controls,0);

       utils::getVariables(problem,decVarMatrix,final_states,controls,problem.nNodes-1);

       problem.events(initial_states,final_states,t0,tf,e,problem);


       for(int k=0; k<nEvents;k++){

           offset=nStates*(nOrder)+k;
           g(offset)=e(k);

       }


       //Time constraint

       g(problem.nCns-1)=(t0-tf)*problem.scale.time;


       //If the function is called from the scaling factors just return the value of g(X)

       //Else scale the constraints
       if(problem.use_constraint_scaling!=false){

           cns=g.cwiseProduct(problem.cnsScaling);

       } else{

           cns=g;
       }





}


void rightHandSideSparsity(Eigen::VectorXd &x,Eigen::VectorXd &cns,Problem &problem, Alg &algorithm){

    //This function exploits the separability in the trapezoidal defect constraint function in order to used the sparsity properties of the functions
    //For more information refer to: Betts(2016)

        //Index variables

        int nSegments=problem.nSegments;
        int nStates=problem.nStates;
        int nControls=problem.nControls;

        int offset=0;

        int pathOffset=nStates*(nSegments+1)+problem.nEvents;

        int nEvents=problem.nEvents;
        int nPath= problem.nPath;

        //Defect constraint variables

        Eigen::VectorXd residk(nStates);
        Eigen::VectorXd residk1(nStates);

        double t0,tf;

        Eigen::MatrixXd decVarMatrix;

        Eigen::VectorXd controls(nControls);
        Eigen::VectorXd states(nStates);

        Eigen::VectorXd dx(nStates);
        Eigen::VectorXd path(nPath);

        Eigen::VectorXd controls_next(nControls);
        Eigen::VectorXd states_next(nStates);

        Eigen::VectorXd dx_next(nStates);
        Eigen::VectorXd path_next(nPath);

        Eigen::VectorXd controls_mid(nControls);
        Eigen::VectorXd states_mid(nStates);
        Eigen::VectorXd dx_mid(nStates);
        Eigen::VectorXd path_mid(nPath);

        //Main Vector of the separability
        Eigen::VectorXd q(nStates*(nSegments+1)+nEvents+( nPath*(nSegments+1) ) +1);
        q.setZero();

        //Obtain the decision variable vector in matrix form

        utils::getDecVarMatrix(problem,algorithm,x,decVarMatrix,t0,tf);

        //Differential defects


        for(int k=0;k<nSegments+1;k++){

            //Obtain the states, controls and time in the k node and evaluate the dae

            utils::getVariables(problem,decVarMatrix,states,controls,k);

            double tk=utils::convert_to_original_time(problem.snodes(k),t0,tf);

            problem.dae(states,controls,tk,dx,path,problem);

            double hk=(tf-t0)/problem.nSegments;

            residk=(hk*dx);
            offset=(k*nStates);

            q.segment(offset,nStates)<<residk;

            //---------------Path Constraints-----------------

            for(int j=0;j<nPath;j++){

                offset=pathOffset+(k*nPath)+j;
                q(offset)=path(j)*problem.scale.path(j,k);

            }

        }


        //---------------Event constraints-----------------

        Eigen::VectorXd initial_states(problem.nStates);
        Eigen::VectorXd final_states(problem.nStates);

        Eigen::VectorXd e(nEvents);

        utils::getVariables(problem,decVarMatrix,initial_states,controls,0);

        utils::getVariables(problem,decVarMatrix,final_states,controls,problem.nNodes-1);


        //Obtain the value of the event constraints and scale it

        problem.events(initial_states,final_states,t0,tf,e,problem);

        for(int k=0; k<nEvents;k++){

            offset=nStates*(nSegments+1)+k;
            q(offset)=e(k)*problem.scale.events(k);

        }


        //Obtain the value of the time constraint and scale it

        q(q.size()-1)=(t0-tf)*problem.scale.time*problem.scale.timeCns;

        cns=q;


}


}//END trapezoidalMethod
}//END NLP namespace
}//END OCSolver namespace

