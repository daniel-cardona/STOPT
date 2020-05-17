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

namespace NLP
{

    void determineScalingFactorsDecVariables(Problem &problem,Alg &algorithm){

            int nControls=problem.nControls;
            int nStates=problem.nStates;

            double zlower, zupper;

            Eigen::VectorXd control_scaling(nControls);
            Eigen::VectorXd state_scaling(nStates);
            double time_scaling;

                //Control scaling

                control_scaling.setOnes();

                if(algorithm.scaling=="automatic"){

                    for(int i=0;i<nControls;i++){

                        zlower=problem.bounds.controls.lower(i);
                        zupper=problem.bounds.controls.upper(i);

                        if(zlower !=-INF && zupper!=INF){
                            if(zlower !=0.0 || zupper !=0.0){
                                control_scaling(i)=1.0/MAX(fabs(zlower),fabs(zupper));
                            }
                        }

                        else if(zlower==-INF && zupper!=INF && zupper!=0.0){
                            control_scaling(i)=1.0/fabs(zupper);
                        }

                        else if(zupper==INF && zlower!=-INF && zlower!=0.0){
                            control_scaling(i)=1.0/fabs(zlower);
                        }

                    }
                }

                //State scaling

                state_scaling.setOnes();

                if(algorithm.scaling=="automatic"){

                    for(int i=0;i<nStates;i++){

                        zlower=problem.bounds.states.lower(i);
                        zupper=problem.bounds.states.upper(i);

                        if(zlower !=-INF && zupper!=INF){
                            if(zlower !=0.0 || zupper !=0.0){
                                state_scaling(i)=1.0/MAX(fabs(zlower),fabs(zupper));
                            }
                        }

                        else if(zlower==-INF && zupper!=INF && zupper!=0.0)
                            state_scaling(i)=1.0/fabs(zupper);

                        else if(zupper==INF && zlower!=-INF && zlower!=0.0)
                            state_scaling(i)=1.0/fabs(zlower);
                    }
                }

                //Time scaling

                time_scaling=1;

                if(algorithm.scaling=="automatic"){

                        zlower=problem.bounds.startTime.lower(0);
                        zupper=problem.bounds.finalTime.upper(0);

                        if(zlower !=-INF && zupper!=INF){
                            if(zlower !=0.0 || zupper !=0.0){
                                time_scaling=1.0/MAX(fabs(zlower),fabs(zupper));
                            }
                        }

                        else if(zlower==-INF && zupper!=INF && zupper!=0.0)
                            time_scaling=1.0/fabs(zupper);

                        else if(zupper==INF && zlower!=-INF && zlower!=0.0)
                            time_scaling=1.0/fabs(zlower);
                    }

                //--------------------------------------------------------------------------------------------------------------------------
                //Obtain a matrix with the scaling factors of the decision variables and obtain the inverse (ANALYTICAL DERIVATIVES FEATURE)

                Eigen::VectorXd scaleVector(problem.nDecVar);
                Eigen::MatrixXd scaleMatrix(problem.nDecVar,problem.nDecVar);

                problem.scale.invScaleMatrix.resize(problem.nDecVar,problem.nDecVar);

                scaleVector(0)=time_scaling;
                scaleVector(1)=time_scaling;

                for(int k=0;k<problem.nNodes;k++){

                   scaleVector.segment(2+(k*(nStates+nControls)),nStates+nControls)<<state_scaling,control_scaling;

                }

                scaleMatrix=scaleVector.asDiagonal();

                //--------------------------------------------------------------------------------------------------------------------------



                //Return the scale factors!
                problem.scale.controls=control_scaling;
                problem.scale.states=state_scaling;
                problem.scale.time=time_scaling;
                problem.scale.invScaleMatrix=scaleMatrix.inverse().sparseView();


        }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    void determineObjectiveScaling(Problem &problem, Alg &algorithm){

         Eigen::VectorXd x0(problem.nDecVar);

         Eigen::VectorXd grad(problem.nDecVar);

         double enorm;

         x0=problem.x0;

         problem.scale.cost_fcn=1;

         if(algorithm.scaling=="automatic"){

            //if(algorithm.derivatives=="numerical"){

                problem.scale.cost_fcn=-1.0;


                if(problem.discretizationMethod=="Trapezoidal")
                    derivatives::numerical::scalarGradient(collocationMethod::trapezoidal::costFunction,x0,problem,algorithm,grad);

                //}

             enorm=grad.norm();

             if(enorm!=0.0 &&enorm<INF){
                 problem.scale.cost_fcn=1/enorm;
                }

              else
                 problem.scale.cost_fcn=1;

          }

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    void determineConstraintScaling(Problem &problem,Alg &algorithm){

        int nDecVar=problem.nDecVar;
        int nCns=problem.nCns;

        Eigen::VectorXd jac_row_norm(nCns);

        problem.cnsScaling.setOnes();

            if (algorithm.scaling=="automatic"){

                //if(algorithm.derivatives=="numerical"){

                        Eigen::VectorXd jacCol(nCns);

                        jac_row_norm.setZero();
                        jacCol.setZero();

                        problem.use_constraint_scaling=false;

                        for(int i=0;i<nDecVar;i++){

                            if(problem.discretizationMethod=="Trapezoidal"){

                                derivatives::numerical::jacobianColumn(collocationMethod::trapezoidal::cnsFunction,problem.x0,problem,algorithm,i,jacCol);

                            }

                            jac_row_norm.array()+=(jacCol.array().square());

                        }

                        jac_row_norm.array()=jac_row_norm.array().sqrt();

                //}

                //Obtain the scaling factors of the constraints

                double sqeps=sqrt(MC_EPSILON);

                for(int i=0;i<nCns;i++){

                    if(jac_row_norm(i)<1e7){

                        problem.cnsScaling(i)=1.0/(jac_row_norm(i)+sqeps);

                    }

                    else {

                        if(jac_row_norm(i)>1.e7){

                            problem.cnsScaling(i)=1.0/1.e7;
                        }

                        if(jac_row_norm(i)==0.0){

                            problem.cnsScaling(i)=1.0;
                        }

                    }


                }

               }

                //---------Set variables special variables used in the analytical differientation methods

               //Obtain the scaling factor for the event constraints

                problem.scale.events=problem.cnsScaling.segment(problem.nDefectCns,problem.nEvents);

                //Obtain the scaling factor for the path constraints

                Eigen::Map<Eigen::MatrixXd> pathCns(problem.cnsScaling.segment(problem.nDefectCns+problem.nEvents,problem.nPath*problem.nNodes).data(),problem.nPath,problem.nNodes);

                problem.scale.path=pathCns;

                //Obtain the scaling factor of the time constraint

                problem.scale.timeCns=problem.cnsScaling(problem.nCns-1);

              //--------------------------------------------------------------------------------------------------


                //Set flag for using scaling constraints

                problem.use_constraint_scaling=true;


                //For the defects constraints is  used a Jacobian based scaling?

                if(algorithm.defect_scaling=="Jacobian-based"){

                    return;
                }


                //By default, use the state scaling factors for the differential defects


                Eigen::VectorXd state_scaling(problem.nStates);

                int nOrder=problem.nSegments;

                state_scaling=problem.scale.states;

                int offset,l=0;

                for(int k=0;k<nOrder+1;k++){

                    if(k!=nOrder){

                        for(int j=0;j<problem.nStates;j++){

                            l=(k*problem.nStates)+j;

                            problem.cnsScaling(l)=state_scaling(j);

                            problem.scale.defects(j,k)=state_scaling(j);

                            problem.scale.defectsV(l)=state_scaling(j);
                       }

                    }

                }







    }


}//END NLP namespace

}//END OCSolver namespace


