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
#include <chrono>


using namespace  std::chrono;

namespace  OCSolver
{

namespace core {


   void robTrajSol(Problem &problem,Alg &algorithm){

       //Let's validate the user input!

       try {
            validateUserInput(problem,algorithm);
           }
       catch (const char* msg) {
           cerr<<msg<<endl;
           return;
           }

        //Display the information of the problem

          displayInformation(problem,algorithm);

        //1. Obtain the initial guess and obtain the scaling factors of the variables, bounds, constraints and cost function

          NLP::setInitialGuess(problem,algorithm);

        //2. Set the variables bounds of the NLP Problem

          NLP::setVariableBounds(problem,algorithm);

        //3. Set the constraints bounds of the NLP Problem

          NLP::setConstraintBounds(problem,algorithm);


            //-------------------------------------------
            NLP::sparsity::generateSparsityTemplates(problem,algorithm);
            NLP::sparsity::detectSparsity(problem);
            NLP::sparsity::getIndexGroups(problem);
            NLP::sparsity::getIndexGroupsDerivativeMatrix(problem);
            NLP::sparsity::getPerturbationMatrix(problem,algorithm);
            //-------------------------------------------




        //9. Solve the problem!

            NLP::ipopt::solve(problem,algorithm);



    }

   void validateUserInput(Problem &problem,Alg &algorithm){

       if(problem.discretizationMethod != "Trapezoidal" && problem.discretizationMethod != "Hermite-Simpson"){

           throw "Error: Invalid discretization method! Available methods: Trapezoidal and Hermite-Simpson";

       }

       if(algorithm.derivatives != "numerical" && algorithm.derivatives!= "analytical"){

           throw "Error: Invalid derivative method! Available options: numerical and analytical";

       }

       if(algorithm.scaling !="automatic" && algorithm.scaling !="none"){

           throw "Error: Invalid option for scaling algorithm! Available options: automatic and none";

       }


       if(problem.nControls>0 && (problem.bounds.controls.lower.array()>problem.bounds.controls.upper.array()).any()){

               throw "Error:Infeasible control variable bounds supplied by the user! Lower bound > Upper bound";
        }


       if(problem.nStates>0 && (problem.bounds.states.lower.array()>problem.bounds.states.upper.array()).any()){

               throw "Error:Infeasible state variable bounds supplied by the user! Lower bound > Upper bound";
        }

       if(problem.nEvents>0 && (problem.bounds.events.lower.array()>problem.bounds.events.upper.array()).any()){

               throw "Error:Infeasible events bounds supplied by the user! Lower bound > Upper bound";
        }


       if(problem.nEvents>0 && (problem.bounds.events.lower.array()>problem.bounds.events.upper.array()).any()){

               throw "Error:Infeasible events bounds supplied by the user! Lower bound > Upper bound";
        }

       if( (problem.bounds.startTime.lower.array()>problem.bounds.startTime.upper.array()).any()){

               throw "Error:Infeasible start time bounds supplied by the user! Lower bound > Upper bound";
        }

       if( (problem.bounds.finalTime.lower.array()>problem.bounds.finalTime.upper.array()).any()){

               throw "Error:Infeasible final time bounds supplied by the user! Lower bound > Upper bound";
        }

   }

   void displayInformation(Problem &problem,Alg &algorithm){

       cout<<"                                                                                      "<<endl;
       cout<<"---------------------------Using direct collocation solver----------------------------"<<endl;
       cout<<"                                                                                      "<<endl;
       cout<<"**************************************************************************************"<<endl;
       cout<<"**************************** Transcription information *******************************"<<endl;
       cout<<"**************************************************************************************"<<endl;
       cout<<endl;
       cout<<"Discretization method: "<<problem.discretizationMethod<<endl;
       cout<<"Number of discretization points set by user:"<<problem.nNodes<<endl;
       cout<<"Number of segments generated: "<<problem.nSegments<<endl;
       cout<<"Number of NLP variables: "<<problem.nDecVar<<endl;
       cout<<"Number of NLP constraints: "<<problem.nCns<<endl;
       cout<<"Using "<<algorithm.derivatives<<" derivatives"<<endl;



   }







}


}


