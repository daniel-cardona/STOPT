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

namespace  OCSolver
{

void Problem::setup(){

    //=========== 1. NLP Setup: Compute the number of NLP variables and constraints generated for the transcription process =================

    nSegments=nNodes-1;                         //Number of segments
    nDecVar=(nControls+nStates)*(nNodes)+2;     //Number of decision variables
    nDefectCns=nStates*nSegments;               //Number of defect constraints
    nCns=nDefectCns+nEvents+nPath*(nNodes)+1;   //Total number of constraints

    //If any HS is used consider the midpoint variables and constraints

    if(discretizationMethod=="Hermite-Simpson"){

      nDecVar+=(nControls)*(nNodes);            //Add the mid point variables

      nCns+=nPath*nSegments;                    //Add the mid point constraints

    }

    //Set NULL pointer to all functions (Just for default)
        endpoint_cost=NULL;
        integrand_cost=NULL;
        dae=NULL;
        events=NULL;


    //===========2. NLP elements initialization: Set the sizes of the OCP data vectors =================


        //State bounds
        bounds.states.lower.setZero(nStates);
        bounds.states.upper.setZero(nStates);

        //Control bounds
        bounds.controls.lower.setZero(nControls);
        bounds.controls.upper.setZero(nControls);

        //Event constraints bounds
        bounds.events.lower.setZero(nEvents);
        bounds.events.upper.setZero(nEvents);

        //Path constraints bounds
        bounds.path.lower.setZero(nPath);
        bounds.path.upper.setZero(nPath);

        //Time bounds
        bounds.startTime.lower.setZero(1);
        bounds.startTime.upper.setZero(1);

        bounds.finalTime.lower.setZero(1);
        bounds.finalTime.upper.setZero(1);

        //Initial guess state vector
        guess.states.setZero(nStates,nNodes);

        //Initial guess control vector
        guess.controls.setZero(nControls,nNodes);

        //Initial guess time vector
        guess.time.setZero(1,nNodes);

        //Scaling vectors
        scale.controls.setZero(nControls);
        scale.states.setZero(nStates);
        scale.defects.setZero(nStates,nNodes-1);
        scale.defectsV.setZero(nDefectCns);
        scale.path.setZero(nPath,nNodes);
        scale.events.setZero(nEvents);

        //Solution vectors
        solution.controls.setZero(nControls,interpolationPoints);
        solution.states.setZero(nStates,interpolationPoints);
        solution.time.setZero();

      //Set the sizes of the NLP vectors

        //Time vector
        snodes.setLinSpaced(nNodes,-1.0,1.0);

        //Initial guess vector
        x0.setZero(nDecVar);

        //Decision variables lower bound
        xlb.setZero(nDecVar);

        //Decision variables upper bounds
        xub.setZero(nDecVar);

        //Constraint scaling
        cnsScaling.setZero(nCns);

        //Lower bounds of the constraints
        glb.setZero(nCns);

        //Upper bounds of the constraints
        gub.setZero(nCns);

        //===========3. Initialization of the function variables =================


        if(discretizationMethod=="Trapezoidal"){

            decVarMatrix.setZero(nStates+nControls,nNodes);
        }

           kStates.setZero(nStates);
           kControls.setZero(nControls);
           dx.setZero(nStates);
           kPath.setZero(nPath);

        //===========4.Initialization of the derivatives variables=============

           dxGradient.setZero(nStates,nStates+nControls);
           pathGradient.setZero(nPath,nStates+nControls);
           eventGradient_t0.setZero(nEvents,nStates+1);
           eventGradient_tF.setZero(nEvents,nStates+1);
           derivativeMatrix.setZero(nNodes*(nStates+nPath)+nEvents+1,nDecVar);

           D.resize(nCns,nDecVar);
           J.resize(nNodes*(nStates+nPath)+nEvents+1,nDecVar);




}

}


