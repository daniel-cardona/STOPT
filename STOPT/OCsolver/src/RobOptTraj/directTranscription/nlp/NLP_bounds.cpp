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


void setVariableBounds(Problem &problem, Alg &algorithm){

    //-------Index Variables-------------

    int nOrder=problem.nSegments;
    int nStates=problem.nStates;
    int nControls=problem.nControls;
    int offset=0;

    //---------Scaling variables---------

    Eigen::VectorXd stateScaling(nStates);
    Eigen::VectorXd controlScaling(nControls);
    double timeScaling;

    stateScaling=problem.scale.states;
    controlScaling=problem.scale.controls;
    timeScaling=problem.scale.time;

    //----User defined variable bounds----

    Eigen::VectorXd states_lb(nStates);
    Eigen::VectorXd controls_lb(nControls);
    Eigen::VectorXd time_lb(2);

    Eigen::VectorXd states_ub(nStates);
    Eigen::VectorXd controls_ub(nControls);
    Eigen::VectorXd time_ub(2);

    states_lb=problem.bounds.states.lower;
    controls_lb=problem.bounds.controls.lower;

    states_ub=problem.bounds.states.upper;
    controls_ub=problem.bounds.controls.upper;

    time_lb(0)=problem.bounds.startTime.lower(0);
    time_lb(1)=problem.bounds.finalTime.lower(0);

    time_ub(0)=problem.bounds.startTime.upper(0);
    time_ub(1)=problem.bounds.finalTime.upper(0);

    Eigen::VectorXd lbk;
    Eigen::VectorXd ubk;

    //Define the structure of the variable bounds depending on the defect chosed by the user

    if(problem.discretizationMethod=="Trapezoidal"){

        offset=problem.nStates+problem.nControls;

        lbk.setZero(offset);
        ubk.setZero(offset);

        lbk<<states_lb.cwiseProduct(stateScaling),controls_lb.cwiseProduct(controlScaling);
        ubk<<states_ub.cwiseProduct(stateScaling),controls_ub.cwiseProduct(controlScaling);

    }

    if(problem.discretizationMethod=="Hermite-Simpson"){

        offset=problem.nStates+problem.nControls*2;

        lbk.setZero(offset);
        ubk.setZero(offset);

        lbk<<states_lb.cwiseProduct(stateScaling),controls_lb.cwiseProduct(controlScaling),controls_lb.cwiseProduct(controlScaling);
        ubk<<states_ub.cwiseProduct(stateScaling),controls_ub.cwiseProduct(controlScaling),controls_ub.cwiseProduct(controlScaling);

    }

    //Set the variable bounds vectors

    Eigen::VectorXd xlb(problem.xlb.size());
    Eigen::VectorXd xub(problem.xlb.size());

    for(int k=0;k<nOrder+1;k++){

        //If a Hermite-Simpson method is used

        if(problem.discretizationMethod=="Hermite-Simpson"){

            if(k!=nOrder){
                xlb.segment(k*offset+2,offset)=lbk;
                xub.segment(k*offset+2,offset)=ubk;
            }

            else{
                xlb.segment(k*offset+2,offset-nControls)<<states_lb.cwiseProduct(stateScaling),controls_lb.cwiseProduct(controlScaling);
                xub.segment(k*offset+2,offset-nControls)<<states_ub.cwiseProduct(stateScaling),controls_ub.cwiseProduct(controlScaling);

            }
        }

        //If a Trapezoidal method is used
        else{

            xlb.segment(k*offset+2,offset)=lbk;
            xub.segment(k*offset+2,offset)=ubk;

        }


    }

    //Set time bounds

    xlb(0)=time_lb(0)*timeScaling;
    xub(0)=time_ub(0)*timeScaling;

    xlb(1)=time_lb(1)*timeScaling;
    xub(1)=time_ub(1)*timeScaling;

    //Return this

    problem.xlb=xlb;
    problem.xub=xub;

}




void setConstraintBounds(Problem &problem,Alg &algorithm){


    Eigen::VectorXd cns_scaling(problem.nCns);

    cns_scaling=problem.cnsScaling;

    int nStates=problem.nStates;
    int nEvents=problem.nEvents;
    int nPath=problem.nPath;
    int nOrder=problem.nSegments;

    double event_sc, path_sc;

    int offset;
    int j;

    //Set the equalities bounds for defects constraint

    for(int i=0;i<problem.nDefectCns;i++){

        problem.glb(i)=0.0;
        problem.gub(i)=0.0;

    }

    offset=problem.nDefectCns;

    //Bound of the event constraints

    for(int k=0;k<nEvents;k++){

        j=offset+k;

        event_sc=cns_scaling(j);

        problem.glb(j)=problem.bounds.events.lower(k)*event_sc;
        problem.gub(j)=problem.bounds.events.upper(k)*event_sc;


    }

    offset+=nEvents;

    //Bounds of the path constraints

    for(int k=0;k<nPath;k++){

        for(int l=0;l<(nOrder+1);l++){

            j=offset+(l*nPath)+k;

            path_sc=cns_scaling(j);

            problem.glb(j)=problem.bounds.path.lower(k)*path_sc;
            problem.gub(j)=problem.bounds.path.upper(k)*path_sc;

        }

    }

    if(problem.discretizationMethod=="Hermite-Simpson"){

        offset+=nPath*(nOrder+1);

        for(int k=0;k<nPath;k++){

            for(int l=0;l<(nOrder);l++){

                j=offset+(l*nPath)+k;

                path_sc=cns_scaling(j);


                problem.glb(j)=problem.bounds.path.lower(k)*path_sc;
                problem.gub(j)=problem.bounds.path.upper(k)*path_sc;

            }

        }

    }


    //Bound for t0<=tf constraint

    double diff_t0Min_tfMax=problem.bounds.startTime.lower(0)-problem.bounds.finalTime.upper(0);

    diff_t0Min_tfMax*=problem.scale.time;

    diff_t0Min_tfMax*=cns_scaling(problem.nCns-1);

    problem.glb(problem.nCns-1)=diff_t0Min_tfMax;
    problem.gub(problem.nCns-1)=0.0;


}



} //END NLP namespace

} //END OCSolver namespace


