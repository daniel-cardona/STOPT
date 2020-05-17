#include "RobOptTraj/directTranscription.h"


int main(){

    OCSolver::Problem problem;
    OCSolver::Alg     algorithm;

    problem.nStates=2;
    problem.nControls=1;
    problem.nEvents=4;
    problem.nPath=2;

    problem.nNodes=10;
    problem.interpolationPoints=20;

    problem.discretizationMethod="Trapezoidal";

    problem.setup();

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&-----------------------Variable bounds-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    double xL = 0.0;
    double yL = 0.0;
    double xU = 2.0;
    double yU = 2.0;

    double thetaL = -10.0;
    double thetaU = 10.0;

    double x0 = 0.0;
    double y0 = 0.0;
    double xf = 1.2;
    double yf = 1.6;


    problem.bounds.states.lower(0)=xL;
    problem.bounds.states.upper(0)=xU;

    problem.bounds.states.lower(1)=yL;
    problem.bounds.states.upper(1)=yU;

    problem.bounds.controls.lower(0)=thetaL;
    problem.bounds.controls.upper(0)=thetaU;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------Boundary and path constraints-------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    problem.bounds.events.lower(0)=x0;
    problem.bounds.events.lower(1)=y0;
    problem.bounds.events.lower(2)=xf;
    problem.bounds.events.lower(3)=yf;

    problem.bounds.events.upper(0)=x0;
    problem.bounds.events.upper(1)=y0;
    problem.bounds.events.upper(2)=xf;
    problem.bounds.events.upper(3)=yf;

    problem.bounds.path.lower(0)=0.1;
    problem.bounds.path.upper(0)=100.0;

    problem.bounds.path.lower(1)=0.1;
    problem.bounds.path.upper(1)=100.0;

    problem.bounds.startTime.lower(0)=0;
    problem.bounds.startTime.upper(0)=0;

    problem.bounds.finalTime.lower(0)=0;
    problem.bounds.finalTime.upper(0)=10;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------------Set initial guess-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    int nNodes=problem.nNodes;

    problem.guess.states.row(0).setLinSpaced(nNodes,x0,xf);
    problem.guess.states.row(1).setLinSpaced(nNodes,y0,yf);

    problem.guess.controls.row(0).setZero();

    problem.guess.time.row(0).setLinSpaced(nNodes,0,5);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------------Function pointers-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    problem.integrand_cost=&OCSolver::integrand_cost;
    problem.endpoint_cost=&OCSolver::endpoint_cost;
    problem.dae=&OCSolver::dae;
    problem.events=&OCSolver::events;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------------Algorithm options-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    algorithm.nlp_iter_max          =1000;
    algorithm.nlp_tolerance         =1.e-4;
    algorithm.scaling               ="none";
    algorithm.derivatives           ="numerical";
    algorithm.defect_scaling        ="state-based";
    algorithm.derivativeChecker     =true;
    algorithm.sparsityOptimization  =true;


    OCSolver::core::robTrajSol(problem,algorithm);



};
