#include "RobOptTraj/directTranscription.h"



double OCSolver::endpoint_cost(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tf,Problem &problem){

    return 0;

}

double OCSolver::integrand_cost(Eigen::VectorXd &states,Eigen::VectorXd &controls,double &tk,Problem &problem){

    double V=2.138;

    double theta=controls(0);

    double dxdt=V*cos(theta);
    double dydt=V*sin(theta);

    double L=pow(dxdt,2.0)+ pow(dydt,2.0);

    return L;

}

void OCSolver::dae(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::VectorXd &derivatives,Eigen::VectorXd &path,Problem &problem){

    double x =states(0);
    double y =states(1);

    double theta=controls(0);

    double V=2.138;


    double dxdt= V*cos(theta);
    double dydt= V*sin(theta);

    derivatives.setZero(states.size());

    derivatives(0)=dxdt;
    derivatives(1)=dydt;

    path.setZero(2);

    path(0)=pow(x-0.4,2.0)+pow(y-0.5,2.0);
    path(1)=pow(x-0.8,2.0)+pow(y-1.5,2.0);


}

void OCSolver::events(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::VectorXd &e,Problem &problem){

    double x0=initial_states(0);
    double y0=initial_states(1);

    double xf=final_states(0);
    double yf=final_states(1);

    e.setZero(4);

    e(0)=x0;
    e(1)=y0;
    e(2)=xf;
    e(3)=yf;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//                          ANALYTIC GRADIENT FUNCTIONS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void OCSolver::NLP::derivatives::analytical::fGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient,Problem &problem){

    //The gradient of the dynamics should be return using the sparsity pattern:

    //grad(f (x, u))=[ df | df ]
    //               [ -- | -- ]
    //               [ dx | du ]

    //dim(grad)=[nStates,nStates+nControls]

    double x =states(0);
    double y =states(1);

    double theta=controls(0);

    double V=2.138;

    //dxdt= V*cos(theta);
    //dydt= V*sin(theta);

    gradient.setZero(2,3);

    gradient<<0,0,-V*sin(theta),    //dxdt w.r.t [x y theta]
              0,0, V*cos(theta);    //dydt w.r.t [x y theta]

}

void OCSolver::NLP::derivatives::analytical::pathGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient,Problem &problem){

    //The gradient of the path constraints should be return using the sparsity pattern:

    //grad(g (x, u))=[ dg | dg ]
    //               [ -- | -- ]
    //               [ dx | du ]

    //dim(grad)=[nPath,nStates+nControls]

    double x =states(0);
    double y =states(1);

     //path(0)=pow(x-0.6,2.0)+pow(y-1,2.0);

    gradient.setZero(2,3);


    gradient<<(2.0*x)-1.2,(2.0*y)-2.0,0,
              (2.0*x)-1.6,(2.0*y)-3.0,0;

}

void OCSolver::NLP::derivatives::analytical::eventGradient(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::MatrixXd &gradient0, Eigen::MatrixXd &gradienttF,Problem &problem){

    //The gradient of the event constraints should be return using the sparsity pattern:

    //grad(e (t0,x0))=[ de  | de  ]
    //                [ --  | --  ]
    //                [ dt0 | dx0 ]


    //grad(e (t0,x0))=[ de  | de  ]
    //                [ --  | --  ]
    //                [ dtF | dxF ]

    //dim(grad)=[nEvents,nStates+1]

    //e(0)=x0;
    //e(1)=y0;
    //e(2)=xf;
    //e(3)=yf;


    gradient0.setZero(4,3);

    gradient0<<0,1.0,0,
               0,0,1.0,
               0,0,0,
               0,0,0;

    gradienttF.setZero(4,3);

    gradienttF<<0,0,0,
                0,0,0,
                0,1.0,0,
                0,0,1.0;

}

