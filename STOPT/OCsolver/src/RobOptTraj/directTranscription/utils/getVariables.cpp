#include "RobOptTraj/directTranscription.h"


namespace OCSolver {

namespace utils {


void getDecVarMatrix(Problem &problem, Alg &algorithm,Eigen::VectorXd &x,Eigen::MatrixXd &decVarMatrix,double &t0, double &tF){

    //Recover and return the initial and final time

    t0=x(0)/problem.scale.time;
    tF=x(1)/problem.scale.time;


    //-----------------TRAPEZOIDAL DEFECT CONSTRAINTS------------------------------

    if(problem.discretizationMethod=="Trapezoidal"){

        //Now recover the matrix with the following structure:
        //[xk xk+1 xk+2 ...xM
        // uk uk+1 uk+2 ...uM]

        Eigen::Map<Eigen::MatrixXd> decVar(x.tail(x.size()-2).data(),problem.nStates+problem.nControls,problem.nNodes);

        //Return this!
        decVarMatrix=decVar;

    }
    //-----------------HERMITE-SIMPSON DEFECT CONSTRAINTS---------------------------

    if(problem.discretizationMethod=="Hermite-Simpson"){

        //Now recover the matrix with the following structure:
        //[xk xk+1 xk+2 ...xM
        // uk uk+1 uk+2 ...uM
        // ük ük+1 ük+2 ... 0]

        Eigen::VectorXd xHs(x.size()+problem.nControls);
        xHs<<x,Eigen::VectorXd::Zero(problem.nControls);

        Eigen::Map<Eigen::MatrixXd> decVar(xHs.tail(xHs.size()-2).data(),problem.nStates+problem.nControls*2,problem.nNodes);

        //Return this!
        decVarMatrix=decVar;

    }


}

void getVariables(Problem &problem,Eigen::MatrixXd &decVark,Eigen::VectorXd &xk, Eigen::VectorXd &uk,int k){

    //Obtain the states and the controls in k and unscaled it

    xk=decVark.col(k).head(problem.nStates).array()/problem.scale.states.array();

    uk=decVark.col(k).segment(problem.nStates,problem.nControls).array()/problem.scale.controls.array();


}


void getMidVariablesHSC(Problem &problem,Eigen::MatrixXd& decVark,Eigen::VectorXd &uk_mid, int k){

    uk_mid=decVark.col(k).tail(problem.nControls).array()/problem.scale.controls.array();

}

}//END utils namespace


}//END OCSolver namespace
