#include "RobOptTraj/directTranscription/derivatives/numerical.h"

namespace OCSolver
{

namespace NLP {

namespace derivatives {

namespace analytical {

void obtainAnalyticalGradientInformation(Problem &problem, Alg &algorithm);

void fGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient, Problem &problem);

void pathGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient, Problem &problem);

void eventGradient(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::MatrixXd &gradient0, Eigen::MatrixXd &gradienttF, Problem &problem);

void computeAnalyticalJacobian(Eigen::VectorXd &z, Problem &problem, Alg& algorithm,Eigen::VectorXd &nzValues);

}//END analytical namespace

}//END derivatives namespace

}//END NLP namespace

}//END OCSOlver namespace
