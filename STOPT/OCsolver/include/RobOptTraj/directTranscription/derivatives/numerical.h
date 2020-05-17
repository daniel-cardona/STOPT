#include "RobOptTraj/directTranscription/nlpSolver/IPOPT_Interface.h"

namespace OCSolver
{

namespace NLP
{

namespace derivatives
{

namespace numerical {


//! Function that obtains the gradient of of the cost function w.r.t. the decision variables
//!
//! This function is stored in: src/derivatives/numericalDerivatives.cpp
//!
/*! \param $fun$               Pointer to the function to evaluate the cost
 *  \param $x$                  Vector with the decision variables
 *  \param $Prob$              Object with the problem information
 *  \param $Alg$               Object with the algorithm options
 *
 *  \return $grad$             Gradient of the cost function w.r.t. to the decision variables
 */

void scalarGradient(double(*fun)(Eigen::VectorXd &x,Problem &problem,Alg &algorithm),Eigen::VectorXd &x,Problem &problem,Alg &algorithm,Eigen::VectorXd &grad);

//! Function that obtains one column of the jacobian of the constraints
//!
/*! \param $fun$               Pointer to the function to evaluate the constraints
 *  \param $x$                 Vector with the decision variables
 *  \param $jCol$              Index of the column of the jacobian that is going to be calculated
 *  \param $Prob$              Object with the problem information
 *  \param $Alg$               Object with the algorithm options
 *
 *  \return $grad$             Gradient of the cost function w.r.t. to the decision variables
 */
void jacobianColumn(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x, Problem &problem, Alg &algorithm,int jCol,Eigen::VectorXd &jacColumn);


//! Function that obtains one column of the jacobian of the constraints
//!
/*! \param $fun$               Pointer to the function to evaluate the constraints
 *  \param $x$                 Vector with the decision variables
 *  \param $jCol$              Index of the column of the jacobian that is going to be calculated
 *  \param $Prob$              Object with the problem information
 *  \param $Alg$               Object with the algorithm options
 *
 *  \return $grad$             Gradient of the cost function w.r.t. to the decision variables
 */

void computeNumericalJacobian(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm, Eigen::VectorXd &nzValues);



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEPRECATED FUNCTIONS

void computeNumericalJacobianv2(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm, Eigen::VectorXd &nzValues);


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION 3 NUMERICAL DERIVATIVES FUNCTIONS %%%%%%%%%%%%%%%%%


void dxJacColumn(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem);

void pathJacColumn(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem);

void eventJacColumn_t0(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem);

void eventJacColumn_tf(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem);

}//END numerical namespace




}//END derivatives namespace

}//END NLP namespace

}//END OCSolver namespace
