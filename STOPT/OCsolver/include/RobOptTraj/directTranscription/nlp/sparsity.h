#include "RobOptTraj/directTranscription/nlp/scaling.h"


namespace OCSolver
{

namespace NLP
{


namespace sparsity {


void detectFullJacobianSparsity(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm);

//! Function that obtain the index groups used for the aproximation of the numerical jacobian
//!
//! This function is stored in: src/derivatives/numericalDerivatives.cpp
//!
/*! \param $Prob$              Object with the problem information
 *
 *  \return void
*/

void getIndexGroups(Problem &problem);

//! Function that obtain the index groups used for the aproximation of the numerical jacobian
//!
//! This function is stored in: src/derivatives/numericalDerivatives.cpp
//!
/*! \param $Prob$              Object with the problem information
 *
 *  \return void
*/
void detectDerivativeMatrixSparsity(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm);


//! Function that obtain the index groups used for the aproximation of the numerical jacobian
//!
//! This function is stored in: src/derivatives/numericalDerivatives.cpp
//!
/*! \param $Prob$              Object with the problem information
 *
 *  \return void
*/
void getIndexGroupsDerivativeMatrix(Problem &problem);


//! Function that obtain the index groups used for the aproximation of the numerical jacobian
//!
//! This function is stored in: src/derivatives/numericalDerivatives.cpp
//!
/*! \param $Prob$              Object with the problem information
 *
 *  \return void
*/
void getPerturbationMatrix(Problem &problem,Alg &algorithm);


//! Function that obtain the index groups used for the aproximation of the numerical jacobian
//!
//! This function is stored in: src/derivatives/numericalDerivatives.cpp
//!
/*! \param $Prob$              Object with the problem information
 *
 *  \return void
*/
void generateSparsityTemplates(Problem &problem,Alg &algorithm);

//%%%%%%%%%%%%%%%%%%% VERSION 3 SPARSITY FUNCTIONS %%%%%%%%%%%%%%%%%%%

void detectSparsity(Problem &problem);

void detectSparsity_dx(Problem &problem,Eigen::MatrixXd &sparsityTemplate);

void detectSparsity_path(Problem &problem, Eigen::MatrixXd &sparsityTemplate_path);

void detectSparsity_events(Problem &problem,Eigen::MatrixXd &sparsityTemplate_t0,Eigen::MatrixXd &sparsityTemplate_tf);

}


}//END NLP namespace
}//END OCSolver namespace
