#include "RobOptTraj/directTranscription/core/userFunctions.h"

namespace OCSolver
{

namespace collocationMethod
{


    namespace trapezoidal {

          //Cost Function

          double costFunction(Eigen::VectorXd &x,Problem &problem,Alg &algorithm);


          //Constraints functions

          void cnsFunction(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg &algorithm);

          void rightHandSideSparsity(Eigen::VectorXd &x,Eigen::VectorXd &cns,Problem &problem, Alg &algorithm);

    }//END trapezoidalMethod


}//END NLP namespace

}//END OCSolver namespace
