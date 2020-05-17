#include "RobOptTraj/directTranscription/utils/utils.h"

namespace OCSolver {

namespace utils {


    void getDecVarMatrix(Problem &problem, Alg &algorithm,Eigen::VectorXd &x,Eigen::MatrixXd &decVarMatrix,double &t0, double &tF);

    void getVariables(Problem &problem,Eigen::MatrixXd &decVark,Eigen::VectorXd &xk, Eigen::VectorXd &uk,int k);

    void getMidVariablesHSC(Problem &problem,Eigen::MatrixXd& decVark,Eigen::VectorXd &uk_mid, int k);
}


}
