#include "RobOptTraj/directTranscription/derivatives/analytical.h"

namespace OCSolver {

namespace utils {

    double convert_to_original_time(double& tbar,double& t0,double &tf);

    void clip_vector_given_bounds(Eigen::VectorXd &xp, Eigen::VectorXd &xlb, Eigen::VectorXd &xub);

}


}
