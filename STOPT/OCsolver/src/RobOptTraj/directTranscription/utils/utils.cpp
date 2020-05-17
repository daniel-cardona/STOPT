#include "RobOptTraj/directTranscription.h"


namespace OCSolver {

namespace utils {

double convert_to_original_time(double& tbar,double& t0,double &tf){

    double retval= (tf-t0)/2.0 + (tf-t0)*tbar/2.0;

    return retval;
}

void clip_vector_given_bounds(Eigen::VectorXd &xp, Eigen::VectorXd &xlb, Eigen::VectorXd &xub){

//This function put the variable in x(i) between the xlb and xub

    for(int i=0;i<xp.size();i++){

        if(xp(i)>xub(i)){

            xp(i)=xub(i);

        }

        else if(xp(i)<xlb(i)){

            xp(i)=xlb(i);
        }

    }


}


}//END utils namespace


}//END OCSolver namespace
