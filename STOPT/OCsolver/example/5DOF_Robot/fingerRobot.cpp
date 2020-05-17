
//Main include of STOPT

#include "RobOptTraj/directTranscription.h"


//Read the urdf data

//The .xml data of the robot is in path/to/STOPT/OCSolver/data 
std::string robotFile = "/home/daniel/Documentos/Maestria/RobOptTraj-OpenHRC/OCsolver/data/5dof.xml";

using std::cout;
using std::endl;

//Set auxiliar functions
void systemDynamics(Eigen::VectorXd &x,Eigen::VectorXd &u, Eigen::VectorXd &xdd);

int main(){

    //Create the main objects of the problem
    OCSolver::Problem problem;
    OCSolver::Alg algorithm;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&-------------- Build the Robot using the rigid body dynamics library OPENHRC ----------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    hr::core::World world = hr::core::World();   		//Create the world
    std::string sFile = robotFile;				//Create a variable for the robot file			
    int robotId = world.getRobotsVector()->size();		//Get an id for the robot identification
    world.loadMultiBodyURDF(sFile,robotId, hr::core::kNAO);	//Upload the robot to the world

    problem.robot = world.getRobot(0);				//Asign the robot object to the problem main information

    //Create useful vector for the states of the robot
	
    Eigen::VectorXd q(problem.robot->getDoF());			
    Eigen::VectorXd qd(problem.robot->getDoF());
    Eigen::VectorXd qdd(problem.robot->getDoF());
    Eigen::VectorXd tau(problem.robot->getDoF());

    //Set the vectors and the states of the robot to zero			

    q.setZero();     problem.robot->setConfiguration(q);
    qd.setZero();    problem.robot->setGeneralizedVelocity(qd);
    qdd.setZero();   problem.robot->setGeneralizedAcceleration(qdd);
    tau.setZero();   problem.robot->setGeneralizedTorques(tau);

    problem.robot->setConfiguration(q);

    //Initialize the forward kinematics of the robot	

    problem.robot->computeForwardKinematics();

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------------OC Problem data --------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    problem.nStates=problem.robot->getDoF()*2; //10 states
    problem.nControls=problem.robot->getDoF(); //5  controls
    problem.nEvents=problem.robot->getDoF()*4; //20 event constraints
    problem.nPath=0;

    problem.nNodes=60;				//Number of collocation points	

    problem.discretizationMethod="Trapezoidal";	//Collocation method (In this version only Trapezoidal is available)

    problem.setup();				//Setup the problem with the OC data

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&-----------------------Variable bounds-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    double q1L  =-pi;     double q1U =pi;
    double q2L  =-pi;     double q2U =pi;
    double q3L  =-pi;     double q3U =pi;
    double q4L  =-pi;     double q4U =pi;
    double q5L  =-pi;     double q5U = pi;

    double q1dL =-10.0;     double q1dU =10.0;
    double q2dL =-10.0;     double q2dU =10.0;
    double q3dL =-10.0;     double q3dU =10.0;
    double q4dL =-10.0;     double q4dU =10.0;
    double q5dL =-10.0;     double q5dU =10.0;


    double T1L = -1.0;     double T1U = 1.0;
    double T2L = -1.0;     double T2U = 1.0;
    double T3L = -1.0;     double T3U = 1.0;
    double T4L = -1.0;     double T4U = 1.0;
    double T5L = -1.0;     double T5U = 1.0;

    //Set the lower and upper bounds of the states
	
    problem.bounds.states.lower(0)=q1L;
    problem.bounds.states.lower(1)=q2L;
    problem.bounds.states.lower(2)=q3L;
    problem.bounds.states.lower(3)=q4L;
    problem.bounds.states.lower(4)=q5L;
    problem.bounds.states.lower(5)=q1dL;
    problem.bounds.states.lower(6)=q2dL;
    problem.bounds.states.lower(7)=q3dL;
    problem.bounds.states.lower(8)=q4dL;
    problem.bounds.states.lower(9)=q5dL;

    problem.bounds.states.upper(0)=q1U;
    problem.bounds.states.upper(1)=q2U;
    problem.bounds.states.upper(2)=q3U;
    problem.bounds.states.upper(3)=q4U;
    problem.bounds.states.upper(4)=q5U;
    problem.bounds.states.upper(5)=q1dU;
    problem.bounds.states.upper(6)=q2dU;
    problem.bounds.states.upper(7)=q3dU;
    problem.bounds.states.upper(8)=q4dU;
    problem.bounds.states.upper(9)=q5dU;

    //Set the lower and upper bounds of the controls

    problem.bounds.controls.lower(0)=T1L;
    problem.bounds.controls.lower(1)=T2L;
    problem.bounds.controls.lower(2)=T3L;
    problem.bounds.controls.lower(3)=T4L;
    problem.bounds.controls.lower(4)=T5L;

    problem.bounds.controls.upper(0)=T1U;
    problem.bounds.controls.upper(1)=T2U;
    problem.bounds.controls.upper(2)=T3U;
    problem.bounds.controls.upper(3)=T4U;
    problem.bounds.controls.upper(4)=T5U;


 //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 //&---------------Boundary and path constraints-------------------&
 //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    //Initial states        //Final states

    double q1_t0  = 0;      double q1_tF  = -pi/5;
    double q2_t0  = 0;      double q2_tF  = pi/4;
    double q3_t0  = 0;      double q3_tF  = pi/3;
    double q4_t0  = 0;      double q4_tF  = pi/4;
    double q5_t0  = 0;      double q5_tF  = pi/8;

    double q1d_t0 = 0;      double q1d_tF = 0;
    double q2d_t0 = 0;      double q2d_tF = 0;
    double q3d_t0 = 0;      double q3d_tF = 0;
    double q4d_t0 = 0;      double q4d_tF = 0;
    double q5d_t0 = 0;      double q5d_tF = 0;

    //Set the bounds of the event constraints
	
    problem.bounds.events.lower(0)=q1_t0;
    problem.bounds.events.lower(1)=q2_t0;
    problem.bounds.events.lower(2)=q3_t0;
    problem.bounds.events.lower(3)=q4_t0;
    problem.bounds.events.lower(4)=q5_t0;
    problem.bounds.events.lower(5)=q1d_t0;
    problem.bounds.events.lower(6)=q2d_t0;
    problem.bounds.events.lower(7)=q3d_t0;
    problem.bounds.events.lower(8)=q4d_t0;
    problem.bounds.events.lower(9)=q5d_t0;

    problem.bounds.events.lower(10)=q1_tF;
    problem.bounds.events.lower(11)=q2_tF;
    problem.bounds.events.lower(12)=q3_tF;
    problem.bounds.events.lower(13)=q4_tF;
    problem.bounds.events.lower(14)=q5_tF;

    problem.bounds.events.lower(15)=q1d_tF;
    problem.bounds.events.lower(16)=q2d_tF;
    problem.bounds.events.lower(17)=q3d_tF;
    problem.bounds.events.lower(18)=q4d_tF;
    problem.bounds.events.lower(19)=q5d_tF;

    problem.bounds.events.upper=problem.bounds.events.lower;

    //Initial and final time boundaries		

    problem.bounds.startTime.lower(0)=0;
    problem.bounds.startTime.upper(0)=0;

    problem.bounds.finalTime.lower(0)=0;
    problem.bounds.finalTime.upper(0)=10;


    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------------Set initial guess-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    int nNodes=problem.nNodes;

    problem.guess.states.row(0).setLinSpaced(nNodes,q1_t0,q1_tF);
    problem.guess.states.row(1).setLinSpaced(nNodes,q2_t0,q2_tF);
    problem.guess.states.row(2).setLinSpaced(nNodes,q3_t0,q3_tF);
    problem.guess.states.row(3).setLinSpaced(nNodes,q4_t0,q4_tF);
    problem.guess.states.row(4).setLinSpaced(nNodes,q5_t0,q5_tF);

    problem.guess.states.row(5).setLinSpaced(nNodes,q1d_t0,q1d_tF);
    problem.guess.states.row(6).setLinSpaced(nNodes,q2d_t0,q2d_tF);
    problem.guess.states.row(7).setLinSpaced(nNodes,q3d_t0,q3d_tF);
    problem.guess.states.row(8).setLinSpaced(nNodes,q4d_t0,q4d_tF);
    problem.guess.states.row(9).setLinSpaced(nNodes,q5d_t0,q5d_tF);

    problem.guess.controls.row(0).setZero();
    problem.guess.controls.row(1).setZero();
    problem.guess.controls.row(2).setZero();
    problem.guess.controls.row(3).setZero();
    problem.guess.controls.row(4).setZero();

    problem.guess.time.row(0).setLinSpaced(nNodes,0,10);

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

    algorithm.scaling               ="none";
    algorithm.derivatives           ="analytical";  //Set the algorithm for compute the first order information! (Numerical OR Analytical)
    algorithm.defect_scaling        ="state-based";
    algorithm.derivativeChecker     =false;
    algorithm.sparsityOptimization  =true;


    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&---------------------CALL THE SOLVER!--------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    OCSolver::core::robTrajSol(problem,algorithm);


    return 0;

}


double OCSolver::endpoint_cost(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tf, Problem &problem){

    return 0;

}

double OCSolver::integrand_cost(Eigen::VectorXd &states,Eigen::VectorXd &controls,double &tk,Problem &problem){

   //Cost function

    double u1=controls(0);  //Tau 1
    double u2=controls(1);  //Tau 2
    double u3=controls(2);  //Tau 3
    double u4=controls(3);  //Tau 4
    double u5=controls(4);  //Tau 5


    double L=pow(u1,2.0)+ pow(u2,2.0) + pow(u3,2.0) + pow(u4,2.0) + pow(u5,2.0);

    //Return the value of the cost function!	

    return L;

}

//Function to evaluate the dynamics of the system
 
void OCSolver::dae(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::VectorXd &derivatives,Eigen::VectorXd &path, Problem &problem){

    //Variables of the dynamic system

        hr::VectorXr q(problem.robot->getDoF());
        hr::VectorXr qd(problem.robot->getDoF());
        hr::VectorXr tau(problem.robot->getDoF());
        hr::VectorXr ddq(problem.robot->getDoF());

        hr::core::ForwardDynamics ForwardDynamics(problem.robot);

        q(0)= states(0);     qd(0)=states(5);       tau(0)=controls(0);
        q(1)= states(1);     qd(1)=states(6);       tau(1)=controls(1);
        q(2)= states(2);     qd(2)=states(7);       tau(2)=controls(2);
        q(3)= states(3);     qd(3)=states(8);       tau(3)=controls(3);
        q(4)= states(4);     qd(4)=states(9);       tau(4)=controls(4);


        problem.robot->setConfiguration(q);
        problem.robot->setGeneralizedVelocity(qd);
        problem.robot->setGeneralizedTorques(tau);

        bool computePartialDerivatives=false;

        ForwardDynamics.computeForwardDynamics(computePartialDerivatives);

        derivatives.setZero(problem.nStates);

        ddq=ForwardDynamics.getJointAcceleration();

	//Return the vector \dot{x}

        derivatives<<qd,ddq;


}

//Function to evaluate the event constraints

void OCSolver::events(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::VectorXd &e, Problem &problem){

               double q1_t0  = initial_states(0);
               double q2_t0  = initial_states(1);
               double q3_t0  = initial_states(2);
               double q4_t0  = initial_states(3);
               double q5_t0  = initial_states(4);

               double q1d_t0 = initial_states(5);
               double q2d_t0 = initial_states(6);
               double q3d_t0 = initial_states(7);
               double q4d_t0 = initial_states(8);
               double q5d_t0 = initial_states(9);


               double l1=0.01;
               double l2=0.05;
               double l3=0.03;
               double l4=0.02;

               double q1_tF  = final_states(0);
               double q2_tF  = final_states(1);
               double q3_tF  = final_states(2);
               double q4_tF  = final_states(3);
               double q5_tF  = final_states(4);

               double q1d_tF = final_states(5);
               double q2d_tF = final_states(6);
               double q3d_tF = final_states(7);
               double q4d_tF = final_states(8);
               double q5d_tF = final_states(9);


               e.setZero(20); //20 event constraints

               e(0) = q1_t0;
               e(1) = q2_t0;
               e(2) = q3_t0;
               e(3) = q4_t0;
               e(4) = q5_t0;
               e(5) = q1d_t0;
               e(6) = q2d_t0;
               e(7) = q3d_t0;
               e(8) = q4d_t0;
               e(9) = q5d_t0;

               e(10) = q1_tF;
               e(11) = q2_tF;
               e(12) = q3_tF;
               e(13) = q4_tF;
               e(14) = q5_tF;

               e(15) = q1d_tF;
               e(16) = q2d_tF;
               e(17) = q3d_tF;
               e(18) = q4d_tF;
               e(19)=  q5d_tF;



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

        hr::VectorXr q(problem.robot->getDoF());
        hr::VectorXr qd(problem.robot->getDoF());
        hr::VectorXr tau(problem.robot->getDoF());

        q=states.head(problem.robot->getDoF());
        qd=states.tail(problem.robot->getDoF());
        tau=controls;

        problem.robot->setConfiguration(q);
        problem.robot->setGeneralizedVelocity(qd);
        problem.robot->setGeneralizedTorques(tau);

        //! Forward dynamics object
        hr::core::ForwardDynamics ForwardDynamics(problem.robot);

        bool computePartialDerivatives=true;

        ForwardDynamics.computeForwardDynamics(computePartialDerivatives);

        ForwardDynamics.computeInverseInertiaMatrix();

        gradient.setZero();

        gradient.block(0,problem.robot->getDoF(),problem.robot->getDoF(),problem.robot->getDoF())=Eigen::MatrixXd::Identity(problem.robot->getDoF(),problem.robot->getDoF());

        gradient.block(problem.robot->getDoF(),0,problem.robot->getDoF(),problem.nStates)=ForwardDynamics.getDiffJointAcceleration();

        gradient.block(problem.robot->getDoF(),problem.nStates,problem.robot->getDoF(),problem.robot->getDoF())=ForwardDynamics.getInverseInertiaMatrix();

}

void OCSolver::NLP::derivatives::analytical::pathGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient,Problem &problem){

    //The gradient of the path constraints should be return using the sparsity pattern:

    //grad(g (x, u))=[ dg | dg ]
    //               [ -- | -- ]
    //               [ dx | du ]

    //dim(grad)=[nPath,nStates+nControls]

}

void OCSolver::NLP::derivatives::analytical::eventGradient(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::MatrixXd &gradient0, Eigen::MatrixXd &gradienttF,Problem &problem){

    //The gradient of the event constraints should be return using the sparsity pattern:

    //grad(e (t0,x0))=[ de  | de  ]
    //                [ --  | --  ]
    //                [ dt0 | dx0 ]


    //grad(e (t0,x0))=[ de  | de  ]
    //                [ --  | --  ]
    //                [ dtF | dxF

    //dim(grad)=[nEvents,nStates+1]

    //e(0)=x0;
    //e(1)=y0;
    //e(2)=xf;
    //e(3)=yf;

        gradient0.setZero();
        gradienttF.setZero();

        gradient0.block(0,1,problem.nStates,problem.nStates)=Eigen::MatrixXd::Identity(problem.nStates,problem.nStates);

        gradienttF.block(problem.nStates,1,problem.nStates,problem.nStates)=Eigen::MatrixXd::Identity(problem.nStates,problem.nStates);

}
