//Main include of STOPT

#include "RobOptTraj/directTranscription.h"


//Read the urdf data

//The .xml data of the robot is in path/to/STOPT/OCSolver/data 

std::string robotFile = "/home/daniel/Documentos/Maestria/RobOptTraj-OpenHRC/OCsolver/data/kukaYoubot.xml";

using std::cout;
using std::endl;

#define infty 100000


//************************************************************************************************
//************************************* MAIN FUNCTION ********************************************
//************************************************************************************************

int main(){

    //Create the main objects of the problem
	
    OCSolver::Problem problem;
    OCSolver::Alg algorithm;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&------------------------ Build KukaYoubot Robot ----------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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

    problem.nStates=problem.robot->getDoF()*2; //16 states
    problem.nControls=problem.robot->getDoF(); //8 controls
    problem.nEvents=problem.robot->getDoF()*4; //32 event constraints
    problem.nPath=0;

    problem.nNodes=60;				//Number of collocation points	

    problem.discretizationMethod="Trapezoidal";	//Collocation method (In this version only Trapezoidal is available)

    problem.setup();				//Setup the problem with the OC data

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&-----------------------Variable bounds-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    //-------------Joint position limits------------

    problem.bounds.states.lower(0)=-1.0;           problem.bounds.states.upper(0)=1.0;
    problem.bounds.states.lower(1)=-1.0;           problem.bounds.states.upper(1)=1.0;
    problem.bounds.states.lower(2)=-0.5;           problem.bounds.states.upper(2)=0.5;

    problem.bounds.states.lower(3)=-2.949;         problem.bounds.states.upper(3)=2.949;
    problem.bounds.states.lower(4)=-1.1344;        problem.bounds.states.upper(4)=1.57;
    problem.bounds.states.lower(5)=-2.61;          problem.bounds.states.upper(5)=2.54;
    problem.bounds.states.lower(6)=-1.78;          problem.bounds.states.upper(6)=1.78;
    problem.bounds.states.lower(7)=-2.91;          problem.bounds.states.upper(7)=2.91;


    //------------Joint velocities limits----------

     problem.bounds.states.lower(8)=-5.0;       problem.bounds.states.upper(8)=5.0;
     problem.bounds.states.lower(9)=-5.0;       problem.bounds.states.upper(9)=5.0;
     problem.bounds.states.lower(10)=-5.0;      problem.bounds.states.upper(10)=5.0;
     problem.bounds.states.lower(11)=-5.0;      problem.bounds.states.upper(11)=5.0;
     problem.bounds.states.lower(12)=-5.0;      problem.bounds.states.upper(12)=5.0;
     problem.bounds.states.lower(13)=-5.0;      problem.bounds.states.upper(13)=5.0;
     problem.bounds.states.lower(14)=-5.0;      problem.bounds.states.upper(14)=5.0;
     problem.bounds.states.lower(15)=-5.0;      problem.bounds.states.upper(15)=5.0;



    //-----------Joint torque limits---------------

     problem.bounds.controls.lower(0)=-10;       problem.bounds.controls.upper(0)=10;
     problem.bounds.controls.lower(1)=-10;       problem.bounds.controls.upper(1)=10;
     problem.bounds.controls.lower(2)=-10;       problem.bounds.controls.upper(2)=10;
     problem.bounds.controls.lower(3)=-10;       problem.bounds.controls.upper(3)=10;
     problem.bounds.controls.lower(4)=-10;       problem.bounds.controls.upper(4)=10;
     problem.bounds.controls.lower(5)=-10;       problem.bounds.controls.upper(5)=10;
     problem.bounds.controls.lower(6)=-10;       problem.bounds.controls.upper(6)=10;
     problem.bounds.controls.lower(7)=-10;       problem.bounds.controls.upper(7)=10;


   //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   //&---------------Boundary and path constraints-------------------&
   //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

     //Initial joint position
     problem.bounds.events.lower(0)=0;
     problem.bounds.events.lower(1)=0;
     problem.bounds.events.lower(2)=0;
     problem.bounds.events.lower(3)=0;
     problem.bounds.events.lower(4)=0;
     problem.bounds.events.lower(5)=0;
     problem.bounds.events.lower(6)=0;
     problem.bounds.events.lower(7)=0;

     //Initial joint velocities
     problem.bounds.events.lower(8)=0;
     problem.bounds.events.lower(9)=0;
     problem.bounds.events.lower(10)=0;
     problem.bounds.events.lower(11)=0;
     problem.bounds.events.lower(12)=0;
     problem.bounds.events.lower(13)=0;
     problem.bounds.events.lower(14)=0;
     problem.bounds.events.lower(15)=0;

     //Final joint position

     problem.bounds.events.lower(16)=0.7;
     problem.bounds.events.lower(17)=-0.5;
     problem.bounds.events.lower(18)=0.5;
     problem.bounds.events.lower(19)=0.0;
     problem.bounds.events.lower(20)=-0.3752;
     problem.bounds.events.lower(21)=1.027;
     problem.bounds.events.lower(22)=-0.11;
     problem.bounds.events.lower(23)=0;


     //Final joint velocities
     problem.bounds.events.lower(24)= 0.0;
     problem.bounds.events.lower(25)= 0.0;
     problem.bounds.events.lower(26)= 0.0;
     problem.bounds.events.lower(27)= 0.0;
     problem.bounds.events.lower(28)= 0.0;
     problem.bounds.events.lower(29)= 0.0;
     problem.bounds.events.lower(30)= 0.0;
     problem.bounds.events.lower(31)= 0.0;

   problem.bounds.events.upper=problem.bounds.events.lower;


    //Initial and final time bounds
   problem.bounds.startTime.lower(0)=0;
   problem.bounds.startTime.upper(0)=0;

   problem.bounds.finalTime.lower(0)=1.0;
   problem.bounds.finalTime.upper(0)=30.0;


  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //&---------------------Set initial guess-------------------------&
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   problem.guess.states.setZero();

   //Set a linear distribution in the joint position states

   problem.guess.states.row(0).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(16));
   problem.guess.states.row(1).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(17));
   problem.guess.states.row(2).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(18));
   problem.guess.states.row(3).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(19));
   problem.guess.states.row(4).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(20));
   problem.guess.states.row(5).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(21));
   problem.guess.states.row(6).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(22));
   problem.guess.states.row(7).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(23));

   problem.guess.controls.setZero();
   problem.guess.time.row(0).setLinSpaced(problem.nNodes,0,10);


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


//************************************************************************************************
//******************************* USER DEFINED FUNCTIONS *****************************************
//************************************************************************************************

double OCSolver::endpoint_cost(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tf, Problem &problem){

    return 0;

}

double OCSolver::integrand_cost(Eigen::VectorXd &states,Eigen::VectorXd &controls,double &tk,Problem &problem){

    //Cost function

    Eigen::VectorXd diff(problem.nStates);

    diff=problem.states_k1-states;

    //Return the value of the cost function!	

    return diff.squaredNorm();


}

//Function to evaluate the dynamics of the system

void OCSolver::dae(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::VectorXd &derivatives,Eigen::VectorXd &path, Problem &problem){

    //Variables of the dynamic system

    hr::VectorXr q(problem.robot->getDoF());
    hr::VectorXr qd(problem.robot->getDoF());
    hr::VectorXr tau(problem.robot->getDoF());
    hr::VectorXr ddq(problem.robot->getDoF());

    hr::core::ForwardDynamics ForwardDynamics(problem.robot);

    q(0)= states(0);     qd(0)=states(8);        tau(0)=controls(0);
    q(1)= states(1);     qd(1)=states(9);        tau(1)=controls(1);
    q(2)= states(2);     qd(2)=states(10);       tau(2)=controls(2);
    q(3)= states(3);     qd(3)=states(11);       tau(3)=controls(3);
    q(4)= states(4);     qd(4)=states(12);       tau(4)=controls(4);
    q(5)= states(5);     qd(5)=states(13);       tau(5)=controls(5);
    q(6)= states(6);     qd(6)=states(14);       tau(6)=controls(6);
    q(7)= states(7);     qd(7)=states(15);       tau(7)=controls(7);

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

void OCSolver::events(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::VectorXd &e , Problem &problem){

    double q1_t0=  initial_states(0);       double q1d_t0=initial_states(8);
    double q2_t0=  initial_states(1);       double q2d_t0=initial_states(9);
    double q3_t0=  initial_states(2);       double q3d_t0=initial_states(10);
    double q4_t0=  initial_states(3);       double q4d_t0=initial_states(11);
    double q5_t0=  initial_states(4);       double q5d_t0=initial_states(12);
    double q6_t0=  initial_states(5);       double q6d_t0=initial_states(13);
    double q7_t0=  initial_states(6);       double q7d_t0=initial_states(14);
    double q8_t0=  initial_states(7);       double q8d_t0=initial_states(15);



    double q1_tf= final_states(0);          double q1d_tf=  final_states(8);
    double q2_tf= final_states(1);          double q2d_tf=  final_states(9);
    double q3_tf= final_states(2);          double q3d_tf=  final_states(10);
    double q4_tf= final_states(3);          double q4d_tf=  final_states(11);
    double q5_tf= final_states(4);          double q5d_tf=  final_states(12);
    double q6_tf= final_states(5);          double q6d_tf=  final_states(13);
    double q7_tf= final_states(6);          double q7d_tf=  final_states(14);
    double q8_tf= final_states(7);          double q8d_tf=  final_states(15);


    //Initial     |  //Final
    //states(q)   |  //states(q)
    e(0)=q1_t0;      e(16)=q1_tf;
    e(1)=q2_t0;      e(17)=q2_tf;
    e(2)=q3_t0;      e(18)=q3_tf;
    e(3)=q4_t0;      e(19)=q4_tf;
    e(4)=q5_t0;      e(20)=q5_tf;
    e(5)=q6_t0;      e(21)=q6_tf;
    e(6)=q7_t0;      e(22)=q7_tf;
    e(7)=q8_t0;      e(23)=q8_tf;


    //Initial     |  //Final
    //states(qd)  |  //states(qd)
    e(8)=q1d_t0;     e(24)=q1d_tf;
    e(9)=q2d_t0;     e(25)=q2d_tf;
    e(10)=q3d_t0;    e(26)=q3d_tf;
    e(11)=q4d_t0;    e(27)=q4d_tf;
    e(12)=q5d_t0;    e(28)=q5d_tf;
    e(13)=q6d_t0;    e(29)=q6d_tf;
    e(14)=q7d_t0;    e(30)=q7d_tf;
    e(15)=q8d_t0;    e(31)=q8d_tf;



}


//************************************************************************************************
//****************************** ANALYTIC GRADIENT FUNCTIONS *************************************
//************************************************************************************************

void OCSolver::NLP::derivatives::analytical::fGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient, Problem &problem){

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

void OCSolver::NLP::derivatives::analytical::pathGradient(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::MatrixXd &gradient, Problem &problem){

    //The gradient of the path constraints should be return using the sparsity pattern:

    //grad(g (x, u))=[ dg | dg ]
    //               [ -- | -- ]
    //               [ dx | du ]

    //dim(grad)=[nPath,nStates+nControls]

}

void OCSolver::NLP::derivatives::analytical::eventGradient(Eigen::VectorXd &initial_states,Eigen::VectorXd &final_states,double &t0,double &tF,Eigen::MatrixXd &gradient0, Eigen::MatrixXd &gradienttF, Problem &problem){

    //The gradient of the event constraints should be return using the sparsity pattern:

    //grad(e (t0,x0))=[ de  | de  ]
    //                [ --  | --  ]
    //                [ dt0 | dx0 ]


    //grad(e (tF,xF))=[ de  | de  ]
    //                [ --  | --  ]
    //                [ dtF | dxF ]

    //dim(grad)=[nEvents,nStates+1]

    gradient0.setZero();
    gradienttF.setZero();

    gradient0.block(0,1,problem.nStates,problem.nStates)=Eigen::MatrixXd::Identity(problem.nStates,problem.nStates);

    gradienttF.block(problem.nStates,1,problem.nStates,problem.nStates)=Eigen::MatrixXd::Identity(problem.nStates,problem.nStates);

}


