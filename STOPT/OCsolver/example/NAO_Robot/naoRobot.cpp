
//Main include of STOPT
#include "RobOptTraj/directTranscription.h"


//Read the urdf data

//The .xml data of the robot is in path/to/STOPT/OCSolver/data 
std::string robotFile = "/home/daniel/Documentos/Maestria/RobOptTraj-OpenHRC/OCsolver/data/NAOURDF_inertial.xml";

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

    problem.nStates=problem.robot->getDoF()*2; //48 states
    problem.nControls=problem.robot->getDoF(); //24 controls
    problem.nEvents=problem.robot->getDoF()*4; //96 event constraints
    problem.nPath=0;

    problem.nNodes=90;				//Number of collocation points	

    problem.discretizationMethod="Trapezoidal";	//Collocation method (In this version only Trapezoidal is available)

    problem.setup();				//Setup the problem with the OC data


    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //&-----------------------Variable bounds-------------------------&
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    //-------------Joint position limits------------

     problem.bounds.states.lower(0)=-0.39270;       problem.bounds.states.upper(0)=0.76969;
     problem.bounds.states.lower(1)=-1.18857;       problem.bounds.states.upper(1)=0.92328;
     problem.bounds.states.lower(2)=-0.08552;       problem.bounds.states.upper(2)=2.10138;
     problem.bounds.states.lower(3)=-1.53065;       problem.bounds.states.upper(3)=0.48346;
     problem.bounds.states.lower(4)=-0.37350;       problem.bounds.states.upper(4)=0.78540;
     problem.bounds.states.lower(5)=-1.13970;       problem.bounds.states.upper(5)=0.74002;

     problem.bounds.states.lower(6)=-2.08218;       problem.bounds.states.upper(6)=2.07694;
     problem.bounds.states.lower(7)=-0.19722;       problem.bounds.states.upper(7)=1.32645;
     problem.bounds.states.lower(8)=-2.08392;       problem.bounds.states.upper(8)=2.08218;
     //problem.bounds.states.lower(9)=-1.53589;       problem.bounds.states.upper(9)=-0.03840;
     problem.bounds.states.lower(9)=-1.53589;       problem.bounds.states.upper(9)=0.03840;
     problem.bounds.states.lower(10)=-2.08392;      problem.bounds.states.upper(10)=2.08218;

     problem.bounds.states.lower(11)=-2.07694;      problem.bounds.states.upper(11)=2.08392;
     problem.bounds.states.lower(12)=-0.66148;      problem.bounds.states.upper(12)=0.50964;

     problem.bounds.states.lower(13)=-2.08218;      problem.bounds.states.upper(13)=2.08567;   
     problem.bounds.states.lower(14)=-1.32645;      problem.bounds.states.upper(14)=0.20420;
     problem.bounds.states.lower(15)=-2.08392;      problem.bounds.states.upper(15)=2.08218;
     //problem.bounds.states.lower(16)=0.03665;       problem.bounds.states.upper(16)=1.53589;
     problem.bounds.states.lower(16)=-0.03665;       problem.bounds.states.upper(16)=1.53589;
     problem.bounds.states.lower(17)=-2.08392;      problem.bounds.states.upper(17)=2.08218;

     problem.bounds.states.lower(18)=-1.14494;      problem.bounds.states.upper(18)=0.74002;
     problem.bounds.states.lower(19)=-0.78714;      problem.bounds.states.upper(19)=0.37874;
     problem.bounds.states.lower(20)=-1.53065;      problem.bounds.states.upper(20)=0.47822;
     problem.bounds.states.lower(21)=-0.08901;      problem.bounds.states.upper(21)=2.10836;
     problem.bounds.states.lower(22)=-1.18682;      problem.bounds.states.upper(22)=0.93201;
     problem.bounds.states.lower(23)=-0.76794;      problem.bounds.states.upper(23)=0.39444;



    //------------Joint velocities limits----------

     problem.bounds.states.lower(24)=-5.0;      problem.bounds.states.upper(24)=5.0;
     problem.bounds.states.lower(25)=-5.0;      problem.bounds.states.upper(25)=5.0;
     problem.bounds.states.lower(26)=-5.0;      problem.bounds.states.upper(26)=5.0;
     problem.bounds.states.lower(27)=-5.0;      problem.bounds.states.upper(27)=5.0;
     problem.bounds.states.lower(28)=-5.0;      problem.bounds.states.upper(28)=5.0;
     problem.bounds.states.lower(29)=-5.0;      problem.bounds.states.upper(29)=5.0;
     problem.bounds.states.lower(30)=-5.0;      problem.bounds.states.upper(30)=5.0;
     problem.bounds.states.lower(31)=-5.0;      problem.bounds.states.upper(31)=5.0;
     problem.bounds.states.lower(32)=-5.0;      problem.bounds.states.upper(32)=5.0;
     problem.bounds.states.lower(33)=-5.0;      problem.bounds.states.upper(33)=5.0;
     problem.bounds.states.lower(34)=-5.0;      problem.bounds.states.upper(34)=5.0;
     problem.bounds.states.lower(35)=-5.0;      problem.bounds.states.upper(35)=5.0;
     problem.bounds.states.lower(36)=-5.0;      problem.bounds.states.upper(36)=5.0;
     problem.bounds.states.lower(37)=-5.0;      problem.bounds.states.upper(37)=5.0;
     problem.bounds.states.lower(38)=-5.0;      problem.bounds.states.upper(38)=5.0;
     problem.bounds.states.lower(39)=-5.0;      problem.bounds.states.upper(39)=5.0;
     problem.bounds.states.lower(40)=-5.0;      problem.bounds.states.upper(40)=5.0;
     problem.bounds.states.lower(41)=-5.0;      problem.bounds.states.upper(41)=5.0;
     problem.bounds.states.lower(42)=-5.0;      problem.bounds.states.upper(42)=5.0;
     problem.bounds.states.lower(43)=-5.0;      problem.bounds.states.upper(43)=5.0;
     problem.bounds.states.lower(44)=-5.0;      problem.bounds.states.upper(44)=5.0;
     problem.bounds.states.lower(45)=-5.0;      problem.bounds.states.upper(45)=5.0;
     problem.bounds.states.lower(46)=-5.0;      problem.bounds.states.upper(46)=5.0;
     problem.bounds.states.lower(47)=-5.0;      problem.bounds.states.upper(47)=5.0;

    //-----------Joint torque limits---------------

     problem.bounds.controls.lower(0)=-10;       problem.bounds.controls.upper(0)=10;
     problem.bounds.controls.lower(1)=-10;       problem.bounds.controls.upper(1)=10;
     problem.bounds.controls.lower(2)=-10;       problem.bounds.controls.upper(2)=10;
     problem.bounds.controls.lower(3)=-10;       problem.bounds.controls.upper(3)=10;
     problem.bounds.controls.lower(4)=-10;       problem.bounds.controls.upper(4)=10;
     problem.bounds.controls.lower(5)=-10;       problem.bounds.controls.upper(5)=10;
     problem.bounds.controls.lower(6)=-10;       problem.bounds.controls.upper(6)=10;
     problem.bounds.controls.lower(7)=-10;       problem.bounds.controls.upper(7)=10;
     problem.bounds.controls.lower(8)=-10;       problem.bounds.controls.upper(8)=10;
     problem.bounds.controls.lower(9)=-10;       problem.bounds.controls.upper(9)=10;
     problem.bounds.controls.lower(10)=-10;      problem.bounds.controls.upper(10)=10;
     problem.bounds.controls.lower(11)=-10;      problem.bounds.controls.upper(11)=10;
     problem.bounds.controls.lower(12)=-10;      problem.bounds.controls.upper(12)=10;
     problem.bounds.controls.lower(13)=-10;      problem.bounds.controls.upper(13)=10;
     problem.bounds.controls.lower(14)=-10;      problem.bounds.controls.upper(14)=10;
     problem.bounds.controls.lower(15)=-10;      problem.bounds.controls.upper(15)=10;
     problem.bounds.controls.lower(16)=-10;      problem.bounds.controls.upper(16)=10;
     problem.bounds.controls.lower(17)=-10;      problem.bounds.controls.upper(17)=10;
     problem.bounds.controls.lower(18)=-10;      problem.bounds.controls.upper(18)=10;
     problem.bounds.controls.lower(19)=-10;      problem.bounds.controls.upper(19)=10;
     problem.bounds.controls.lower(20)=-10;      problem.bounds.controls.upper(20)=10;
     problem.bounds.controls.lower(21)=-10;      problem.bounds.controls.upper(21)=10;
     problem.bounds.controls.lower(22)=-10;      problem.bounds.controls.upper(22)=10;
     problem.bounds.controls.lower(23)=-10;      problem.bounds.controls.upper(23)=10;


   //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   //&---------------Boundary and path constraints-------------------&
   //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

     Eigen::VectorXd initial_configuration(problem.robot->getDoF());
     Eigen::VectorXd end_configuration(problem.robot->getDoF());

     initial_configuration<< 0.06109,
			     -1.03847,
			     2.09090,
			     -1.03847,
			     0.41539,
			     -0.00175,
			     0.18850,
			     0.10297,
			     -0.00349,
			     -1.46084,
			     0.00000,
			     -1.18333,
			     0.00000,
			     1.04545,
			     -1.32470,
			     1.56731,
			     0.07854,
			     0.00000,
			     -0.00175,
			     -0.54629,
			     -0.10647,
			     0,
			     -0.00873,
			     0.39444;

     initial_configuration.setZero();

     //Pose 1
     end_configuration<<   -0.0070,
                           -0.8639,
                            1.9967,
                           -1.1345,
                            0.0052,
                            0.0000,
                            0.0018,
                            0.8919,
                           -2.0822,
                           -0.0436,
                           -0.0000,
                           -0.9861,
                           -0.0000,
                            0.0018,
                           -0.9756,
                            2.0734,
                            0.0384,
                           -0.0000,
                            0.0000,
                           -0.0018,
                           -1.1345,
                            1.9967,
                           -0.8639,
                            0.0035;

     /*/Pose 2

     end_configuration<<      0.0611,
                             -1.0385,
                              2.0909,
                             -1.0385,
                              0.4154,
                             -0.0018,
                              0.1885,
                              0.1030,
                             -0.0035,
                             -1.4608,
                             -0.0000,
                             -1.1833,
                             -0.0000,
                              1.0454,
                             -1.3247,
                              1.5673,
                              0.0785,
                              0.0000,
                             -0.0018,
                             -0.5463,
                             -0.1065,
                             -0.0000,
                             -0.0087,
                              0.3944;*/

     //Initial position joints                  			//Final joint positions
     problem.bounds.events.lower(0)=initial_configuration(0);          problem.bounds.events.lower(48)=end_configuration(0);
     problem.bounds.events.lower(1)=initial_configuration(1);          problem.bounds.events.lower(49)=end_configuration(1);
     problem.bounds.events.lower(2)=initial_configuration(2);          problem.bounds.events.lower(50)=end_configuration(2);
     problem.bounds.events.lower(3)=initial_configuration(3);          problem.bounds.events.lower(51)=end_configuration(3);
     problem.bounds.events.lower(4)=initial_configuration(4);          problem.bounds.events.lower(52)=end_configuration(4);
     problem.bounds.events.lower(5)=initial_configuration(5);          problem.bounds.events.lower(53)=end_configuration(5);

     problem.bounds.events.lower(6)=initial_configuration(6);          problem.bounds.events.lower(54)= end_configuration(6);
     problem.bounds.events.lower(7)=initial_configuration(7);          problem.bounds.events.lower(55)= end_configuration(7);
     problem.bounds.events.lower(8)=initial_configuration(8);          problem.bounds.events.lower(56)=end_configuration(8);
     problem.bounds.events.lower(9)=initial_configuration(9);          problem.bounds.events.lower(57)=end_configuration(9);
     problem.bounds.events.lower(10)=initial_configuration(10);        problem.bounds.events.lower(58)=end_configuration(10);

     problem.bounds.events.lower(11)=initial_configuration(11);        problem.bounds.events.lower(59)= end_configuration(11);
     problem.bounds.events.lower(12)=initial_configuration(12);        problem.bounds.events.lower(60)= end_configuration(12);

     problem.bounds.events.lower(13)=initial_configuration(13);        problem.bounds.events.lower(61)= end_configuration(13);
     problem.bounds.events.lower(14)=initial_configuration(14);        problem.bounds.events.lower(62)=end_configuration(14);
     problem.bounds.events.lower(15)=initial_configuration(15);        problem.bounds.events.lower(63)= end_configuration(15);
     problem.bounds.events.lower(16)=initial_configuration(16);        problem.bounds.events.lower(64)= end_configuration(16);
     problem.bounds.events.lower(17)=initial_configuration(17);        problem.bounds.events.lower(65)= end_configuration(17);

     problem.bounds.events.lower(18)=initial_configuration(18);        problem.bounds.events.lower(66)=end_configuration(18);
     problem.bounds.events.lower(19)=initial_configuration(19);        problem.bounds.events.lower(67)= end_configuration(19);
     problem.bounds.events.lower(20)=initial_configuration(20);        problem.bounds.events.lower(68)=end_configuration(20);
     problem.bounds.events.lower(21)=initial_configuration(21);        problem.bounds.events.lower(69)= end_configuration(21);
     problem.bounds.events.lower(22)=initial_configuration(22);        problem.bounds.events.lower(70)=end_configuration(22);
     problem.bounds.events.lower(23)=initial_configuration(23);        problem.bounds.events.lower(71)=end_configuration(23);


     //Initial joint velocities                 //Final joint velocities
     problem.bounds.events.lower(24)=0;         problem.bounds.events.lower(72)=0;
     problem.bounds.events.lower(25)=0;         problem.bounds.events.lower(73)=0;
     problem.bounds.events.lower(26)=0;         problem.bounds.events.lower(74)=0;
     problem.bounds.events.lower(27)=0;         problem.bounds.events.lower(75)=0;
     problem.bounds.events.lower(28)=0;         problem.bounds.events.lower(76)=0;
     problem.bounds.events.lower(29)=0;         problem.bounds.events.lower(77)=0;
     problem.bounds.events.lower(30)=0;         problem.bounds.events.lower(78)=0;
     problem.bounds.events.lower(31)=0;         problem.bounds.events.lower(79)=0;
     problem.bounds.events.lower(32)=0;         problem.bounds.events.lower(80)=0;
     problem.bounds.events.lower(33)=0;         problem.bounds.events.lower(81)=0;
     problem.bounds.events.lower(34)=0;         problem.bounds.events.lower(82)=0;
     problem.bounds.events.lower(35)=0;         problem.bounds.events.lower(83)=0;
     problem.bounds.events.lower(36)=0;         problem.bounds.events.lower(84)=0;
     problem.bounds.events.lower(37)=0;         problem.bounds.events.lower(85)=0;
     problem.bounds.events.lower(38)=0;         problem.bounds.events.lower(86)=0;
     problem.bounds.events.lower(39)=0;         problem.bounds.events.lower(87)=0;
     problem.bounds.events.lower(40)=0;         problem.bounds.events.lower(88)=0;
     problem.bounds.events.lower(41)=0;         problem.bounds.events.lower(89)=0;
     problem.bounds.events.lower(42)=0;         problem.bounds.events.lower(90)=0;
     problem.bounds.events.lower(43)=0;         problem.bounds.events.lower(91)=0;
     problem.bounds.events.lower(44)=0;         problem.bounds.events.lower(92)=0;
     problem.bounds.events.lower(45)=0;         problem.bounds.events.lower(93)=0;
     problem.bounds.events.lower(46)=0;         problem.bounds.events.lower(94)=0;
     problem.bounds.events.lower(47)=0;         problem.bounds.events.lower(95)=0;


   problem.bounds.events.upper=problem.bounds.events.lower;

   //Initial and final time boundaries

   problem.bounds.startTime.lower(0)=0;
   problem.bounds.startTime.upper(0)=0;

   problem.bounds.finalTime.lower(0)=1.0;
   problem.bounds.finalTime.upper(0)=20.0;


  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //&---------------------Set initial guess-------------------------&
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   problem.guess.states.setZero();

   //Set a linear distribution in the joint position states

   problem.guess.states.row(0).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(48));
   problem.guess.states.row(1).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(49));
   problem.guess.states.row(2).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(50));
   problem.guess.states.row(3).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(51));
   problem.guess.states.row(4).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(52));
   problem.guess.states.row(5).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(53));
   problem.guess.states.row(6).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(54));
   problem.guess.states.row(7).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(55));
   problem.guess.states.row(8).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(56));
   problem.guess.states.row(9).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(57));
   problem.guess.states.row(10).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(58));
   problem.guess.states.row(11).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(59));
   problem.guess.states.row(12).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(60));
   problem.guess.states.row(13).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(61));
   problem.guess.states.row(14).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(62));
   problem.guess.states.row(15).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(63));
   problem.guess.states.row(16).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(64));
   problem.guess.states.row(17).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(65));
   problem.guess.states.row(18).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(66));
   problem.guess.states.row(19).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(67));
   problem.guess.states.row(20).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(68));
   problem.guess.states.row(21).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(69));
   problem.guess.states.row(22).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(70));
   problem.guess.states.row(23).setLinSpaced(problem.nNodes,0, problem.bounds.events.lower(71));

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

void OCSolver::dae(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::VectorXd &derivatives,Eigen::VectorXd &path, Problem &problem){

    //Variables of the dynamic system

    hr::VectorXr q(problem.robot->getDoF());
    hr::VectorXr qd(problem.robot->getDoF());
    hr::VectorXr tau(problem.robot->getDoF());
    hr::VectorXr ddq(problem.robot->getDoF());

    hr::core::ForwardDynamics ForwardDynamics(problem.robot);

    q(0)= states(0);     qd(0)=states(24);       tau(0)=controls(0);
    q(1)= states(1);     qd(1)=states(25);       tau(1)=controls(1);
    q(2)= states(2);     qd(2)=states(26);       tau(2)=controls(2);
    q(3)= states(3);     qd(3)=states(27);       tau(3)=controls(3);
    q(4)= states(4);     qd(4)=states(28);       tau(4)=controls(4);
    q(5)= states(5);     qd(5)=states(29);       tau(5)=controls(5);
    q(6)= states(6);     qd(6)=states(30);       tau(6)=controls(6);
    q(7)= states(7);     qd(7)=states(31);       tau(7)=controls(7);
    q(8)= states(8);     qd(8)=states(32);       tau(8)=controls(8);
    q(9)= states(9);     qd(9)=states(33);       tau(9)=controls(9);
    q(10)=states(10);    qd(10)=states(34);      tau(10)=controls(10);
    q(11)=states(11);    qd(11)=states(35);      tau(11)=controls(11);
    q(12)=states(12);    qd(12)=states(36);      tau(12)=controls(12);
    q(13)=states(13);    qd(13)=states(37);      tau(13)=controls(13);
    q(14)=states(14);    qd(14)=states(38);      tau(14)=controls(14);
    q(15)=states(15);    qd(15)=states(39);      tau(15)=controls(15);
    q(16)=states(16);    qd(16)=states(40);      tau(16)=controls(16);
    q(17)=states(17);    qd(17)=states(41);      tau(17)=controls(17);
    q(18)=states(18);    qd(18)=states(42);      tau(18)=controls(18);
    q(19)=states(19);    qd(19)=states(43);      tau(19)=controls(19);
    q(20)=states(20);    qd(20)=states(44);      tau(20)=controls(20);
    q(21)=states(21);    qd(21)=states(45);      tau(21)=controls(21);
    q(22)=states(22);    qd(22)=states(46);      tau(22)=controls(22);
    q(23)=states(23);    qd(23)=states(47);      tau(23)=controls(23);

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

    double q1_t0=  initial_states(0);       double q1d_t0=initial_states(24);
    double q2_t0=  initial_states(1);       double q2d_t0=initial_states(25);
    double q3_t0=  initial_states(2);       double q3d_t0=initial_states(26);
    double q4_t0=  initial_states(3);       double q4d_t0=initial_states(27);
    double q5_t0=  initial_states(4);       double q5d_t0=initial_states(28);
    double q6_t0=  initial_states(5);       double q6d_t0=initial_states(29);
    double q7_t0=  initial_states(6);       double q7d_t0=initial_states(30);
    double q8_t0=  initial_states(7);       double q8d_t0=initial_states(31);
    double q9_t0=  initial_states(8);       double q9d_t0=initial_states(32);
    double q10_t0= initial_states(9);       double q10d_t0=initial_states(33);
    double q11_t0= initial_states(10);      double q11d_t0=initial_states(34);
    double q12_t0= initial_states(11);      double q12d_t0=initial_states(35);
    double q13_t0= initial_states(12);      double q13d_t0=initial_states(36);
    double q14_t0= initial_states(13);      double q14d_t0=initial_states(37);
    double q15_t0= initial_states(14);      double q15d_t0=initial_states(38);
    double q16_t0= initial_states(15);      double q16d_t0=initial_states(39);
    double q17_t0= initial_states(16);      double q17d_t0=initial_states(40);
    double q18_t0= initial_states(17);      double q18d_t0=initial_states(41);
    double q19_t0= initial_states(18);      double q19d_t0=initial_states(42);
    double q20_t0= initial_states(19);      double q20d_t0=initial_states(43);
    double q21_t0= initial_states(20);      double q21d_t0=initial_states(44);
    double q22_t0= initial_states(21);      double q22d_t0=initial_states(45);
    double q23_t0= initial_states(22);      double q23d_t0=initial_states(46);
    double q24_t0= initial_states(23);      double q24d_t0=initial_states(47);


    double q1_tf= final_states(0);          double q1d_tf=  final_states(24);
    double q2_tf= final_states(1);          double q2d_tf=  final_states(25);
    double q3_tf= final_states(2);          double q3d_tf=  final_states(26);
    double q4_tf= final_states(3);          double q4d_tf=  final_states(27);
    double q5_tf= final_states(4);          double q5d_tf=  final_states(28);
    double q6_tf= final_states(5);          double q6d_tf=  final_states(29);
    double q7_tf= final_states(6);          double q7d_tf=  final_states(30);
    double q8_tf= final_states(7);          double q8d_tf=  final_states(31);
    double q9_tf= final_states(8);          double q9d_tf=  final_states(32);
    double q10_tf= final_states(9);         double q10d_tf= final_states(33);
    double q11_tf= final_states(10);        double q11d_tf= final_states(34);
    double q12_tf= final_states(11);        double q12d_tf= final_states(35);
    double q13_tf= final_states(12);        double q13d_tf= final_states(36);
    double q14_tf= final_states(13);        double q14d_tf= final_states(37);
    double q15_tf= final_states(14);        double q15d_tf= final_states(38);
    double q16_tf= final_states(15);        double q16d_tf= final_states(39);
    double q17_tf= final_states(16);        double q17d_tf= final_states(40);
    double q18_tf= final_states(17);        double q18d_tf= final_states(41);
    double q19_tf= final_states(18);        double q19d_tf= final_states(42);
    double q20_tf= final_states(19);        double q20d_tf= final_states(43);
    double q21_tf= final_states(20);        double q21d_tf= final_states(44);
    double q22_tf= final_states(21);        double q22d_tf= final_states(45);
    double q23_tf= final_states(22);        double q23d_tf= final_states(46);
    double q24_tf= final_states(23);        double q24d_tf= final_states(47);

    //Initial     |  //Final
    //states(q)   |  //states(q)
    e(0)=q1_t0;      e(48)=q1_tf;
    e(1)=q2_t0;      e(49)=q2_tf;
    e(2)=q3_t0;      e(50)=q3_tf;
    e(3)=q4_t0;      e(51)=q4_tf;
    e(4)=q5_t0;      e(52)=q5_tf;
    e(5)=q6_t0;      e(53)=q6_tf;
    e(6)=q7_t0;      e(54)=q7_tf;
    e(7)=q8_t0;      e(55)=q8_tf;
    e(8)=q9_t0;      e(56)=q9_tf;
    e(9)=q10_t0;     e(57)=q10_tf;
    e(10)=q11_t0;    e(58)=q11_tf;
    e(11)=q12_t0;    e(59)=q12_tf;
    e(12)=q13_t0;    e(60)=q13_tf;
    e(13)=q14_t0;    e(61)=q14_tf;
    e(14)=q15_t0;    e(62)=q15_tf;
    e(15)=q16_t0;    e(63)=q16_tf;
    e(16)=q17_t0;    e(64)=q17_tf;
    e(17)=q18_t0;    e(65)=q18_tf;
    e(18)=q19_t0;    e(66)=q19_tf;
    e(19)=q20_t0;    e(67)=q20_tf;
    e(20)=q21_t0;    e(68)=q21_tf;
    e(21)=q22_t0;    e(69)=q22_tf;
    e(22)=q23_t0;    e(70)=q23_tf;
    e(23)=q24_t0;    e(71)=q24_tf;

    //Initial     |  //Final
    //states(qd)  |  //states(qd)
    e(24)=q1d_t0;    e(72)=q1d_tf;
    e(25)=q2d_t0;    e(73)=q2d_tf;
    e(26)=q3d_t0;    e(74)=q3d_tf;
    e(27)=q4d_t0;    e(75)=q4d_tf;
    e(28)=q5d_t0;    e(76)=q5d_tf;
    e(29)=q6d_t0;    e(77)=q6d_tf;
    e(30)=q7d_t0;    e(78)=q7d_tf;
    e(31)=q8d_t0;    e(79)=q8d_tf;
    e(32)=q9d_t0;    e(80)=q9d_tf;
    e(33)=q10d_t0;   e(81)=q10d_tf;
    e(34)=q11d_t0;   e(82)=q11d_tf;
    e(35)=q12d_t0;   e(83)=q12d_tf;
    e(36)=q13d_t0;   e(84)=q13d_tf;
    e(37)=q14d_t0;   e(85)=q14d_tf;
    e(38)=q15d_t0;   e(86)=q15d_tf;
    e(39)=q16d_t0;   e(87)=q16d_tf;
    e(40)=q17d_t0;   e(88)=q17d_tf;
    e(41)=q18d_t0;   e(89)=q18d_tf;
    e(42)=q19d_t0;   e(90)=q19d_tf;
    e(43)=q20d_t0;   e(91)=q20d_tf;
    e(44)=q21d_t0;   e(92)=q21d_tf;
    e(45)=q22d_t0;   e(93)=q22d_tf;
    e(46)=q23d_t0;   e(94)=q23d_tf;
    e(47)=q24d_t0;   e(95)=q24d_tf;


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


