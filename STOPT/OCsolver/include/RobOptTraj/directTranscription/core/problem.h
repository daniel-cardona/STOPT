#include "iostream"
#include <string>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <fstream>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"


//----OPEN_HRC dependences
#include "openhrc/core.h"
#include "openhrc/dynamics.h"
#include <iomanip>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
//-----------------------

#include "IpTNLP.hpp"

using namespace std;

#define INF (1e19)
#define MC_EPSILON 2.221e-16

#ifndef MAX
#define MAX(a, b) ( (a)>(b)?  (a):(b) )
#endif
#ifndef MIN
#define MIN(a, b) ( (a)<(b)?  (a):(b) )
#endif


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//----------------------------------------STRUCTURES---------------------------------------------------
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//******************Bounds structures******************

struct ulb_str{

    Eigen::VectorXd lower;
    Eigen::VectorXd upper;

};

typedef struct ulb_str ulbounds;


//******************Variable bounds******************

struct bnd_str{

  ulbounds states;
  ulbounds controls;
  ulbounds events;
  ulbounds path;
  ulbounds startTime;
  ulbounds finalTime;

};

typedef struct bnd_str Bounds;

//******************Guess structure******************


struct str_guess{

    Eigen::MatrixXd states;
    Eigen::MatrixXd controls;
    Eigen::MatrixXd time;

};

typedef  struct str_guess Guess;

//******************Scaling structure******************

struct str_scale{

    Eigen::VectorXd controls;
    Eigen::VectorXd states;
    Eigen::MatrixXd defects;
    Eigen::VectorXd defectsV;
    Eigen::VectorXd events;
    Eigen::MatrixXd path;

    Eigen::SparseMatrix<double> invScaleMatrix;     //Inverse of the scaling factors matrix

    double timeCns;
    double time;
    double cost_fcn;

};
typedef  struct str_scale Scale;


//******************Sparsity templates***************

struct str_sparsity{

    Eigen::SparseMatrix<double> A;

    Eigen::SparseMatrix<double> B;

    Eigen::SparseMatrix<double> B2;

    Eigen::MatrixXd positivePerturbationMatrix;

    Eigen::MatrixXd negativePerturbationMatrix;


};

typedef  struct str_sparsity sparsityTemplates;


//******************Solution structure******************

struct str_solution{

    Eigen::MatrixXd controls;
    Eigen::MatrixXd states;
    Eigen::Vector2d time;

};
typedef struct str_solution Solution;


//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

namespace OCSolver
{

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//--------------------------------------PROBLEM CLASS--------------------------------------------------
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//!  Optimal Control Problem information class definition
/*!  \ingroup core_module
        The Optimal Control Problem information class is used as a representation of the problem data of an OPC
    */

class Problem{

    //---------------------------------------------------------------------
    //------------------------------- MEMBERS -----------------------------
    //---------------------------------------------------------------------


public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //*********************************************************************
    //******************Optimal control problem variables******************
    //*********************************************************************

    int nStates;             //Number of states of the system
    int nControls;           //Number of control inputs in the systems
    int nEvents;             //Number of event constraints
    int nPath;               //Number of path constraints
    int nNodes;              //Number of discrete points

    int interpolationPoints; //Mesh Variables

    Bounds bounds;           //Bound of the variables and the constraints

    Guess guess;             //Initial guess of the problem states and controls

    Scale scale;             //Scaling factors of the problem states,controls and cns

    string discretizationMethod; //Discretization method used for the defect cns and the cost function

    Solution solution;       //Solution of the OCP

    bool use_constraint_scaling; //Flag used in the scaling process


    //**************Pointers to the OCP user-defined functions********

    double (*endpoint_cost)(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tf, Problem &problem);

    double (*integrand_cost)(Eigen::VectorXd &states,Eigen::VectorXd &controls,double &tk, Problem &problem);

    void (*dae)(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &t, Eigen::VectorXd &derivatives,Eigen::VectorXd &path, Problem &problem);

    void (*events)(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tF,Eigen::VectorXd &e, Problem &problem);


    //*********************************************************************
    //********************** NLP problem variables ************************
    //*********************************************************************

    int nDecVar;            //Number of decision variables generated by the transcription process
    int nCns;               //Number of constraints generated by the transcription process
    int nDefectCns;         //Number of defect constraints generated by the dynamics of the system
    int nSegments;          //Number of segments given by the discretization method selected by the user

    Eigen::VectorXd snodes; //Vector used for the computation of the time vector
    Eigen::VectorXd x0;     //Initial guess vector for the NLP Problem
    Eigen::VectorXd xlb;    //Lower bound of the decision variables
    Eigen::VectorXd xub;    //Upper bound of the decision variables

    Eigen::VectorXd glb;    //Lower bounds of the NLP constraints
    Eigen::VectorXd gub;    //Upper bounds of the NLP constraints

    Eigen::VectorXd cnsScaling; //Vector with the scaling factors of the constraints

    //*********************************************************************
    //*************************Sparsity structures ************************
    //*********************************************************************

    int nnz;                        //Number of nonzero elements in the jacobian matrix
    Eigen::MatrixXd nzG;            //Matrix with the information of the non-zero elements of the jacobian matrix [row,column]

    Eigen::MatrixXd idxGroups;
    Eigen::VectorXd sizeGroups;

    int nnzD;                       //Number of nonzero elements in the derivative matrix
    Eigen::MatrixXd nzGD;           //Matrix with the information of the non-zero elements of the derivative matrix [row,column]

    Eigen::MatrixXd idxGroupsD;     //Index groups of the derivative matrix
    Eigen::VectorXd sizeGroupsD;    //Size of the index groups

    sparsityTemplates sparsity;     //Sparsity templates needed for the problem to compute the jacobian matrix using sparse algebra


    //*********************************************************************
    //********************** Function variables ***************************
    //*********************************************************************

    Eigen::MatrixXd decVarMatrix;    //Matrix with the unscaled decision variables of the problem
    Eigen::VectorXd kStates;         //Vector with the states of the system in k
    Eigen::VectorXd kControls;       //Vector with the controls of the system in k
    Eigen::VectorXd dx;             //Vector with the derivatives of the states of the system in k
    Eigen::VectorXd kPath;           //Vector with the value of the path constraints in k


    //----New version variables-----

    Eigen::VectorXd states_k1;      //Vector to store the state of the system in k+1 (Useful for cost functions that requires the use of the next state values)


    //*********************************************************************
    //************* Analytical derivatives variables **********************
    //*********************************************************************

    Eigen::MatrixXd dxGradient;
    Eigen::MatrixXd pathGradient;
    Eigen::MatrixXd eventGradient_t0;
    Eigen::MatrixXd eventGradient_tF;
    Eigen::MatrixXd derivativeMatrix;

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> J;


    //*********************************************************************
    //************* Open_HRC variables **********************
    //*********************************************************************

    std::shared_ptr<hr::core::MultiBody> robot;


    //---------------------------------------------------------------------
    //------------------------------- METHODS -----------------------------
    //---------------------------------------------------------------------

    //!  Compute the size of the NLP problem and initialize the size of the vectors and matrix needed for the transcription process.
    /*! \param none
     * \returns void
     */
    void setup();



};

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//------------------------------------ALGORITHM CLASS--------------------------------------------------
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    class Alg{

    public:

      int nlp_iter_max;
      double nlp_tolerance;
      double EPS;
      string scaling;
      string defect_scaling;
      string derivatives;

      bool sparsityOptimization; //This feature is temporal (Only for test version)

      bool derivativeChecker;

    };

}

static const double pi = 3.141592653589793;
