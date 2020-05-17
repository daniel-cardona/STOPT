/*
 *
 * Copyright (C) 2019
 * Daniel S. Cardona <ingdanielcardona@gmail.com>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx>
 * CINVESTAV - Saltillo Campus
 *
 * This file is part of OpenHRC
 * OpenHRC is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * OpenHRC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

/**
 *	\file src/problem.h
 *	\author Daniel S. Cardona
 *	\version 1.0
 *	\date 2019
 *
 * Problem and algorithm class implementation.
 */

#include "RobOptTraj/directTranscription.h"
#include "IpIpoptApplication.hpp"

namespace OCSolver
{

namespace NLP
{

namespace ipopt {


// constructor
IPOPT_INTERFACE::IPOPT_INTERFACE(Problem &data,Alg &algdata)
{
    problem=data;
    algorithm=algdata;
}

// destructor
IPOPT_INTERFACE::~IPOPT_INTERFACE()
{ }



//------------------Set up the size of the problem!-----------------------

bool IPOPT_INTERFACE::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
   )
{
   // Obtain the number of decision variables
   int nDecVar=problem.nDecVar;
   int nConstraints=problem.nCns;
   int nzJac=problem.nnz;


   //Number of decision variables

   n = nDecVar;

   // Number of equality contraints
   m = nConstraints;

   // Number of non-zero elements
   nnz_jac_g = nzJac;

   //Using the Quasi-Newton Approximation of Second Derivatives (this option should be disabled)
   nnz_h_lag = 10;

   // Use the C style indexing (0-based)
   index_style = TNLP::C_STYLE;

   return true;
}


//----------------Set the variable and constraints bounds-------------------------


bool IPOPT_INTERFACE::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
   )
{
   // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
   // If desired, we could assert to make sure they are what we think they are.
   assert(n == problem.nDecVar);
   assert(m == problem.nCns);

   //Obtain the defined lower and upper bounds

   Eigen::VectorXd lb(problem.nDecVar);
   Eigen::VectorXd ub(problem.nDecVar);

   lb=problem.xlb;
   ub=problem.xub;


   // Setting the lower bounds
   for( Index i = 0; i < problem.nDecVar; i++ )
   {
      x_l[i] = lb(i);


   }
   // Setting the upper bounds
   for( Index i = 0; i < problem.nDecVar; i++ )
   {
      x_u[i] = ub(i);

   }


   // Set the upper and lower bounds of the contraints

   Eigen::VectorXd cns_lb(problem.nCns);
   Eigen::VectorXd cns_ub(problem.nCns);

   cns_lb=problem.glb;
   cns_ub=problem.gub;

   for( Index i = 0; i < problem.nCns; i++ )
   {
       g_l[i] = cns_lb(i);
       g_u[i] = cns_ub(i);

   }
   return true;
}

//------------Set the initial guess of the problem--------------------

bool IPOPT_INTERFACE::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
   )
{
   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   // initialize to the given starting point
   Eigen::VectorXd iGuess(problem.nDecVar);

   iGuess<<problem.x0;

   for( Index i = 0; i < problem.nDecVar; i++ )
   {
    x[i]=iGuess(i);
   }

   return true;
}

//----------------- Cost function------------------------

// returns the value of the objective function
bool IPOPT_INTERFACE::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
   )
{
   assert(n == problem.nDecVar);

   Eigen::VectorXd z(problem.nDecVar);

   //Let's obtain the decVar vector in such a way that can be used for EIGEN!
   for( Index i = 0; i < problem.nDecVar; i++ )
   {
    z(i)=x[i];

   }

   double costObjValue=0;

   if(problem.discretizationMethod=="Trapezoidal")
   {
       costObjValue=collocationMethod::trapezoidal::costFunction(z,problem,algorithm);
    }

   obj_value =costObjValue;

   return true;
}

//--------------------Gradient of the objective function--------------

bool IPOPT_INTERFACE::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
   )
{
   assert(n == problem.nDecVar);

   //Let's obtain the decVar vector in such a way that can be used for EIGEN!

   Eigen::VectorXd z(problem.nDecVar);
   for( Index i = 0; i < problem.nDecVar; i++ )
   {
    z(i)=x[i];
   }

   //Obtain the gradient of the cost function
   Eigen::VectorXd costObjValueGrad;

   if(problem.discretizationMethod=="Trapezoidal")
   {
       derivatives::numerical::scalarGradient(collocationMethod::trapezoidal::costFunction,z,problem,algorithm,costObjValueGrad);
   }


   for( Index i = 0; i < problem.nDecVar; i++ )
   {
    grad_f[i]=costObjValueGrad(i);

   }


   return true;
}

//-----------------------Value of the constraints------------------------

// return the value of the constraints: g(x)
bool IPOPT_INTERFACE::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
   )
{
   assert(n == problem.nDecVar);
   assert(m == problem.nCns);


   Eigen::VectorXd z(problem.nDecVar);

   //Let's obtain the decVar vector in such a way that can be used for EIGEN!
   for( Index i = 0; i < problem.nDecVar; i++ )
   {
    z(i)=x[i];

   }

   Eigen::VectorXd cns(problem.nCns);

   if(problem.discretizationMethod=="Trapezoidal"){

          collocationMethod::trapezoidal::cnsFunction(z,cns,problem,algorithm);
   }

   for( Index i = 0; i < problem.nCns; i++ )
   {
    g[i]=cns(i);
   }

   return true;
}

//-------------- Jacobian of the constraints ----------------------

bool IPOPT_INTERFACE::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
   )
{
    assert(nele_jac == problem.nnz);
    assert(n == problem.nDecVar);
    assert(m == problem.nCns);

   if( values == NULL )
   {

       //------------NEW VERSION MODIFICATION----------

      // return the structure of the Jacobian

        int nnzG=problem.nnz;

        int row, col;
       for(int i = 0; i < nnzG; i++)
       {
         row=problem.nzG(i,0);
         col=problem.nzG(i,1);

        iRow[i]=row;
        jCol[i]=col;
       }


   }
   else
   {


       Eigen::VectorXd z(problem.nDecVar);

       //int nnzA=problem.nnzA;
       //int nnzG=problem.nnzG;

       //Let's obtain the decVar vector in such a way that can be used for EIGEN!
       for( Index i = 0; i < problem.nDecVar; i++ )
       {

        z(i)=x[i];

       }

       //Obtain the gradient of the jacobian
       Eigen::VectorXd nzValues(problem.nnz);

       if(problem.discretizationMethod=="Trapezoidal"){

           if(algorithm.derivatives=="numerical"){

               derivatives::numerical::computeNumericalJacobian(collocationMethod::trapezoidal::rightHandSideSparsity,z,problem,algorithm,nzValues);

           }

           if(algorithm.derivatives=="analytical"){

               derivatives::analytical::computeAnalyticalJacobian(z,problem,algorithm,nzValues);
           }

       }

      // return the values of the Jacobian of the constraints

       for( Index i = 0; i < problem.nnz; i++ )
       {
            values[i]=nzValues(i);
       }


   }
   return true;
}

//----------------------Evaluation of the Hessian------------------------------

//return the structure or values of the Hessian
bool IPOPT_INTERFACE::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
   )
{
   return true;
}

//-------------------------Print the solution-----------------------------------

void IPOPT_INTERFACE::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              x,
   const Number*              z_L,
   const Number*              z_U,
   Index                      m,
   const Number*              g,
   const Number*              lambda,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
   )
{
   // here is where we would store the solution to variables, or write to a file, etc
   // so we could use the solution.

   // For this example, we write the solution to the console


   //Let's return this solution and return the control to the solver!

   Eigen::VectorXd sol(problem.nDecVar);

   for( Index i = 0; i < problem.nDecVar; i++ )
   {
      sol(i)=x[i];
   }

   //Obtain the real solution
    /*
   std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;

   std::cout << "t(0) = " << sol(0) << std::endl;
   std::cout << "t(F) = " << sol(1)<< std::endl;

   std::cout << std::endl ;

   for( Index i = 2; i < n; i++ )
   {
      std::cout << "x[" << i << "] = " << sol(i) << std::endl;;
   }

   std::cout << std::endl << std::endl << "Objective value" << std::endl;
   std::cout << "f(x*) = " << obj_value/problem.scale.cost_fcn << std::endl;

   std::cout << std::endl << "Final value of the constraints:" << std::endl;
   for( Index i = 0; i < m; i++ )
   {
      std::cout << "g(" << i << ") = " << g[i] << std::endl;
S
   }
*/
   //Save the solution

      std::ofstream file("solution.txt");
        if (file.is_open())
        {
          file << sol << '\n';
        }

      Eigen::MatrixXd decVarMatrix;

      Eigen::VectorXd xk(problem.nStates);
      Eigen::VectorXd uk(problem.nControls);

      Eigen::MatrixXd states(problem.nStates,problem.nNodes);
      Eigen::MatrixXd controls(problem.nControls,problem.nNodes);

      double t0,tf;

      utils::getDecVarMatrix(problem,algorithm,sol,decVarMatrix,t0,tf);

      for(int k=0;k<problem.nNodes;k++){

          utils::getVariables(problem,decVarMatrix,xk,uk,k);

          states.col(k)=xk;

          controls.col(k)=uk;

      }

      std::ofstream file2("states_solution.txt");
        if (file2.is_open())
        {
          file2 << states << '\n';
        }

       std::ofstream file3("controls_solution.txt");
          if (file3.is_open())
          {
            file3 << controls << '\n';
          }

}






int solve(Problem &problem,Alg &algorithm){

    double sqreps=sqrt(MC_EPSILON);

    double delj=sqreps;


    // Create a new instance of your nlp
        //  (use a SmartPtr, not raw)
        SmartPtr<TNLP> mynlp = new IPOPT_INTERFACE(problem,algorithm);

        // Create a new instance of IpoptApplication
        //  (use a SmartPtr, not raw)
        // We are using the factory, since this allows us to compile this
        // example with an Ipopt Windows DLL
        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->RethrowNonIpoptException(true);

        // Change some options
        // Note: The following choices are only examples, they might not be
        //       suitable for your optimization problem.
        app->Options()->SetNumericValue("tol", 1e-4);
        app->Options()->SetIntegerValue("max_iter",4000);
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("output_file", "ipopt.out");
        app->Options()->SetStringValue("hessian_approximation","limited-memory");
        app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
        //app->Options()->SetStringValue("jacobian_approximation","finite-difference-values");
        app->Options()->SetIntegerValue("print_level",5);
        //app->Options()->SetStringValue("derivative_test_print_all","yes");

        //Set the derivative test

        if(algorithm.derivativeChecker==true){

            app->Options()->SetNumericValue("derivative_test_perturbation",delj);
            app->Options()->SetStringValue("derivative_test","first-order");
        }


        // Initialize the IpoptApplication and process the options
        ApplicationReturnStatus status;
        status = app->Initialize();
        if( status != Solve_Succeeded )
        {
           std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
           return (int) status;
        }

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        if( status == Solve_Succeeded )
        {
           std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
        }
        else
        {
           std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
        }



}



}


}//END NLP namespace

}//END OCSolver namespace


