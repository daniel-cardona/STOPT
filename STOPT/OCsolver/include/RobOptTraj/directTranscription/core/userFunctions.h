#include "trajSolver.h"


namespace OCSolver {

    //! User defined function that must return the Mayer cost(if defined) fo the objective function
    //!
    //! This function is stored in: src/examples/.../OCP_definition.cpp
    //!
    /*! \param $x0$      Vector with initial states
     *  \param $xF$      Vector with final states
     *  \param $t0$      Initial time
     *  \param $tF$      Final time
     *
     *  \return Mayer cost
     */

    double endpoint_cost(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tf,Problem &problem);


    //! User defined function that must return the Lagrange cost(if defined) of the objective function
    //!
    //! This function is stored in: src/examples/.../OCP_definition.cpp
    //!
    /*! \param $states$      Vector with states in tk
     *  \param $controls$    Vector with controls in tk
     *  \param $tk$          Time at the k discretization node
     *  \param $w$           Weight vector
     *
     *  \return Lagrange cost
     */

    double integrand_cost(Eigen::VectorXd &states,Eigen::VectorXd &controls,double &tk,Problem &problem);


    //! User defined function that must return the derivatives of the state vector(this is generally defined by the dynamics of the system)
    //! and the path constraints
    //!
    //! This function is stored in: src/examples/.../OCP_definition.cpp
    //!
    /*! \param $states$         Vector with states in tk
     *  \param $controls$       Vector with controls in tk
     *  \param $tk$             Time at the k discretization node
     *
     *  \return $derivatives$   Vector with the derivative of the states at tk
     *  \return $path$          Vector with the values of the path contraints at tk
     */

    void dae(Eigen::VectorXd &states, Eigen::VectorXd &controls,double &tk, Eigen::VectorXd &derivatives,Eigen::VectorXd &path,Problem &problem);


    //! User defined function that must return a vector with the value of the event constraints
    //!
    //! This function is stored in: src/examples/.../OCP_definition.cpp
    //!
    /*! \param $x0$      Vector with initial states
     *  \param $xF$      Vector with final states
     *  \param $t0$      Initial time
     *  \param $tF$      Final time
     *
     *  \return $e$   Vector with the value of the event constraints
     *
     */

    void events(Eigen::VectorXd &x0,Eigen::VectorXd &xF,double &t0,double &tF, Eigen::VectorXd &e,Problem &problem);




}



