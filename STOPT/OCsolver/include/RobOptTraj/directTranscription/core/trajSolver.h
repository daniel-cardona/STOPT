#include "problem.h"

namespace OCSolver {

    namespace core {

    //! Main function of the solver that call the functions to solve the problem
    //!
    //! This function is stored in: src/solver/trajSolver.cpp
    //!
    /*! \param $Prob$              Object with the problem information
     *  \param $Alg$               Object with the algorithm options
     *
     *  \return $solution.txt$     File with the NLP solution
     */
    void robTrajSol(Problem &problem,Alg &algorithm);

    //! Function that reads the .txt file with the NLP solution and
    //! make the interpolation to generate a continuos solution
    //!
    //! This function is stored in: src/solver/solution.cpp
    //!
    /*! \param $Prob$   Object with the problem information
     *  \param $Alg$    Object with the algorithm options
     *
     *  \return .txt file with the final solution
     */

    void generateSol(Problem &problem, Alg & algorithm);

    //! Function that validates the basic information of the problem defined by the user, if some infeasible options
    //! are detected by this function an exception is thrown
    //!
    //! This function is stored in: src/solver/trajSolver.cpp
    //!
    /*! \param $Prob$   Object with the problem information
     *  \param $Alg$    Object with the algorithm options
     *
     *  \return void
     */

    void validateUserInput(Problem &problem, Alg &algorithm);




    void displayInformation(Problem &problem,Alg &algorithm);

    }

}
