#include "RobOptTraj/directTranscription/collocationMethods/hermiteSimpson.h"

namespace OCSolver
{

namespace NLP {

    //! Function that transcribes and sets the initial guess in such a way that can be
    //! used by the NLP solver
    //!
    /*! \param $Prob$   Object with the problem information
     *  \param $Alg$    Object with the algorithm options
     *
     *  \return void
     */

       void setInitialGuess(Problem &problem, Alg &algorithm);

    //! Function that transcribes and sets the variable bounds in such a way that can be
    //! used by the NLP solver
    //!
    //! This function is stored in: src/transcription/NLP_Bounds.cpp
    //!
    /*! \param $Prob$   Object with the problem information
     *  \param $Alg$    Object with the algorithm options
     *
     *  \return void
     */

       void setVariableBounds(Problem &problem, Alg &algorithm);


    //! Function that transcribes and sets the constraints bounds(sides) in such a way that can be
    //! used by the NLP solver (IPOPT)
    //!
    //! This function is stored in: src/transcription/NLP_Bounds.cpp
    //!
    /*! \param $Prob$   Object with the problem information
     *  \param $Alg$    Object with the algorithm options
     *
     *  \return void
     */

       void setConstraintBounds(Problem &problem,Alg &algorithm);






}

}
