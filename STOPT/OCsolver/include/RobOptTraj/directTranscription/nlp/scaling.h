#include "RobOptTraj/directTranscription/nlp/setup.h"


namespace OCSolver
{

namespace NLP
{

        //! Function that obtains and set the scaling factor of the decision variables for the NLP problem
        //!
        //!
        /*! \param $Prob$              Object with the problem information
         *  \param $Alg$               Object with the algorithm options
         *
         *  \return void
         */

        void determineScalingFactorsDecVariables(Problem &problem,Alg &algorithm);


        //! Function that obtains and set the scaling factor of the cost function for the NLP problem
        //!
        //!
        /*! \param $Prob$              Object with the problem information
         *  \param $Alg$               Object with the algorithm options
         *
         *  \return void
         */

        void determineObjectiveScaling(Problem &problem,Alg &algorithm);

        //! Function that obtains and set the scaling factor of the constraints
        //!
        //!
        /*! \param $Prob$              Object with the problem information
         *  \param $Alg$               Object with the algorithm options
         *
         *  \return void
         */

        void determineConstraintScaling(Problem &problem,Alg &algorithm);

}

}
