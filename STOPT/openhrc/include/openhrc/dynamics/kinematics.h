/*
 *
 * Copyright (C) 2018
 * Gustavo Arechavaleta <garechav@cinvestav.edu.mx>, Alvaro Paz
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
 *	\file include/openhrc/dynamics/kinematics.h
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2018
 *
 *	Class to implement the multibody kinematics.
 */

#ifndef HR_DYNAMICS_KINEMATICS_H
#define HR_DYNAMICS_KINEMATICS_H

#include <memory>
#include "openhrc/types.h"
#include "openhrc/core.h"
#include "openhrc/dynamics/Lie_operators.h"

namespace hr{
namespace core{


class Kinematics : public LieOperators
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------

    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Custom constructor.
    Kinematics(std::shared_ptr< MultiBody > robot){

        this->robot = robot;
        this->n = robot->getDoF();

        this->q = VectorXr::Zero(n,1);
        this->dq = VectorXr::Zero(n,1);
        this->ddq = VectorXr::Zero(n,1);

        this->Identity3x3.setIdentity();
        this->e_sq.setIdentity();
        this->bodies = robot->getPubBodies();

        //! Variables for differentiation //!<$ Default <- wrt_state

        this->D_X = MatrixXr::Identity(2*n,2*n);
        this->D_q = D_X.block(0,0,n,2*n);
        this->D_dq = D_X.block(n,0,n,2*n);
        this->D_ddq = MatrixXr::Zero(n,2*n);

        SK.clear();    SK_2.clear();
        //! Forward recursion
        for ( bodyIterator = bodies.begin()+1; bodyIterator != bodies.end(); bodyIterator++ )
        {
            tempBody = (DynamicBody *)(*bodyIterator);

            sk = skew(tempBody->getScrewAxes().segment(3,3));

            SK.push_back(sk);
            SK_2.push_back(sk*sk);

        }

    }


    ~Kinematics(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------

    protected:

    //! Dof
    short int n;

    //! MultiBody
    std::shared_ptr< MultiBody > robot;


    ///! Initialize variables
    VectorXr q, dq, ddq;
    Matrix3r sk, sk_2, e_wq, T_wq, Identity3x3;

    Matrix4r e_sq;

    typename std::vector< Body* >::iterator bodyIterator;
    DynamicBody* tempBody;
    DynamicBody* parentBody;
    std::vector<Body*> bodies;
    std::vector<short int> bodyChildren;
    DynamicBody* children;
    typename std::vector< short int >::iterator childrenIterator;

    //! Variables for differentiation

    MatrixXr D_X, D_q, D_dq, D_ddq;

    std::vector< Eigen::Matrix<real_t, 3, 3> > SK, SK_2;

    public:



    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:

    //! Compute the whole multibody forward kinematics
         /*! \param null
         * \return void
         */
    void computeForwardKinematics( );    

    //! Compute the predeccessors forward kinematics
         /*! \param DynamicBody*
         * \return void
         */
    void computeForwardKinematics( DynamicBody* currentBody );

    //! Set generalized coordinates
         /*! \param configuration, generalized velocity, generalized acceleration
         * \return void
         */
    void setGeneralizedCoordinates( const VectorXr &q, const VectorXr &dq, const VectorXr &ddq );

    //! Set generalized coordinates differentiation
         /*! \param D configuration, D generalized velocity, D generalized acceleration
         * \return void
         */
    void setGeneralizedCoordinatesDifferentiation( const MatrixXr &D_q, const MatrixXr &D_dq, const MatrixXr &D_ddq );


protected:


};
} // end of namespace core
} // end of namespace hr

#endif // HR_DYNAMICS_KINEMATICS_H
