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
 *	\file include/openhrc/dynamics/Lie_operators.h
 *	\author Alvaro Paz, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2018
 *
 *	Definition of Lie operators.
 */

#ifndef HR_DYNAMICS_LIE_OPERATORS_H
#define HR_DYNAMICS_LIE_OPERATORS_H

#include "openhrc/types.h"
#include "openhrc/core.h"

#define PI 3.141592653589793

using std::cout;
using std::cin;
using std::endl;
using std::string;

namespace hr{
namespace core{

class LieOperators
{

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------


    public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW //128-bit alignment

    //! Default constructor.
    /*! \param .
        \param */
    LieOperators(){}


    // --------------------------------------------
    // Members
    // --------------------------------------------


    // --------------------------------------------
    // Methods
    // --------------------------------------------

    public:
        /*! Description: Computes the 3x3 skew symmetric matrix from a 3D vector
             * \parameters: 3D vector
             * \return: Skew symmetric matrix
             */
        Matrix3r skew(const Vector3r &vec);

        /*! Description: Upgrades the screw vector into its matrix form
             * \parameters: Spatial screw vector
             * \return: Screw matrix representation
             */
        Matrix4r upScrew(const SpatialVector &screw);

        /*! Description: Computes the Adjoint operator Ad: se(3) -> se(3)
             * \parameters: An element of the group SE(3)
             * \return: Adjoint operator
             */
        SpatialMatrix Ad(const Matrix4r &G);

        /*! Description: Computes the dual Adjoint operator Ad*: se*(3) -> se*(3)
             * \parameters: An element of the group SE(3)
             * \return: Dual Adjoint operator
             */
        SpatialMatrix AdDual(const Matrix4r &G);

        /*! Description: Computes the inverse dual Adjoint operator Ad*: se*(3) -> se*(3)
             * \parameters: An element of the group SE(3)
             * \return: Inverse dual Adjoint operator
             */
        SpatialMatrix AdDualInv(const Matrix4r &G);

        /*! Description: Computes the adjoint operator ad: se(3) -> se(3)
             * \parameters: An element of the algebra se(3)
             * \return: adjoint operator
             */
        SpatialMatrix ad(const SpatialVector &V);

        /*! Description: Computes the dual adjoint operator ad*: se*(3) -> se*(3)
             * \parameters: An element of the algebra se(3)
             * \return: adjoint operator
             */
        SpatialMatrix adDual(const SpatialVector &V);

        /*! Description: Computes the bar adjoint operator ad^{_}
             * \parameters: An element of the algebra se*(3)
             * \return: adjoint operator
             */
        SpatialMatrix adBar(const SpatialVector &F);

};
} // end of namespace core
} // end of namespace hr

#endif // HR_DYNAMICS_LIE_OPERATORS_H
