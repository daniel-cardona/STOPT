/*
 *
 * Copyright (C) 2015
 * Julio Jarquin <jjarquin.ct@gmail.com>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx>
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
 *	\file include/openhrc/types.h
 *	\author Julio Jarquin, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2015
 *
 *	Definition of Types
 */

#ifndef HR_TYPES_H
#define HR_TYPES_H

#include "Eigen/Dense"
#include "Eigen/../unsupported/Eigen/KroneckerProduct"
//#include "Eigen/src/Core/Map.h"
#include <chrono>


//! Namespace hr. Main namespace of OpenHRC
/*!
    More information about hr.
*/
namespace hr{

//! Defines TYPE real_t for facilitating the switching between double and float.
//#ifdef  __USE_SINGLE_PRECISION__
//        using real_t = float;  //c++11
//! Defines real_t TYPE as float.
//typedef float real_t;
//#else
//! Defines real_t TYPE as double.
typedef double real_t;
//#endif

//! Spatial algebra Vector, 6-d vector for spatial transformations.
typedef Eigen::Matrix< real_t, 6, 1> SpatialVector;

//! Spatial algebra Matrix, 6x6 matrix for spatial transformations.
typedef Eigen::Matrix< real_t, 6, 6> SpatialMatrix;

//! Derivative of Spatial algebra Matrix, 6x6 matrix for spatial transformations.
typedef Eigen::Matrix< real_t, 6, Eigen::Dynamic> D_SpatialVector;

//! Derivative of scalar real_t
typedef Eigen::Matrix< real_t, 1, Eigen::Dynamic> D_real_t;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t row major matrix. RowMajor defines the way the matrix is stored
 */
typedef  Eigen::Matrix<real_t,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXr;
typedef  Eigen::Matrix<real_t,Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXrColMajor;


//! Defines a new eigen-based vector type that is based on the real_t type.
/*!
  Defines a real_t vector.
 */
typedef  Eigen::Matrix<real_t,Eigen::Dynamic, 1> VectorXr;

typedef  Eigen::Matrix<real_t, 1, Eigen::Dynamic> RowVectorXr;


//! Defines a new eigen-based vector type that is based on the real_t type.
/*!
  Defines a real_t vector of size 3x1.
 */
typedef Eigen::Matrix<real_t , 3 , 1> Vector3r;


//! Defines a new eigen-based vector type that is based on the real_t type.
/*!
  Defines a real_t vector of size 2x1.
 */
typedef Eigen::Matrix<real_t , 2 , 1> Vector2r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 4x4.
 */
typedef Eigen::Matrix<real_t , 4 , 4> Matrix4r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 3x3.
 */
typedef Eigen::Matrix<real_t , 3 , 3> Matrix3r;


//! Defines a new eigen-based matrix type that is based on the real_t type.
/*!
  Defines a real_t matrix of size 2x2.
 */
typedef Eigen::Matrix<real_t , 2 , 2> Matrix2r;

}

#endif // TYPES_H
