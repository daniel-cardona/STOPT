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
 *	\file include/openhrc/core/types.h
 *	\author Julio Jarquin, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2015
 *
 *	Definition of Types
 */

#ifndef HR_CORE_TYPES_H
#define HR_CORE_TYPES_H

#include "openhrc/types.h"


namespace hr{

//! Namespace core. Namespace of the core Module
/*!
    More information about core.
*/
typedef Eigen::Matrix<real_t, 3 , Eigen::Dynamic , Eigen::RowMajor> Matrix3Xr;

namespace core{

    //! Type definition of Eigen::AngleAxis< Scalar > as Eigen::AngleAxis< real_t >
    typedef Eigen::AngleAxis<real_t> AngleAxisr;


//! Namespace nao. Namespace for nao settings
/*!
    More information about nao.
*/
namespace nao{

enum body{

            kLeftFoot = 24,
            kRightFoot = 30,
            kLeftHand = 13,
            kRightHand = 18,
            kTorso = 6,
            kHead = 8
};

enum Pose{
    kStandZero,
    kStandInit
};


} // end of namespace nao

} // end of namespace core

} // end of namespace hr


#endif // HR_CORE_TYPES_H
