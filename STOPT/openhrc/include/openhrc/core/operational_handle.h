/*
 *
 * Copyright (C) 2015
 * Julio Jarquin <jjarquin.ct@gmail.com>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx> , Gerardo Jarquin
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
 *	\file include/openhrc/core/operational_handle.h
 *	\author Gerardo Jarquin, Gustavo Arechavaleta, Julio Jarquin.
 *	\version 1.0
 *	\date 2015
 *
 * OperationalHandle class.
 */

#ifndef HR_CORE_OPERATIONAL_HANDLE_H
#define HR_CORE_OPERATIONAL_HANDLE_H

#include <string>
#include "spatial_algebra.h"

namespace hr {
namespace core {

//! OperationalHandle class definition
/*! \ingroup core_module
    This class represents a coordinate frame that can be attached to
    a given Body or and object of the environment. \sa Body, Joint
    */
class OperationalHandle {

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:

    //! Custom constructor.
    /*!
        \param id The id of the operational handle.
        \param name The name of the operational handle.
      */
    OperationalHandle(const unsigned int id, const std::string &name = "End-effector");

    //! Custom constructor.
    /*!
        \param id The id of the operational handle.
        \param pos The position of the operational handle, relative to the associated body.
        \param name The name of the operational handle.

      */
    OperationalHandle(const unsigned int id, const Vector3r &pos,
                      const std::string &name = "End-effector");

    //! Custom constructor.
    /*!
        \param id The id of the operational handle.
        \param name The name of the operational handle.
        \param pos The position of the operational handle, relative to the associated body
        \param rot The corresponding rotation (Rx,Ry,Rz) of the operational handle.
      */
    OperationalHandle(const unsigned int id, const std::string &name,
                      const Vector3r &pos, const Matrix3r &rot);

    //! Default destructor.
    ~OperationalHandle();

    // --------------------------------------------
    // Members
    // --------------------------------------------
protected:
    //! OperationalHandle identifier
    unsigned int id;

    //! OperationalHandle name
    std::string name;

    //! OperationalHandle position (x,y,z) relative to the associated body
    Vector3r position;

    //! OperationalHandle rotation matrix (Rx,Ry,Rz)
    /*! This matrix represents the orientation relative to the associated body
     */
    Matrix3r orientation;

    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:
    //! Get OperationalHandle ID
     /*! \param none
       \returns short int identifying the opHandle.
     */
    short int getId(){ return id; }

    //! Set OperationalHandle ID
     /*! \param short int identifying the opHandle
         \returns void
     */
    void setId(const int &id){ this->id = id; }

    //! Get OperationalHandle name
     /*! \param name
      \returns string with the opHandle name.
     */
    std::string getName(){ return this->name; }

    //!  Set OperationalHandle name
     /*! \param string with the opHandle name.
         \returns void
     */
    void setName(const std::string &name){ this->name = name; }

    //! Get the position of the OperationalHandle.
     /*! \param none
      \returns 3-dimensional real_t precision array (3D position)
     */
    Vector3r getPosition(){ return position; }

    //! Set the position of the OperationalHandle
    /*! \param a 3D point of real_t precision numbers.
      \returns void
     */
    void setPosition(const Vector3r &position){ this->position = position; }

    //! Get the orientation of the OperationalHandle.
     /*! \param none
       \returns a 3x3 real_t precision matrix R in SO(3)
     */
    Matrix3r getOrientation(){ return orientation; }

    //! \brief Set the orientation of the OperationalHandle
     /*! \param a 3x3 real_t precision matrix R in SO(3).
       \returns void
     */
    void setOrientation(const Matrix3r &orientation){ this->orientation = orientation; }
};

} // end of namespace core
} // end of namespace hr

#endif // OPERATIONALHANDLE_H_
