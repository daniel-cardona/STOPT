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
 *	\file include/openhrc/core/body.h
 *	\author Gerardo Jarquin, Gustavo Arechavaleta, Julio Jarquin.
 *	\version 1.0
 *	\date 2015
 *
 * Body class. It has functions to build a body from a xml file.
 */

#ifndef HR_CORE_BODY_H
#define HR_CORE_BODY_H

#include <list>
#include "operational_handle.h"
#include "geometry.h"
//#include "PQP/PQP.h"


class PQP_Model; // Forward declaration of the class PQP_Model.

typedef double PQP_REAL; // Forward declaration of the PQP_REAL type used by the PQP library

namespace hr {
namespace core {


//!  Body class definition
/*!  \ingroup core_module
        The Body class is used as a representation of a link in a robot. \sa Joint, OperationalHandle
    */
class Body {
    //! MultiBody class.
    friend class MultiBody;
    //! MultiBodyNAO class.
    friend class MultiBodyNAO;

    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Custom constructor. Initializes the body with the given parameters.
    /*! \param id The id of the body.
        \param name The name of the body.*/
    Body(const unsigned int id, const std::string &name = "End-effector");

    //! Custom constructor. Initializes the body with the given parameters.
    /*! \param id The id of the body.*/
    /*! \param name The name of the body.*/
    /*! \param pos The position of the body.*/
    /*! \param rot The rotation matrix (Rx, Ry, Rz).*/
    Body(const unsigned int id, const std::string &name,
         const Vector3r &pos, const Matrix3r &rot);

    //! Default destructor.
    ~Body();

    // --------------------------------------------
    // Members
    // --------------------------------------------
protected:
    //! Attach an OperationHandle to the body for defining its reference frame.
    OperationalHandle frame;

    //! Body geometric shape.
    Geometry shape; // TODO 03/2015: use a pointer and use forward declaration for the Geometry class. Put a compiler flag to determine if PQP should be used or not

    //! The list of pointers to other operational handles attached to the body.
    std::list<OperationalHandle*> opHandles;

    //! Iterator for the list of operational handles.
    std::list<OperationalHandle*>::iterator opHandles_iter;

    // X_T(i) the position of the frame of handle i relative to the frame of the body
    //!  Relative body position vector.
    std::vector<Vector3r> relativeHandlePosition;

    //! PQP model asociated to the body.
    PQP_Model* pqpBodyModel;   //TODO 03/2015: use a smart pointer.

    //! PQP rotation of the body relative to the inertial frame.
    PQP_REAL pqpRotation[3][3];

    //! PQP position of the body relative to the inertial frame.
    PQP_REAL pqpTranslation[3];



    /* The values of the position and orientation members inherited
              from perationalHandle class are relative to the inertial frame */

    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:

    //!  Set body ID
    /*! \param short int identifying the body
     * \returns void
     */
    void setId(const int &id){ frame.setId(id); }


    //!  Get Body ID
    /*! \param none
     * \returns short int identifying the body.
     */
    short int getId(){ return frame.getId(); }

    //!  Set body name
    /*! \param string with the body name.
    * \returns void
    */
    void setName(const std::string &name){ frame.setName(name); }


    //!  Get body name
    /*! \param name
    * \returns string with the body name.
    */
    std::string getName(){ return frame.getName(); }


    //! Set the geometry of the body
    /*! \param a pointer to an object Geometry
             * \return void
             */
    void setGeometry(const Geometry &shape){ this->shape = shape; }

    //! Get the geometry of the body.
    /*! \param none
             * \return a pointer to an object Geometry
             */
    Geometry* getGeometry(){ return &shape; }

    //! Load the geometry of the body
    /*! \param a string with the name of the xml file
             * \return bool status
             */
    bool loadGeometry(const std::string &xmlFileName, Vector3r &relativeBodyPosition);

    //! Set the position of the Body
    /*! \param a 3D point of real_t precision numbers
    * \return void
    */
    void setPosition(const Vector3r &position );

    //!  Get the position of the body frame.
    /*! \param none
     * \returns 3-dimensional real_t precision array (3D position)
     */
    Vector3r getPosition(){ return frame.getPosition(); }

    //! Set the orientation of the Body
    /*! \param a 3x3 real_t precision matrix R in SO(3).
             * \return void
             */
    void setOrientation(const Matrix3r &orientation);

    //!  Get the orientation of the body frame.
    /*! \param none
     * \returns a 3x3 real_t precision matrix R in SO(3)
     */
    Matrix3r getOrientation(){ return frame.getOrientation(); }



    //! Add an opHandle to the list of pointers of opHandles associated to the body
    /*! \param The 3D position of the opHandle to be added
             * \return void
             */
    int addOperationalHandle(const Vector3r &position);



    //! Remove an opHandle to the list of opHandles associated to the body
    /*! \param the id of the opHandle to be removed
             * \return void
             */
    void removeOperationalHandle(short int id);

    //! Find and get an opHandle associated to the body
    /*! \param opHandle id
             * \return the pointer to the opHandle
             */
    OperationalHandle* getOperationalHandle(short int id);

    //! Update the position of the list of opHandles associated to the body
    /*! \param none
             * \return void
             */
    void updateOperationalHandlesPose();

    //! Get the position of the body i relative to its parent.
    /*! \param the id of the body
         * \return the relative position vector
         */
    Vector3r getRelativeHandlePosition(short int handleId)
    {
        return relativeHandlePosition.at(handleId);
    }


    //! true if there are opHandles associated to the body
    /*! \param none
             * \return true if the list of opHandles is not empty
             */
    bool hasOperationalHandles();

    //! Get the PQP model associated to the body
    /*! \param none
             * \return a pointer to the PQP model of the body
             */
    PQP_Model* getPQPBodyModel(){ return pqpBodyModel; }

    //! Create the PQP model associated to the body
    /*! \param none
             * \return void
             */
    void createPQPBodyModel();






private:

    //! Get the orientation of the Body as PQP object
    /*! \param none.
              \return a PQP_real rotation matrix
             */
    PQP_REAL (*getPQProtation())[3]{ return pqpRotation; }

    //! Get the position of the Body as PQP object
    /*! \param none.
              \return a PQP_real translation vector
             */
    PQP_REAL* getPQPtranslation(){ return pqpTranslation; }

};
} // end of namespace core
} // end of namespace hr

#endif // HR_CORE_BODY_H
