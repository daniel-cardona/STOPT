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
 *	\file include/openhrc/core/world.h
 *	\author Julio Jarquin, Gustavo Arechavaleta
 *	\version 1.0
 *	\date 2015
 *
 *	World class.
 */
#ifndef HR_CORE_WORLD_H
#define HR_CORE_WORLD_H

#include <memory>
#include "multibody.h"

namespace hr{
namespace core{

//! Enum used to specify the robot type. A robot can be NAO, HRP2, YOUBOT.
/*! \ingroup core_module  */
enum RobotType
{
    kNAO = 0, /*!< The multibody is a humanoid robot NAO */
    kHRP2 = 1, /*!< The multibody is a humanoid robot HRP2 */
    kYOUBOT = 2, /*!< The multibody is a YOUBOT */
    kPA10 = 3 /*!< The multibody is a PA10 */
};

//! World class. A world can hold many bodies and multibodies, and test the collision between the bodies.
/*! \ingroup core_module  */
class World
{
    // --------------------------------------------
    // Constructors and Destructor
    // --------------------------------------------
public:
    //! Default constructor.
    World(){}   

    // --------------------------------------------
    // Members
    // --------------------------------------------
private:

    //! Vector of pointers to Multibody.
    /*! This multibodies are the robots in the world. */
    std::vector< std::shared_ptr< MultiBody > > robots;

    //! Vector of pointers to Body.
    /*! This bodies can be considered obstacles in the world. */
    std::vector< Body* > obstacles;

    //! Iterator for the vector of pointers to robots.
    std::vector< std::shared_ptr< MultiBody > >::iterator robots_iter;

    //! Iterator for the vector of pointers to bodies.
    std::vector< Body* >::iterator robotBodies_iter;


    // --------------------------------------------
    // Methods
    // --------------------------------------------
public:
    //! Vector of pointers to Body.
    /*! These bodies can be considered obstacles in the world. */
    /*! \param sFileName The XML file name of the robot geometry specficiation.
        \param robotID The robot identifier.
        \param robotType The robot type (NAO, HRP2, etc.).
        \return void. */
    void loadMultiBody(std::string sFileName, int robotID, const RobotType robotType);

    //! Return a robot.
    /*! \param id The robot identifier.
        \return A pointer to the Multibody with the corresponding name.*/
    std::shared_ptr< MultiBody > getRobot(short int id);

    //! Return a robot.
    /*! \param robotName The name of the robot.
        \return A pointer to the Multibody with the corresponding name.*/
    std::shared_ptr< MultiBody > getRobot(std::string robotName);

    //! Return the vector of robots.
    std::vector< std::shared_ptr< MultiBody > >* getRobotsVector(){return &robots;}

    //! The the collision of the robots and bodies in the world.
    void collisionQuery();

    //! Vector of pointers to Body.
    /*! These bodies can be considered obstacles in the world. */
    /*! \param sFileName The URDF_XML file name of the robot specficiation.
        \param robotID The robot identifier.
        \param robotType The robot type (NAO, HRP2, etc.).
        \return void. */
    void loadMultiBodyURDF(std::string sFileName, int robotID, const RobotType robotType);


    //bool inCollision(Body* body1, Body* body2);
    //bool collisionQueryAllowed(Body* body1, int robot1_ID, Body* body2, int robot2_ID);
    //void partialCollisionQuery(Body* body, int robotID);

};
}   // end of namespace core
}   // end of namesapce hr

#endif // HR_CORE_WORLD_H

