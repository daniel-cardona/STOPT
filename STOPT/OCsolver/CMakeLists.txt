cmake_minimum_required(VERSION 2.8.9)
project(Sparse_Trapezoidal_Optimization)

####################### SET MAIN DEPENDENCES ###########################

#--------INCLUDE THE MAIN HEADERS

include_directories(include)

#--------EIGEN

include_directories("../eigen3")

#-------IPOPT

list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../ipopt/cmake/linux-64)
find_package(IPOPT)
include_directories(${IPOPT_INCLUDE_DIRS})


#-------PQP

list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../PQP/cmake/linux-64)
find_package(PQP)
include_directories(${PQP_INCLUDE_DIRS})

#-------OPEN HRC

list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../openhrc/cmake/linux-64)
find_package(openhrc)
include_directories(${OpenHRC_INCLUDE_DIRS})

##################### CONFIGURATION OF BUILD TYPE ####################

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)


#####################################################################
###################### Configure additional packages ################
#####################################################################

#However, the file(GLOB...) allows for wildcard additions:

file(GLOB SOLVER 	"src/RobOptTraj/directTranscription/core/*.cpp")
file(GLOB NLP 	 	"src/RobOptTraj/directTranscription/nlp/*.cpp")
file(GLOB NLPSOLVER 	"src/RobOptTraj/directTranscription/nlpSolver/*.cpp")
file(GLOB COLMETHOD 	"src/RobOptTraj/directTranscription/collocationMethods/trapezoidal/*.cpp")
file(GLOB DERIVATIVES 	"src/RobOptTraj/directTranscription/derivatives/*.cpp")
file(GLOB UTILS 	"src/RobOptTraj/directTranscription/utils/*.cpp")


########################################################################
######################## Compile the test###############################
########################################################################

file(GLOB obstacle "example/test/*.cpp")
add_executable(test ${obstacle} ${SOLVER} ${NLP} ${NLPSOLVER} ${COLMETHOD} ${DERIVATIVES} ${UTILS})
target_link_libraries(test ${IPOPT_LIBRARIES} ${OpenHRC_LIBRARIES} ${PQP_LIBRARIES})

########################################################################
######################## Compile the examples ##########################
########################################################################


file(GLOB 5_dof "example/5DOF_Robot/*.cpp")
add_executable(fingerRobot ${5_dof} ${SOLVER} ${NLP} ${NLPSOLVER} ${COLMETHOD} ${DERIVATIVES} ${UTILS})
set_target_properties(fingerRobot PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../example/5DOF_Robot")
target_link_libraries(fingerRobot ${IPOPT_LIBRARIES} ${OpenHRC_LIBRARIES} ${PQP_LIBRARIES})


file(GLOB Kuka "example/Kuka_Youbot/*.cpp")
add_executable(kukaYoubot ${Kuka} ${SOLVER} ${NLP} ${NLPSOLVER} ${COLMETHOD} ${DERIVATIVES} ${UTILS})
set_target_properties(kukaYoubot PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../example/Kuka_Youbot")
target_link_libraries(kukaYoubot ${IPOPT_LIBRARIES} ${OpenHRC_LIBRARIES} ${PQP_LIBRARIES})


file(GLOB nao "example/NAO_Robot/*.cpp")
add_executable(naoRobot ${nao} ${SOLVER} ${NLP} ${NLPSOLVER} ${COLMETHOD} ${DERIVATIVES} ${UTILS})
set_target_properties(naoRobot PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../example/NAO_Robot")
target_link_libraries(naoRobot ${IPOPT_LIBRARIES} ${OpenHRC_LIBRARIES} ${PQP_LIBRARIES})

