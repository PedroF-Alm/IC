cmake_minimum_required(VERSION 3.5.1)
project(Doutorado LANGUAGES CXX C)

# Set the names of the executable
SET(EXE_NAME  exec)

# Specify path to SUNDIALS and MY header files
SET(SUNDIALS_INC_DIR
  /usr/local/sundials/include
  CACHE STRING
  "Location of SUNDIALS header files")
SET(SUNDIALS_INC_DIR ${SUNDIALS_INC_DIR} ./header)

# Add path to ALL header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries 
SET(SUNDIALS_LIB_DIR /usr/local/sundials/lib)

# Find the SUNDIALS solver's library
# FIND_LIBRARY(SUNDIALS_SOLVER_LIB
#   sundials_cvode ${SUNDIALS_LIB_DIR}
#   DOC "CVODE library")

# Find the NVECTOR library
# FIND_LIBRARY(SUNDIALS_NVEC_LIB
#   sundials_nvecserial ${SUNDIALS_LIB_DIR}
#   DOC "NVECTOR library")

# IF(SUNDIALS_SOLVER_LIB AND SUNDIALS_NVEC_LIB)
# 	option(USE_CVODE "Using CVODE" ON)
# 	message(STATUS "CVODE Found!")
# ELSE(SUNDIALS_SOLVER_LIB AND SUNDIALS_NVEC_LIB) 
# 	option(USE_CVODE "Not Using CVODE" OFF)
# 	message(STATUS "CVODE NOT Found!")
# ENDIF(SUNDIALS_SOLVER_LIB AND SUNDIALS_NVEC_LIB)

#OpenMP-----------------------
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (OPENMP_OK true)
	option(USE_OPENMP "Usar OpenMP" ON)
	message(STATUS "OpenMP Encontrado!")
else(OPENMP_FOUND) 
	message(STATUS "OpenMP NÃO encontrado!")
	set (OPENMP_OK false)
	option(USE_OPENMP "Usar OpenMP" OFF)
endif(OPENMP_FOUND)
#-----------------------------

#MPI-----------------------
find_package(MPI)
if (MPI_FOUND)
	INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
	LIST(APPEND SCR_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES})
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
	set (MPI_OK true)
	option(USE_MPI "Usar MPI" ON)
	message(STATUS "MPI Encontrado!")
else(MPI_FOUND) 
	message(STATUS "MPI NÃO encontrado!")
	set (MPI_OK false)
	option(USE_MPI "Usar MPI" OFF)
endif(MPI_FOUND)
#-----------------------------
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR} ${MPI_INCLUDE_PATH})

#Random-----------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#-----------------------------

# List of all libraries
# SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB})

file(GLOB SOURCES ./src/*.cpp ./main.cpp)
#configure_file (./header/prototypes_cmake/CVODE.h.in ../header/CVODE.h)
#configure_file (./header/prototypes_cmake/Model.h.in ../header/Model.h)

ADD_EXECUTABLE(${EXE_NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${SUNDIALS_LIBS} ${MPI_LIBRARIES})
