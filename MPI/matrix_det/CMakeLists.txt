cmake_minimum_required(VERSION 3.5.1)
project(matrix_det LANGUAGES CXX C)

# Set the names of the executable
SET(EXE_NAME matdet)

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

file(GLOB SOURCES ./src/*.cpp ./main.cpp)

ADD_EXECUTABLE(${EXE_NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_LIBRARIES})