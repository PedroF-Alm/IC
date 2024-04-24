cmake_minimum_required(VERSION 3.5.1)
project(matrix_multiply LANGUAGES CXX C)

# Set the names of the executable
SET(EXE_NAME mamult)

#OpenMP-----------------------
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (OPENMP_OK true)
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	option(USE_OPENMP "Usar OpenMP" ON)
	message(STATUS "OpenMP Encontrado!")
else(OPENMP_FOUND) 
	message(STATUS "OpenMP NÃO encontrado!")
	set (OPENMP_OK false)
	option(USE_OPENMP "Usar OpenMP" OFF)
endif(OPENMP_FOUND)
#-----------------------------

#Random-----------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#-----------------------------

file(GLOB SOURCES ./src/*.cpp ./main.cpp)

ADD_EXECUTABLE(${EXE_NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_LIBRARIES})