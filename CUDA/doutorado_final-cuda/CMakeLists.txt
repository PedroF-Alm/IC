cmake_minimum_required(VERSION 3.5.1)
project(Doutorado LANGUAGES CUDA CXX)

# Set the names of the executable
SET(EXE_NAME  exec)

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
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

#CUDA-----------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK3 REQUIRED gtk+-3.0)
include_directories(${GTK3_INCLUDE_DIRS})
link_directories(${GTK3_LIBRARY_DIRS})

add_definitions(${GTK3_CFLAGS_OTHER})

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 75)
endif()

#Random-----------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON
set(CMAKE_CUDA_STANDARD 14))
#-----------------------------

file(GLOB SOURCES ./src/*.cu ./main.cu)

ADD_EXECUTABLE(${EXE_NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_LIBRARIES} ${GTK3_LIBRARIES})
