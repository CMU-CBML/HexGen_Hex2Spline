cmake_minimum_required (VERSION 2.8)
project (ModelGen)

# add our path and macros to cmake module path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include (Common)
##############################################################################
# Mathematical libraries
##############################################################################
find_package (Eigen3 REQUIRED 3.0.0) 
if (NOT EIGEN3_INCLUDE_DIR)
  message (FATAL_ERROR "Eigen3 not found!")
endif ()
include_directories(${EIGEN3_INCLUDE_DIR})
add_definitions(${EIGEN3_DEFINITIONS})

#find_package (Matlab)
#if (NOT Matlab_FOUND )
#  message (FATAL_ERROR "Matlab not found!")
#endif ()
#include_directories(${Matlab_INCLUDE_DIR})

#find_package (MKL)
#if (NOT MKL_FOUND )
#  message (FATAL_ERROR "MKL not found!")
#endif ()
#include_directories(${MKL_INCLUDE_DIR})

if (WIN32)
  add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX)
endif ()

 

#find_package (SuiteSparse)
#if (NOT SuiteSparse_FOUND )
#  message (FATAL_ERROR "SuiteSparse not found!")
#endif ()
#include_directories(${SuiteSparse_INCLUDE_DIR})
#
#find_package (GoTools)
#if (NOT GoTools_FOUND )
#  message (FATAL_ERROR "GoTools not found!")
#endif ()
#include_directories(${GoTools_INCLUDE_DIR})

##############################################################################
# Boost
##############################################################################
#SET(Boost_ADDITIONAL_VERSIONS "1.41" "1.47.0"
#  "1.49" "1.54" )
#FIND_PACKAGE(Boost)
#if (NOT Boost_FOUND)
#	message (FATAL_ERROR "Boost not found.")
#endif()

#include_directories(${Boost_INCLUDE_DIRS})
#include_directories(${Boost_INCLUDE_LIB_DIRS})
#add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})


##############################################################################
# for the compatibility of GoTools
# Check if compiler supports c++-0x
#INCLUDE(CheckCXXCompilerFlag)
#IF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
#  CHECK_CXX_COMPILER_FLAG("-std=gnu++0x" HAVE_0x)
#  IF(HAVE_0x)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
#  ELSE(HAVE_0x)
#    # C++0x is not supported - check for Boost?
#    IF(Boost_FOUND)
#      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
#      SET(GoTools_COMMON_INCLUDE_DIRS
#        ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#    ELSE(Boost_FOUND)
#      MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
#    ENDIF(Boost_FOUND)
#  ENDIF(HAVE_0x)
#ENDIF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
#IF(CMAKE_CXX_COMPILER_ID MATCHES Intel)
#  # icpc's c++0x is lacking
#  IF(Boost_FOUND)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
#    SET(GoTools_COMMON_INCLUDE_DIRS
#      ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#  ELSE(Boost_FOUND)
#    MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
#  ENDIF(Boost_FOUND)
#ENDIF(CMAKE_CXX_COMPILER_ID MATCHES Intel)
#IF(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
#  IF(MSVC90) # Visual Studio 2008 support (c++0x is missing)
#    IF(Boost_FOUND)
#      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=1")
#      SET(GoTools_COMMON_INCLUDE_DIRS
#      ${GoTools_COMMON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
#    ELSE(Boost_FOUND)
#      MESSAGE(FATAL_ERROR "Either Boost or a compiler with c++0x support is needed")
#    ENDIF(Boost_FOUND)
#  ENDIF(MSVC90)
#ENDIF(CMAKE_CXX_COMPILER_ID MATCHES MSVC)


#add_subdirectory (trilib)
#INCLUDE_DIRECTORIES(trilib)
##############################################################################
# OpenMP
##############################################################################
# disable openmp if using with intel mkl_rt.lib. 
# if linking to four mkl layer libs, then openmp is ok.
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	add_definitions(-DUSE_OPENMP)
else ()
	set (OPENMP_NOTFOUND 1)
endif ()
 

##############################################################################
# CoMISo
##############################################################################
#set (CoMISo_DIR 
#  ${CMAKE_SOURCE_DIR}/CoMISo/Solver 
#  ${CMAKE_SOURCE_DIR}/CoMISo/NSolver 
#  ${CMAKE_SOURCE_DIR}/CoMISo/EigenSolver
#  ${CMAKE_SOURCE_DIR}/CoMISo/Config
#  ${CMAKE_SOURCE_DIR}/CoMISo/Utils
#)
#gather_list_by_extension (comiso_Tcc "*T.cc" ${CoMISo_DIR})
#gather_list_by_extension (comiso_cc "*.cc" ${CoMISo_DIR})
#foreach (_ccfile ${comiso_cc})
#			get_filename_component (_filename ${_ccfile} NAME)
#			# message(STATUS "file " ${_filename})
#			if (_filename MATCHES "\\T.cc$$")
#				# message(STATUS "removed")
#				list (REMOVE_ITEM comiso_cc ${_ccfile})
#			endif ()
#endforeach ()
#gather_list_by_extension (comiso_hh "*.hh" ${CoMISo_DIR})
#include_directories(${CMAKE_SOURCE_DIR})
#add_definitions(-DINCLUDE_TEMPLATES) # Tell CoMISo to include *T.hpp into *.h


macro(source_group_by_dir source_files)
 if(MSVC)
  set(sgbd_cur_dir ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(sgbd_file ${${source_files}})
   string(REGEX REPLACE ${sgbd_cur_dir}/\(.*\) \\1 sgbd_fpath ${sgbd_file})
   string(REGEX REPLACE "\(.*\)/.*" \\1 sgbd_group_name ${sgbd_fpath})
   string(COMPARE EQUAL ${sgbd_fpath} ${sgbd_group_name} sgbd_nogroup)
   string(REPLACE "/" "\\" sgbd_group_name ${sgbd_group_name})
   if(sgbd_nogroup)
    set(sgbd_group_name "\\")
   endif(sgbd_nogroup)
   source_group(${sgbd_group_name} FILES ${sgbd_file})
  endforeach(sgbd_file)
 endif(MSVC)
endmacro(source_group_by_dir)


AUX_SOURCE_DIRECTORY(cpp SOURCES)
#AUX_SOURCE_DIRECTORY(MeshIO SOURCES)
#AUX_SOURCE_DIRECTORY(Constructor SOURCES)
#AUX_SOURCE_DIRECTORY(MeshStructure SOURCES)
#AUX_SOURCE_DIRECTORY(SparseSolver SOURCES)
#AUX_SOURCE_DIRECTORY(NURBS SOURCES)
#AUX_SOURCE_DIRECTORY(ConstructInclude SOURCES)
#AUX_SOURCE_DIRECTORY(Wap SOURCES)

SET(SOURCES ${SOURCES})
#FILE(GLOB_RECURSE project_headers cpp/*.h MeshStructure/*.h SparseSolver/*.h ANN/*.h MeshIO/*.h Analysis/*.h Wap/*.h* ConstructInclude/*.h NURBS/*.h Constructor/*.h)
#FILE(GLOB_RECURSE project_cpps cpp/*.c* MeshStructure/*.cpp SparseSolver/*.cpp ANN/*.cpp MeshIO/*.cpp Analysis/*.cpp Wap/*.cpp ConstructInclude/*.c* NURBS/*.c* Constructor/*.cpp)
FILE(GLOB_RECURSE project_headers cpp/*.h)
FILE(GLOB_RECURSE project_cpps cpp/*.c*)
set(all_files ${project_headers} ${project_cpps})
source_group_by_dir(all_files)

###################################################################################################
# Add special directory of files which cannot be generally added by the previous method
###################################################################################################
FILE (GLOB PCHFiles "${CMAKE_SOURCE_DIR}/stdafx.*")
source_group ("Precompiled Headers" FILES ${PCHFiles})
#source_group ("CoMISo Files" FILES ${comiso_cc} ${comiso_hh} )
	
##############################################################################
# Application specific preprocessor definitions
##############################################################################
IF(MSVC)
add_definitions(-wd4018 -wd4267 -wd4101 -wd4190 -wd4244 -wd4305 -wd4099) # Turn off warnings
ENDIF(MSVC)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
#INCLUDE_DIRECTORIES(Entrance Constructor MeshStructure SparseSolver ANN CoMISo NURBS Analysis ConstructInclude Wap MeshIO)
ADD_EXECUTABLE(ModelGen ${all_files})
#target_link_libraries(TTSP_3D ${SuiteSparse_LIBS} ${GoTools_LIBS} ${Boost_LIBRARIES} ${Matlab_LIBS} ) 
#target_link_libraries(conv3D ${Matlab_LIBS} ) 
IF(MSVC)
SET_TARGET_PROPERTIES(ModelGen PROPERTIES COMPILE_FLAGS "/MP")
ELSE(MSVC)
add_definitions(-boost_system -boost_filesystem)
ENDIF(MSVC)
 

