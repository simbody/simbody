/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/08 00:48:24 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 *------------------------------------------------------------------
 * SUNDIALS configuration header file
 *------------------------------------------------------------------
 */
#ifndef SimTK_SIMMATH_SUNDIALS_CONFIG_H_
#define SimTK_SIMMATH_SUNDIALS_CONFIG_H_ 

/***** ADDED BY SHERM 20061121 */
#include "SimTKcommon.h"

/* Keeps MS VC++ 8 quiet about sprintf, strcpy, etc. (sherm) */
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

// No need to expose any of this; SimTK users must access CPodes
// through the C++ API.
#if defined(_WIN32) && defined (_MSC_VER)
    #if defined(SimTK_SIMMATH_BUILDING_SHARED_LIBRARY)
        #define SUNDIALS_EXPORT __declspec(dllexport)
    #elif defined(SimTK_SIMMATH_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SUNDIALS_EXPORT
    #else
        /* i.e., a client of a shared library */
        #define SUNDIALS_EXPORT __declspec(dllimport)
    #endif
#else
    /* Linux, Mac */
    #define SUNDIALS_EXPORT
#endif

/* Define precision of SUNDIALS data type 'realtype' 
 * Depending on the precision level, one of the following 
 * three macros will be defined:
 *     #define SUNDIALS_SINGLE_PRECISION 1
 *     #define SUNDIALS_DOUBLE_PRECISION 1
 *     #define SUNDIALS_EXTENDED_PRECISION 1
 */
#ifdef SimTK_DEFAULT_PRECISION
    #if   (SimTK_DEFAULT_PRECISION == 1)
        #define SUNDIALS_SINGLE_PRECISION 1
    #elif (SimTK_DEFAULT_PRECISION == 2)
        #define SUNDIALS_DOUBLE_PRECISION 1
    #elif (SimTK_DEFAULT_PRECISION == 4)
        #define SUNDIALS_EXTENDED_PRECISION 1
    #endif
#else
    #define SUNDIALS_DOUBLE_PRECISION 1
#endif

/***** END OF SHERM'S ADDITIONS */

/* Define SUNDIALS version number */
#define SUNDIALS_PACKAGE_VERSION "2.3.0"

/* FCMIX: Define Fortran name-mangling macro 
 * Depending on the inferred scheme, one of the following 
 * six macros will be defined:
 *     #define F77_FUNC(name,NAME) name
 *     #define F77_FUNC(name,NAME) name ## _
 *     #define F77_FUNC(name,NAME) name ## __
 *     #define F77_FUNC(name,NAME) NAME
 *     #define F77_FUNC(name,NAME) NAME ## _
 *     #define F77_FUNC(name,NAME) NAME ## __
 */
#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _


/* Use generic math functions 
 * If it was decided that generic math functions can be used, then
 *     #define SUNDIALS_USE_GENERIC_MATH 1
 * otherwise
 *     #define SUNDIALS_USE_GENERIC_MATH 0
 */
#define SUNDIALS_USE_GENERIC_MATH 1

/* FNVECTOR: Allow user to specify different MPI communicator
 * If it was found that the MPI implementation supports MPI_Comm_f2c, then
 *      #define SUNDIALS_MPI_COMM_F2C 1
 * otherwise
 *      #define SUNDIALS_MPI_COMM_F2C 0
 */

#endif // SimTK_SIMMATH_SUNDIALS_CONFIG_H_
