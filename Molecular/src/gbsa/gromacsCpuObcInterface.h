
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __GromacsCpuObcInterface_H__
#define __GromacsCpuObcInterface_H__

// Include gromacs types

#include "SimTKOpenMMGromacsUtilities.h"

/*
#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif
*/

/**---------------------------------------------------------------------------------------

   Setup for Obc calculations

   @param top                      Gromacs t_topology (as in md.c)
   @param log                      log reference (stdlog in md.c)
   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculataions
   @param soluteDielectric         solute dielectric
   @param solventDielectric        solvent dielectric

   Method creates a CpuObc instance w/ OBC parameters set

   The created object is a static member of the
   class CpuObc that is referenced to calculate the  implicit solvation
   forces and energy 

  @return 0

   --------------------------------------------------------------------------------------- */

extern "C" 
int gromacsCpuInitialSetup( const t_topology* top, FILE* log, int includeAceApproximation,
                            float soluteDielectric, float solventDielectric );

/**---------------------------------------------------------------------------------------

   Calculate Obc forces and energy

   @param atomCoordinates   Gromacs atom coordinates ('x' in md.c)
   @param partialCharges    Gromacs charges ('mdatoms->chargeA' in md.c)
   @param forces            output forces in kJ/mol.A; the computed forces
                            are added to the entries in the array

   Function calls a static method in CpuObc class to calculate forces/energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsCpuCalculateObcForces( const rvec* atomCoordinates, const float* partialChargesIn, 
                                  rvec* forces );

/**---------------------------------------------------------------------------------------

   Write Tinker Xyz file \n

   @param atomCoordinates   Gromacs atom coordinates ('x' in md.c)
   @param top               Gromacs t_topology (as in md.c)
   @param inputHeader       xyz file header
   @param inputXyzFileName  xyz file name
   @param log               log file reference (may be NULL)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C" 
int cpuWriteXyzFile( const rvec *atomCoordinates, const t_topology* top,
                     const char* inputHeader, const char* inputXyzFileName, FILE* log );

#endif
