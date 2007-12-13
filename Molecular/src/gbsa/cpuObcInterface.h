
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

#ifndef __CpuObcInterface_H__
#define __CpuObcInterface_H__

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif

#include "SimTKOpenMMRealType.h"
#include <stdio.h>

/**---------------------------------------------------------------------------------------

   Retrieve the calculated implicit solvation energy from the static class member

   @return the calculated energy from the static class member

   --------------------------------------------------------------------------------------- */

externC float cpuGetImplicitSolventEnergy( void );

/**---------------------------------------------------------------------------------------

   Delete the Obc associated object(s)

   @return 0 if static CpuObc object was set; else return -1

   --------------------------------------------------------------------------------------- */

externC int cpuDeleteObcParameters( void );

/**---------------------------------------------------------------------------------------

	Setup for Obc calculations from Gromacs

   @param numberOfAtoms            number of atoms
   @param obcScaleFactors          array of OBC scale factors (one entry each atom)
   @param atomicRadii              atomic radii in Angstrom (one entry each atom)
   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculations
   @param soluteDielectric         solute dielectric
   @param solventDielectric        solvent dielectric
   @param log                      log reference -- if NULL, then errors/warnings
                                   output to stderr

   The method creates a CpuObc instance -- currently the OBC type II model is the
   default (see paper). If the OBC type I model is desired change

      ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );
   to
      ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeI  );

   The created object is a static member of the class CpuObc; 
   when the force routine, cpuCalculateObcForces(), is called, 
   the static object is used to compute the forces and energy 

   @return 0

   --------------------------------------------------------------------------------------- */

externC 
int cpuSetObcParameters( int numberOfAtoms, RealOpenMM* atomicRadii, RealOpenMM* obcScaleFactors,
                         int includeAceApproximation, RealOpenMM soluteDielectric,
                         RealOpenMM solventDielectric, FILE* log );

/**---------------------------------------------------------------------------------------

   Calculate implicit solvent forces and energy

   @param atomCoordinates   atom coordinates in Angstrom; format of array is
                            atomCoordinates[atom][3]
   @param partialCharges    partial atom charges
   @param forces            output forces in kcal/mol.A; format of array is 
                            forces[atom][3]
   @param energy            energy

   Function calls a static method in CpuImplicitSolvent class to calculate forces/energy

   @return result from CpuImplicitSolvent::computeImplicitSolventForces

   --------------------------------------------------------------------------------------- */

externC int cpuCalculateImplicitSolventForces( RealOpenMM** atomCoordinates,
                                               const RealOpenMM* partialChargesIn,
                                               RealOpenMM** forces, RealOpenMM* energy );

/**---------------------------------------------------------------------------------------

   Get OBC scale factors given masses

   @param numberOfAtoms number of atoms
   @param masses        input masses 
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

externC int getObcScaleFactorsGivenAtomMasses( int numberOfAtoms, const RealOpenMM* masses,
                                               RealOpenMM* scaleFactors );

/**---------------------------------------------------------------------------------------

   Get OBC scale factors given atomic numbers

   @param numberOfAtoms number of atoms
   @param atomicNumber  input atomic number for each atom
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

externC int getObcScaleFactors( int numberOfAtoms, const int* atomicNumber, RealOpenMM* scaleFactors );

/**---------------------------------------------------------------------------------------

   Get GBSA radii

   @param numberOfAtoms             number of atoms
   @param atomicNumber              input atomic number for each atom
   @param numberOfCovalentPartners  input number of covalent partners for each atom
                                    1 for H, ...; only used for C, N & O
   @param indexOfCovalentPartner    index of covalent partner -- only used for H
                                    e.g. if atom 22 is a H and it is bonded to atom 24
                                    then indexOfCovalentPartner[22] = 24
   @param gbsaRadii                 output GBSA radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

externC int getGbsaRadii( int numberOfAtoms, const int* atomicNumber, 
                          const int* numberOfCovalentPartners, 
                          const int* indexOfCovalentPartner, RealOpenMM* gbsaRadii );


#undef externC

#endif
