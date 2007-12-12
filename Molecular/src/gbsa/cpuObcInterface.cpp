
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

#include "cpuObcInterface.h"

#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
// #include "SimTKOpenMMGromacsUtilities.h"
#include "ObcParameters.h"
#include "CpuObc.h"

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

extern "C" int 
cpuSetObcParameters( int numberOfAtoms, Real* atomicRadii, Real* obcScaleFactors,
                     int includeAceApproximation,
                     Real soluteDielectric, Real solventDielectric, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ncpuSetObcParameters: ";

   // ---------------------------------------------------------------------------------------
   
   // set log file if not NULL

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   // set OBC parameters (Type II)

   ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );
   obcParameters->setScaledRadiusFactors( obcScaleFactors );
   obcParameters->setAtomicRadii( atomicRadii, SimTKOpenMMCommon::KcalAngUnits );

   // dielectric constants

   obcParameters->setSolventDielectric( solventDielectric );
   obcParameters->setSoluteDielectric( soluteDielectric );

   // ---------------------------------------------------------------------------------------

   // create CpuObc instance that will calculate forces
  
   CpuObc* cpuObc = new CpuObc( obcParameters );

   // set static member for subsequent calls to calculate forces/energy 

   CpuImplicitSolvent::setCpuImplicitSolvent( cpuObc );

   // set base file name, ...

   //cpuObc->readInfoFile( "CpuImplicitSolventInfo" );

   // include/do not include ACE approximation (nonpolar solvation)

   cpuObc->setIncludeAceApproximation( includeAceApproximation );

   // ---------------------------------------------------------------------------------------

   // diagnostics
 
   if( log ){
      std::string state = cpuObc->getStateString( methodName );
      (void) fprintf( log, "\n%s\nDone w/ setup\n", state.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Calculate implicit solvent forces and energy

   @param atomCoordinates   atom coordinates in Angstrom; format of array is
                            atomCoordinates[atom][3] in Angstrom

   @param partialCharges    Gromacs charges ('mdatoms->chargeA' in md.c)

   @param forces            output forces in kcal/mol.A; format of array is 
                            forces[atom][3]

   @param energy            energy in kcal/mol

   Function calls a static method in CpuImplicitSolvent class to calculate forces/energy

   @return result from CpuImplicitSolvent::computeImplicitSolventForces

   --------------------------------------------------------------------------------------- */

extern "C" int
cpuCalculateImplicitSolventForces( Real** atomCoordinates, const Real* partialCharges,
                                   Real** forces, Real* energy ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName   = "\ncpuCalculateImplicitSolventForces: ";

   // ---------------------------------------------------------------------------------------

   int status = CpuImplicitSolvent::computeImplicitSolventForces( atomCoordinates, partialCharges,
                                                                  forces );

   *energy = CpuImplicitSolvent::getCpuImplicitSolvent()->getEnergy(); 
   // printf( "\ncpuCalculateImplicitSolventForcesE=%.5e", *energy );

   return status;

}

/**---------------------------------------------------------------------------------------

   Retrieve the calculated energy from the static class member

   @return the calculated energy from the static class member

   --------------------------------------------------------------------------------------- */

extern "C" float cpuGetImplicitSolventEnergy( void ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName   = "\ncpuGetImplicitSolventEnergy: ";

   // ---------------------------------------------------------------------------------------

   float energy =  CpuImplicitSolvent::getCpuImplicitSolvent()->getEnergy();
   // printf( "\ncpuGetImplicitSolventEnergy E=%.5e", energy );

   return energy;
}

/**---------------------------------------------------------------------------------------

   Delete the Obc associated object(s)

   @return 0 if static CpuObc object was set; else return -1

   --------------------------------------------------------------------------------------- */

extern "C" int cpuDeleteObcParameters( void ){
   return CpuImplicitSolvent::deleteCpuImplicitSolvent();
}

/**---------------------------------------------------------------------------------------

   Get OBC scale factors given masses

   @param numberOfAtoms number of atoms
   @param masses        input masses 
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getObcScaleFactorsGivenAtomMasses( int numberOfAtoms, const Real* masses, Real* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetObcScaleFactorsGivenAtomMasses";

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      double scaleFactor = 0.8;
      Real mass          = masses[atomI];

      if ( mass < 1.2 && mass >= 1.0 ){        // hydrogen
         scaleFactor  = 0.85; 
      } else if( mass > 11.8 && mass < 12.2 ){ // carbon
         scaleFactor  = 0.72; 
      } else if( mass > 14.0 && mass < 15.0 ){ // nitrogen
         scaleFactor  = 0.79;
      } else if( mass > 15.5 && mass < 16.5 ){ // oxygen
         scaleFactor  = 0.85; 
      } else if( mass > 31.5 && mass < 32.5 ){ // sulphur
         scaleFactor  = 0.96;
      } else if( mass > 29.5 && mass < 30.5 ){ // phosphorus
         scaleFactor  = 0.86;
      } else {
         std::stringstream message;
         message << methodName;
         message << " Warning: mass for atom " << atomI << " mass=" << mass << "> not recognized.";
         SimTKOpenMMLog::printMessage( message );
      }

      scaleFactors[atomI] = (Real) scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get OBC scale factors given atomic numbers

   @param numberOfAtoms number of atoms
   @param atomicNumber  input atomic number for each atom
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getObcScaleFactors( int numberOfAtoms, const int* atomicNumber, Real* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetObcScaleFactors";

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      double scaleFactor;
      switch( atomicNumber[atomI] ){

         case 1: // hydrogen

            scaleFactor  = 0.85; 
            break;

         case 6: // carbon

            scaleFactor  = 0.72; 
            break;

         case 7: // nitrogen

            scaleFactor  = 0.79;
            break;

         case 8: // oxygen

            scaleFactor  = 0.85;
            break;

         case 15: // phosphorus

            scaleFactor  = 0.86;
            break;

         case 16: // sulphur

            scaleFactor  = 0.85;
            break;

         default:

            scaleFactor = 0.8;

            std::stringstream message;
            message << methodName;
            message << " Warning: atom number=" << atomicNumber[atomI] << " for atom " << atomI << "> not handled -- useing default value.";
            SimTKOpenMMLog::printMessage( message );
            break;
      }

      scaleFactors[atomI] = (Real) scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get GBSA radii

   @param numberOfAtoms             number of atoms

   @param atomicNumber              input atomic number for each atom

   @param numberOfCovalentPartners  input number of covalent partners for each atom
                                    1 for H, 2,3,4 for C, ...; 
                                    the values are only used for C, N & O

   @param indexOfCovalentPartner    index of covalent partner -- used only for H
                                    e.g. if atom 22 is a H and it is bonded to atom 24,
                                    then indexOfCovalentPartner[22] = 24

   @param gbsaRadii                 output GBSA radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getGbsaRadii( int numberOfAtoms, const int* atomicNumber, 
                             const int* numberOfCovalentPartners, 
                             const int* indexOfCovalentPartner,
                             Real* gbsaRadii ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetGbsaRadii";

   // ---------------------------------------------------------------------------------------

   // loop over atoms

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // branch based on atomic number

      double radius;
      int atomIndex;
      switch ( atomicNumber[atomI] ){

         case 1: // H

            // get index of covalent partner and validate
            // radius is modified if heavy atom is N or O

            atomIndex = indexOfCovalentPartner[atomI];
            if( atomIndex < 0 || atomIndex >= numberOfAtoms ){
               std::stringstream message;
               message << methodName;
               message << " Warning: atomIndex for covalent partner=" << atomIndex << " for atom " << atomI << " is invalid.";
               SimTKOpenMMLog::printMessage( message );
            } else if( atomicNumber[atomIndex] == 7 ){
               radius = 1.15;
            } else if( atomicNumber[atomIndex] == 8 ){
               radius = 1.05;
            } else {
               radius = 1.25;
            }
            break;

         case 3: // Li

            radius = 2.432;
            break;

         case 6: // C

            if( numberOfCovalentPartners[atomI] == 2 ){
               radius = 1.825;
            } else if( numberOfCovalentPartners[atomI] == 3 ){
               radius = 1.875;
            } else {
               radius = 1.90;
            }
            break;

         case 7: // N

            if( numberOfCovalentPartners[atomI] == 4 ){
               radius = 1.625;
            } else if( numberOfCovalentPartners[atomI] == 1 ){
               radius = 1.60;
            } else {
               radius = 1.7063;
            }
            break;

         case 8: // O

            if( numberOfCovalentPartners[atomI] == 1 ){
               radius = 1.48;
            } else {
               radius = 1.535;
            }
            break;

         case 9:
            radius = 1.47;
            break;
         case 10:
            radius = 1.39;
            break;
         case 11:
            radius = 1.992;
            break;
         case 12:
            radius = 1.70;
            break;
         case 14:
            radius = 1.80;
            break;
         case 15:
            radius = 1.87;
            break;
         case 16:
            radius = 1.775;
            break;
         case 17:
            radius = 1.735;
            break;
         case 18:
            radius = 1.70;
            break;
         case 19:
            radius = 2.123;
            break;
         case 20:
            radius = 1.817;
            break;
         case 35:
            radius = 1.90;
            break;
         case 36:
            radius = 1.812;
            break;
         case 37:
            radius = 2.26;
            break;
         case 53:
            radius = 2.10;
            break;
         case 54:
            radius = 1.967;
            break;
         case 55:
            radius = 2.507;
            break;
         case 56:
            radius = 2.188;
            break;
         default:
            radius = 2.0;
            break;
      }
         
      gbsaRadii[atomI] = (Real) radius;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}
