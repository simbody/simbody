
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

#include "gromacsCpuObcInterface.h"
#include "cpuObcInterface.h"

#include "CpuObc.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "SimTKOpenMMRealType.h"

/**---------------------------------------------------------------------------------------

	Setup for Obc calculations from Gromacs

   @param top                      Gromacs t_topology data struct
   @param log                      log reference (stdlog in md.c)
   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculations
   @param soluteDielectric         solute dielectric
   @param solventDielectric        solvent dielectric

   The method creates a CpuObc instance

   The created object is a static member of the
   class CpuObc; when the force routine, cpuCalculateObcForces, is called, 
   the object is used to compute the forces and energy 

   The arrays atomicRadii & obcScaleFactors are allocated off heap
 
   @return 0

   --------------------------------------------------------------------------------------- */

extern "C"
int gromacsCpuInitialSetup( const t_topology* top, FILE* log, int includeAceApproximation,
                            float soluteDielectric, float solventDielectric ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ngromacsCpuInitialSetup: ";
   static int debug                = true;

   // ---------------------------------------------------------------------------------------

   // set log file, if available

/*
fprintf( log, "%s %.3f %.3f %d", methodName, soluteDielectric, solventDielectric, includeAceApproximation );
fflush( log );
if( soluteDielectric < 1.0 ){
   soluteDielectric = 1.0f;
}
if( solventDielectric < 50.0 ){
   solventDielectric = 78.3f;
} */

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   // OBC scale factors based on atomic mass

   int numberOfAtoms             = top->atoms.nr;
   RealOpenMM* obcScaleFactors   = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
   SimTKOpenMMGromacsUtilities::getObcScaleFactors( top, obcScaleFactors );

   // get/set atomic radii (from file)

   RealOpenMM* atomicRadii       = (RealOpenMM*) malloc( sizeof( RealOpenMM )*numberOfAtoms );
   std::string parameterFileName = ObcParameters::ParameterFileName;
   int returnStatus = SimTKOpenMMGromacsUtilities::getMacroModelAtomicRadii( numberOfAtoms,
                                                                             parameterFileName,
                                                                             top, atomicRadii );

   // abort if problem reading parameter file

   if( returnStatus != SimTKOpenMMCommon::DefaultReturn ){
      (void) fprintf( (log ? log : stderr), "%s problem getting atomic radii from file=<%s>",
                      methodName, parameterFileName.c_str() );
      (void) fflush( log );
      exit(-1);
   } else if( log ){
      (void) fprintf( log, "%s obtained atomic radii from file=<%s>",
                      methodName, parameterFileName.c_str() );
      (void) fflush( log );
   }

   int status = cpuSetObcParameters( numberOfAtoms, atomicRadii, obcScaleFactors,
                                     includeAceApproximation, soluteDielectric,
                                     solventDielectric, log );

   // set flags so that atomicRadii & scaleFactor arrays are deleted
   // if parameter array is deleted

   CpuImplicitSolvent* cpuImplicitSolvent = CpuImplicitSolvent::getCpuImplicitSolvent();
   CpuObc* cpuObc                         = static_cast<CpuObc*> (cpuImplicitSolvent);
   ObcParameters* obcParameters           = cpuObc->getObcParameters();
   obcParameters->setOwnAtomicRadii( true );
   obcParameters->setOwnScaleFactors( true );

   return status;
}

/**---------------------------------------------------------------------------------------

   Calculate Obc forces and energy

   @param atomCoordinates   Gromacs atom coordinates ('x' in md.c)
   @param partialCharges    Gromacs charges ('mdatoms->chargeA' in md.c)
   @param forces            output forces in kJ/mol.A; the computed forces
                            are added to the entries in the array

   Function calls a static method in CpuObc class to calculate forces/energy

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C" int gromacsCpuCalculateObcForces( const rvec* atomCoordinates,
                                             const RealOpenMM* partialChargesIn, 
                                             rvec* forces ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName         = "\ngromacsCpuCalculateObcForces: ";
   static RealOpenMM* partialCharges     = NULL;
   static RealOpenMM** coordinates       = NULL;
   static RealOpenMM** localForces       = NULL;

   static const int updateBornRadii      = 0;

   // ---------------------------------------------------------------------------------------

   CpuImplicitSolvent* cpuImplicitSolvent = CpuImplicitSolvent::getCpuImplicitSolvent();
   if( cpuImplicitSolvent == NULL ){
      (void) fprintf( stderr, "%s implicit solvent not initialized -- aborting" );
      (void) fflush( stderr );
      exit(-1);
   }

   int numberOfAtoms = cpuImplicitSolvent->getNumberOfAtoms();

   // first time through, allocate memory for coordinates, ..
   // Note: memory never freed!
 
   if( coordinates == NULL ){

      coordinates    = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );
      localForces    = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfAtoms, 3, NULL, true, (RealOpenMM) 0.0 );

      // required if Gromacs float != RealOpenMM (i.e., RealOpenMM = double)

      partialCharges = new RealOpenMM[numberOfAtoms];
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         partialCharges[ii] = (RealOpenMM) partialChargesIn[ii];         
      }
   }

   // convert coordinates from nm -> Angstrom

   RealOpenMM distanceScaleFactor = (RealOpenMM) 10.0;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      coordinates[ii][0] = (RealOpenMM) (distanceScaleFactor*atomCoordinates[ii][0]);         
      coordinates[ii][1] = (RealOpenMM) (distanceScaleFactor*atomCoordinates[ii][1]);         
      coordinates[ii][2] = (RealOpenMM) (distanceScaleFactor*atomCoordinates[ii][2]);         
   }

   CpuImplicitSolvent::computeImplicitSolventForces( coordinates, partialCharges,
                                                     localForces, updateBornRadii );

   float forceScaleFactor = 0.4184f;
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      forces[ii][0] += (float) forceScaleFactor*localForces[ii][0];         
      forces[ii][1] += (float) forceScaleFactor*localForces[ii][1];         
      forces[ii][2] += (float) forceScaleFactor*localForces[ii][2];         
   }
//RealOpenMM energy = cpuImplicitSolvent->getEnergy(); 
//printf( "\ngromacsCpuCalculateObcForces E=%.5e", energy );

   return 0;

}

/**---------------------------------------------------------------------------------------

   Write Tinker Xyz file

   @param atomCoordinates   Gromacs atom coordinates ('x' in md.c)
   @param top               Gromacs t_topology (as in md.c)
   @param inputHeader       xyz file header
   @param inputXyzFileName  xyz file name
   @param log               log file reference (may be NULL)

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C" int cpuWriteXyzFile( const rvec *atomCoordinates, const t_topology* top,
                                const char* inputHeader, const char* inputXyzFileName, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName   = "\ncpuWriteXyzFile: ";

   // ---------------------------------------------------------------------------------------

   // set logging, if log != NULL

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   // set header & file name

   std::string header;
   if( inputHeader ){
      header      = std::string( inputHeader );
   } else {
      header      = "Test";
   }

   std::string xyzFileName;
   if( inputXyzFileName ){
      xyzFileName = std::string( inputXyzFileName );
   } else {
      xyzFileName = "Test";
   }

   // write file

   SimTKOpenMMGromacsUtilities::writeTinkerXyzFile( top->atoms.nr, atomCoordinates, 
                                                    header, xyzFileName, top );
   return 0;

}
