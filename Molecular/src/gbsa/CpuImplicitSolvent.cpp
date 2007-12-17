
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

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "CpuImplicitSolvent.h"

//#define UseGromacsMalloc 1

// Replacement new/delete w/ Gromac's smalloc() and sfree()

// static data member created by call to cpuSetImplicitSolventParameters() 
// stores parameter settings, ...
// used to calculate GBSA forces/energy

CpuImplicitSolvent* CpuImplicitSolvent::_cpuImplicitSolvent                     = NULL;

// info file related-stuff

std::string CpuImplicitSolvent::_defaultInfoFileName                            = std::string( "CpuImplicitSolventInfo" );

// key for info file: base file name 

const std::string CpuImplicitSolvent::CpuImplicitSolventBaseFileName            = std::string( "CpuImplicitSolventBaseFileName" );

// key for info file: file generation frequency 

const std::string CpuImplicitSolvent::CpuImplicitSolventFileGenerationFrequency = std::string( "CpuImplicitSolventFileGenerationFrequency" );


/**---------------------------------------------------------------------------------------

   CpuImplicitSolvent constructor

   @param implicitSolventParameters      implicitSolventParameters object
   
   --------------------------------------------------------------------------------------- */

CpuImplicitSolvent::CpuImplicitSolvent( ImplicitSolventParameters* implicitSolventParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::CpuImplicitSolvent(2)";

   // ---------------------------------------------------------------------------------------

   _initializeDataMembers( );
   _implicitSolventParameters = implicitSolventParameters;

}

/**---------------------------------------------------------------------------------------

   CpuImplicitSolvent destructor

   --------------------------------------------------------------------------------------- */

CpuImplicitSolvent::~CpuImplicitSolvent( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::~CpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

//#ifdef UseGromacsMalloc
//      save_free( "_bornForce", __FILE__, __LINE__, _bornForce );
//#else
//      delete[] _bornForce;
//#endif

   delete _implicitSolventParameters;

   delete[] _bornForce;
   delete[] _bornRadii;
   delete[] _tempBornRadii;

}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuImplicitSolvent::_initializeDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::_initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _bornRadii                           = NULL;
   _tempBornRadii                       = NULL;
   _bornForce                           = NULL;
 
   _includeAceApproximation             = 0;

   _forceConversionFactor               = (RealOpenMM) 1.0;
   _energyConversionFactor              = (RealOpenMM) 1.0;

   _forceCallIndex                      = 0;

   _implicitSolventEnergy               = (RealOpenMM) 0.0;

   _baseFileName                        = SimTKOpenMMCommon::NotSet;
   _outputFileFrequency                 = 1;
}

/**---------------------------------------------------------------------------------------

   Return number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::getNumberOfAtoms( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getNumberOfAtoms";

   // ---------------------------------------------------------------------------------------

   return _implicitSolventParameters->getNumberOfAtoms();

}

/**---------------------------------------------------------------------------------------

   Return energy 

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuImplicitSolvent::getEnergy( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getEnergy";

   // ---------------------------------------------------------------------------------------

   return _implicitSolventEnergy;

}

/**---------------------------------------------------------------------------------------

   Delete static _cpuImplicitSolvent object if set

   @return SimTKOpenMMCommon::DefaultReturn if _cpuImplicitSolvent was set; 
           otherwise return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::deleteCpuImplicitSolvent( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::deleteCpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   if( _cpuImplicitSolvent != NULL ){
      delete _cpuImplicitSolvent;
      _cpuImplicitSolvent = NULL;
      return SimTKOpenMMCommon::DefaultReturn;
   } else {
      return SimTKOpenMMCommon::ErrorReturn;
   }
}

/**---------------------------------------------------------------------------------------

   Set static member _cpuImplicitSolvent

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setCpuImplicitSolvent( CpuImplicitSolvent* cpuImplicitSolvent ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setCpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   _cpuImplicitSolvent = cpuImplicitSolvent;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get static member cpuImplicitSolvent

   @return static member cpuImplicitSolvent

   --------------------------------------------------------------------------------------- */

CpuImplicitSolvent* CpuImplicitSolvent::getCpuImplicitSolvent( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getCpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   return _cpuImplicitSolvent;

}
/**---------------------------------------------------------------------------------------

   Set energy 

   @param energy

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setEnergy( RealOpenMM energy ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setEnergy";

   // ---------------------------------------------------------------------------------------

   _implicitSolventEnergy = energy;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Return ImplicitSolventParameters

   @return ImplicitSolventParameters

   --------------------------------------------------------------------------------------- */

ImplicitSolventParameters* CpuImplicitSolvent::getImplicitSolventParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getImplicitSolventParameters";

   // ---------------------------------------------------------------------------------------

   return _implicitSolventParameters;

}

/**---------------------------------------------------------------------------------------

   Set ImplicitSolventParameters

   @param ImplicitSolventParameters

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setImplicitSolventParameters( ImplicitSolventParameters* implicitSolventParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setImplicitSolventParameters";

   // ---------------------------------------------------------------------------------------

   _implicitSolventParameters = implicitSolventParameters;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return flag signalling whether AceApproximation for nonpolar term is to be included

   @return flag

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::includeAceApproximation( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::includeAceApproximation";

   // ---------------------------------------------------------------------------------------

   return _includeAceApproximation;

}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether AceApproximation is to be included

   @param includeAceApproximation new includeAceApproximation value

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setIncludeAceApproximation( int includeAceApproximation ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setImplicitSolventParameters";

   // ---------------------------------------------------------------------------------------

   _includeAceApproximation = includeAceApproximation;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return ForceConversionFactor for units

   @return ForceConversionFactor

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuImplicitSolvent::getForceConversionFactor( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getForceConversionFactor";

   // ---------------------------------------------------------------------------------------

   return _forceConversionFactor;

}

/**---------------------------------------------------------------------------------------

   Set ForceConversionFactor

   @param ForceConversionFactor (units conversion)

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setForceConversionFactor( RealOpenMM forceConversionFactor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setForceConversionFactor";

   // ---------------------------------------------------------------------------------------

   _forceConversionFactor = forceConversionFactor;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return EnergyConversionFactor for units

   @return EnergyConversionFactor

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuImplicitSolvent::getEnergyConversionFactor( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getEnergyConversionFactor";

   // ---------------------------------------------------------------------------------------

   return _energyConversionFactor;

}

/**---------------------------------------------------------------------------------------

   Set EnergyConversionFactor

   @param EnergyConversionFactor (units conversion)

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setEnergyConversionFactor( RealOpenMM energyConversionFactor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setEnergyConversionFactor";

   // ---------------------------------------------------------------------------------------

   _energyConversionFactor = energyConversionFactor;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return ForceCallIndex -- number of times forces have been calculated

   @return ForceCallIndex

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::getForceCallIndex( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getForceCallIndex";

   // ---------------------------------------------------------------------------------------

   return _forceCallIndex;

}

/**---------------------------------------------------------------------------------------

   Increment ForceCallIndex

   @return incremented forceCallIndex

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::incrementForceCallIndex( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::incrementForceCallIndex";

   // ---------------------------------------------------------------------------------------

   _forceCallIndex++;

   return _forceCallIndex;

}

/**---------------------------------------------------------------------------------------

   Return bornForce, a work array of size _implicitSolventParameters->getNumberOfAtoms()*sizeof( RealOpenMM )
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

RealOpenMM* CpuImplicitSolvent::getBornForce( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getImplicitSolventBornForce";

   // ---------------------------------------------------------------------------------------

   if( _bornForce == NULL ){
      _bornForce = new RealOpenMM[_implicitSolventParameters->getNumberOfAtoms()];
   }
   return _bornForce;

}

/**---------------------------------------------------------------------------------------

   Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if it is not set

   @return array

   --------------------------------------------------------------------------------------- */

RealOpenMM* CpuImplicitSolvent::getBornRadii( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getBornRadii";

   // ---------------------------------------------------------------------------------------

   if( _bornRadii == NULL ){
      _bornRadii = new RealOpenMM[_implicitSolventParameters->getNumberOfAtoms()];
   }
   return _bornRadii;
}

/**---------------------------------------------------------------------------------------

   Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()

   @return array

   --------------------------------------------------------------------------------------- */

const RealOpenMM* CpuImplicitSolvent::getBornRadiiConst( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getBornRadii";

   // ---------------------------------------------------------------------------------------

   return _bornRadii;
}

/**---------------------------------------------------------------------------------------

   Return Born radii temp work array of size=_implicitSolventParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

RealOpenMM* CpuImplicitSolvent::getBornRadiiTemp( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getImplicitSolventBornRadiiTemp";

   // ---------------------------------------------------------------------------------------

   if( _tempBornRadii == NULL ){
      _tempBornRadii = new RealOpenMM[_implicitSolventParameters->getNumberOfAtoms()];
   }
   return _tempBornRadii;
}

/**---------------------------------------------------------------------------------------

   Compute Born radii

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii
   @param obcChain            output array of Obc chain derivatives

   @return SimTKOpenMMCommon::DefaultReturn or SimTKOpenMMCommon::ErrorReturn 
           if problems encountered

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii,
                                          RealOpenMM* obcChain ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuImplicitSolvent::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << methodName;
   message << " Error: calling from base class.";
   SimTKOpenMMLog::printError( message );
   return SimTKOpenMMCommon::ErrorReturn;

}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces

   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces (output)

   @return SimTKOpenMMCommon::DefaultReturn; abort if cpuImplicitSolvent is not set

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeImplicitSolventForces( RealOpenMM** atomCoordinates,
                                                      const RealOpenMM* partialCharges,
                                                      RealOpenMM** forces ){

   // ---------------------------------------------------------------------------------------

   int printSampleOutput         = 0;
   static const char* methodName = "\nCpuImplicitSolvent::computeImplicitSolventForces";

   // ---------------------------------------------------------------------------------------

   // check CpuImplicitSolvent initialized

   CpuImplicitSolvent* cpuImplicitSolvent = getCpuImplicitSolvent();
   if( cpuImplicitSolvent == NULL ){
      std::stringstream message;
      message << methodName;
      message << " CpuImplicitSolvent has not been initialized!";
      SimTKOpenMMLog::printError( message );
      return SimTKOpenMMCommon::ErrorReturn; 
   }

   int callId = cpuImplicitSolvent->incrementForceCallIndex();

   // check parameters have been initialized

   ImplicitSolventParameters* implicitSolventParameters =  cpuImplicitSolvent->getImplicitSolventParameters();
   if( !implicitSolventParameters ){
      std::stringstream message;
      message << methodName;
      message << " implicitSolventParameters has not been initialized!";
      SimTKOpenMMLog::printError( message );
      return SimTKOpenMMCommon::ErrorReturn; 
   }

   // check if parameters good-to-go
   // radii and scalefactors (OBC) need to be set

   if( callId == 1 ){
      if( implicitSolventParameters->isNotReady() ){
         std::stringstream message;
         message << methodName;
         message << " implicitSolventParameters are not set for force calculations!";
         SimTKOpenMMLog::printError( message );
         return SimTKOpenMMCommon::ErrorReturn; 
      } else {
         std::stringstream message;
         message << methodName;
         message << " implicitSolventParameters appear to be set.";
         // SimTKOpenMMLog::printMessage( message );
      }
   }

   // check to see if Born radii have been previously calculated
   // if not, then calculate;
   // logic here assumes that the radii are intitialized to zero 
   // and then once computed, always greater than zero.

   // after first iteration Born radii are updated in force calculation (computeBornEnergyForces())

   RealOpenMM* bornRadii = cpuImplicitSolvent->getBornRadii();
   if( bornRadii[0] < (RealOpenMM) 0.0001 || callId == 1 ){

      cpuImplicitSolvent->computeBornRadii( atomCoordinates, bornRadii );

      // diagnostics

      if( printSampleOutput ){

         RealOpenMM* atomicRadii = implicitSolventParameters->getAtomicRadii();
         std::stringstream message;
         message.precision( 6 );
         message.width( 12 );
         message << methodName;
         int numberOfAtoms = implicitSolventParameters->getNumberOfAtoms();
         message << " initialize Born radii for " << numberOfAtoms << " atoms on call=" << callId; 

         for( int ii = 0; ii < printSampleOutput && ii < numberOfAtoms; ii++ ){
		      message << "\n   " << ii << " rad=" << atomicRadii[ii] << " q=" << partialCharges[ii] << " bR=" << bornRadii[ii] << " X[";
            SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii] );
		      message << "]";
         }

		   message << "\n";

         int startIndex = implicitSolventParameters->getNumberOfAtoms() - printSampleOutput > 0 ?
                          implicitSolventParameters->getNumberOfAtoms() - printSampleOutput : numberOfAtoms;
         for( int ii = startIndex; ii < numberOfAtoms; ii++ ){
		      message << "\n   " << ii << " " << atomicRadii[ii] << " " << bornRadii[ii] << " X[";
            SimTKOpenMMUtilities::formatRealStringStream( message, atomCoordinates[ii] );
		      message << "]";
         }
         SimTKOpenMMLog::printMessage( message );
      }
   }

   // compute forces

   cpuImplicitSolvent->computeBornEnergyForces( cpuImplicitSolvent->getBornRadii(), atomCoordinates,
                                                partialCharges, forces );

   // diagnostics

   if( printSampleOutput && callId == 1 ){

      std::stringstream message;
      message.precision( 6 );
      message.width( 12 );
      int numberOfAtoms = implicitSolventParameters->getNumberOfAtoms();

      message << methodName;
      message << " call=" << callId << " E=" << cpuImplicitSolvent->getEnergy() << " " << numberOfAtoms << " atoms.";

      for( int ii = 0; ii < printSampleOutput && ii < numberOfAtoms; ii++ ){
         message << "\n   " << ii << " [ ";
         SimTKOpenMMUtilities::formatRealStringStream( message, forces[ii] );
         message << "] bRad=" << bornRadii[ii]; 
      }
		message << "\n";
      int startIndex = implicitSolventParameters->getNumberOfAtoms() - printSampleOutput > 0 ?
                       implicitSolventParameters->getNumberOfAtoms() - printSampleOutput : numberOfAtoms;
      for( int ii = startIndex; ii < numberOfAtoms; ii++ ){
         message << "\n   " << ii << " [ ";
         SimTKOpenMMUtilities::formatRealStringStream( message, forces[ii] );
         message << "] bRad=" << bornRadii[ii]; 
      }
      SimTKOpenMMLog::printMessage( message );

      // write Born forces

      std::stringstream resultsFileName;
      resultsFileName << cpuImplicitSolvent->getBaseFileName() << "." << callId << ".gbsa0";
      cpuImplicitSolvent->writeBornEnergyForces( atomCoordinates, partialCharges, forces, resultsFileName.str() );

   }

   return SimTKOpenMMCommon::DefaultReturn; 
}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   @param bornRadii           Born radii -- optional; if NULL, then ImplicitSolventParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return SimTKOpenMMCommon::ErrorReturn since the call should be implemented 
           in a derived class

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeBornEnergyForces( RealOpenMM* bornRadii,
                                                 RealOpenMM** atomCoordinates,
                                                 const RealOpenMM* partialCharges,
                                                 RealOpenMM** forces ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuImplicitSolvent::computeBornEnergyForces";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << methodName;
   message << " Error: calling from base class.";
   SimTKOpenMMLog::printError( message );
   return SimTKOpenMMCommon::ErrorReturn; 

}

/**---------------------------------------------------------------------------------------

   Get nonpolar solvation force constribution via ACE approximation

   @param implicitSolventParameters parameters
   @param vdwRadii                  Vdw radii
   @param bornRadii                 Born radii
   @param energy                    energy (output): value is incremented from input value 
   @param forces                    forces: values are incremented from input values

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeAceNonPolarForce( const ImplicitSolventParameters* implicitSolventParameters,
                                                 const RealOpenMM* bornRadii, RealOpenMM* energy,
                                                 RealOpenMM* forces ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::computeAceNonPolarForce";

   static const RealOpenMM minusSix = -6.0;

   // ---------------------------------------------------------------------------------------

   // compute the nonpolar solvation via ACE approximation

   const RealOpenMM probeRadius          = implicitSolventParameters->getProbeRadius();
   const RealOpenMM surfaceAreaFactor    = implicitSolventParameters->getPi4Asolv();
   const RealOpenMM* atomicRadii         = implicitSolventParameters->getAtomicRadii();
   int numberOfAtoms                     = implicitSolventParameters->getNumberOfAtoms();

   // 1 + 1 + pow + 3 + 1 + 2 FLOP

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      if( bornRadii[atomI] > 0.0 ){
         RealOpenMM r            = atomicRadii[atomI] + probeRadius;
         RealOpenMM ratio6       = POW( atomicRadii[atomI]/bornRadii[atomI], (RealOpenMM) 6.0 );
         RealOpenMM saTerm       = surfaceAreaFactor*r*r*ratio6;
         *energy                += saTerm;
         forces[atomI]          += minusSix*saTerm/bornRadii[atomI]; 
      }
   }

   return SimTKOpenMMCommon::DefaultReturn; 

}

/**---------------------------------------------------------------------------------------

   Write Born energy and forces (Simbios)

   @param atomCoordinates     atomic coordinates
   @param forces              forces
   @param resultsFileName     output file name

   @return SimTKOpenMMCommon::DefaultReturn unless
           file cannot be opened
           in which case return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::writeBornEnergyForces( RealOpenMM** atomCoordinates,
                                               const RealOpenMM* partialCharges,
                                               RealOpenMM** forces,
                                               const std::string& resultsFileName ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName  = "\nCpuImplicitSolvent::writeBornEnergyForces";

   // ---------------------------------------------------------------------------------------

   ImplicitSolventParameters* implicitSolventParameters = getImplicitSolventParameters();

   int numberOfAtoms              = implicitSolventParameters->getNumberOfAtoms();
   const RealOpenMM* atomicRadii  = implicitSolventParameters->getAtomicRadii();
   const RealOpenMM* bornRadii    = getBornRadiiConst();

   // ---------------------------------------------------------------------------------------

   // open file -- return if unsuccessful

   FILE* implicitSolventResultsFile = NULL;
#ifdef WIN32
   fopen_s( &implicitSolventResultsFile, resultsFileName.c_str(), "w" );
#else
   implicitSolventResultsFile = fopen( resultsFileName.c_str(), "w" );
#endif

   // diganostics

   std::stringstream message;
   message << methodName;
   if( implicitSolventResultsFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << resultsFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << "  could not open file=<" << resultsFileName << "> -- abort output.";
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;
   }

   // header

   (void) fprintf( implicitSolventResultsFile, "# %d atoms format: coords bornRadii q atomicRadii forces\n", numberOfAtoms );

   RealOpenMM forceConversion  = 1.0;
   RealOpenMM lengthConversion = 1.0;

   // output

   if( forces != NULL && atomCoordinates != NULL && partialCharges != NULL && atomicRadii != NULL ){
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
            (void) fprintf( implicitSolventResultsFile, "%.7e %.7e %.7e %.7e %.5f %.5f %.7e %.7e %.7e\n",
                            lengthConversion*atomCoordinates[ii][0],
                            lengthConversion*atomCoordinates[ii][1], 
                            lengthConversion*atomCoordinates[ii][2],
                           (bornRadii != NULL ? lengthConversion*bornRadii[ii] : 0.0),
                            partialCharges[ii], lengthConversion*atomicRadii[ii],
                            forceConversion*forces[ii][0],
                            forceConversion*forces[ii][1],
                            forceConversion*forces[ii][2]
                          );
      }
   }
   (void) fclose( implicitSolventResultsFile );

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get BaseFileName

   @return    baseFileName

   --------------------------------------------------------------------------------------- */

const std::string& CpuImplicitSolvent::getBaseFileName( void ) const {

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::getBaseFileName";

// ---------------------------------------------------------------------------------------

   return _baseFileName;

}

/**---------------------------------------------------------------------------------------

   Set BaseFileName

   @param    input baseFileName

   @return   SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setBaseFileName( const std::string& baseFileName ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::setBaseFileName";

// ---------------------------------------------------------------------------------------

   _baseFileName = baseFileName;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Set OutputFileFrequency

   @param    input outputFileFrequency

   @return   SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setOutputFileFrequency( int outputFileFrequency ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::setOutputFileFrequency";

// ---------------------------------------------------------------------------------------

   _outputFileFrequency = outputFileFrequency;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get OutputFileFrequency

   @return   output file frequency

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::getOutputFileFrequency( void ) const {

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::getOutputFileFrequency";

// ---------------------------------------------------------------------------------------

   return _outputFileFrequency;

}

/**---------------------------------------------------------------------------------------

   Read info file

   @param infoFileName       file name to read

   @return SimTKOpenMMCommon::DefaultReturn if ok; else SimTKOpenMMCommon::ErrorReturn; if major 
           problem, program will exit

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::readInfoFile( const std::string infoFileName ){

   // ---------------------------------------------------------------------------------------

   const int bufferSize = 2048;
   char buffer[bufferSize];

   static const std::string methodName = "\nCpuImplicitSolvent::readInfoFile";

   // ---------------------------------------------------------------------------------------

   FILE* infoFile = NULL;
#ifdef WIN32
   fopen_s( &infoFile, infoFileName.c_str(), "r" );
#else
   infoFile = fopen( infoFileName.c_str(), "r" );
#endif

   if( infoFile == NULL ){

#ifdef WIN32
      fopen_s( &infoFile, _defaultInfoFileName.c_str(), "r" );
#else
      infoFile = fopen( _defaultInfoFileName.c_str(), "r" );
#endif

      if( infoFile == NULL ){
         std::stringstream message;
         message << methodName << " fopen failed on info file=<" << _defaultInfoFileName << "> and file=<" << infoFileName << ">";
         SimTKOpenMMLog::printMessage( message );
         return SimTKOpenMMCommon::ErrorReturn;
      } else {
         std::stringstream message;
         message << methodName << " opened info file=<" << _defaultInfoFileName << ">.";
         SimTKOpenMMLog::printMessage( message );
      }
   } else {
      std::stringstream message;
      message << methodName << " opened info file=<" << infoFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   }

   /* Sample file:

   */

   int errors = 0;
   while( fgets( buffer, bufferSize, infoFile ) ){

      std::stringstream message;
      std::string fileName   = std::string( );
      size_t bufferLen       = strlen( buffer );
      if( bufferLen > 0 ){
         buffer[bufferLen-1] = '\0'; 
      }
      
      message << "\n<" << buffer << ">";

//fprintf( amoebaLog->getLogFile(), "\nreadInfoFile <%s>", buffer );
//fflush( amoebaLog->getLogFile() );

      StringVector tokens;
      SimTKOpenMMUtilities::tokenizeString( buffer, tokens );

      bool done = false;
      for( StringVectorI ii = tokens.begin(); ii != tokens.end() && !done; ii++ ){

         std::string token        = *ii;
         size_t tokenLength       = token.length();
         bool recognized          = false;

         // skip comments and blank lines

         if( tokenLength < 2 || !token.compare( "#" ) ){
            done       = true;
            recognized = true;

         // base file name

         } else if( !token.compare( CpuImplicitSolvent::CpuImplicitSolventBaseFileName ) ){
            recognized = true;
            setBaseFileName( *(++ii) );

         // file frequency

         } else if( !token.compare( CpuImplicitSolvent::CpuImplicitSolventFileGenerationFrequency ) ){
            recognized        = true;
            std::string value = *(++ii);
            if( SimTKOpenMMUtilities::isValidInteger( value ) ){ 
               setOutputFileFrequency( atoi( value.c_str() ) );
            } else {
               message << "\nToken=<" << token << "> is not a valid integer for key=" << CpuImplicitSolvent::CpuImplicitSolventFileGenerationFrequency;
            }

         // keys used by objects other than 'CpuImplicitSolvent'

         } else {
/*
            recognized      = true;
            std::string key = *ii;
            ii++;
            if( ii != tokens.end() ){
               _inputArguments[key] = *ii;
            } else {
               _inputArguments[key] = CpuImplicitSolventCommon::Comment;
            }
*/
         }
         if( !recognized ){
            message << "\nToken=<" << token << "> not recognized.";
         }
      }
      SimTKOpenMMLog::printMessage( message );
   }

   (void) fclose( infoFile );

   // report if errors

   if( errors ){
      std::stringstream message;
      message << "\nErrors -- aborting.\n";
      SimTKOpenMMLog::printMessage( message );
      exit(-1);
   } else {
      std::stringstream message;
      message << "\nNo errors detected parsing info file " << infoFileName << "\n";
/*
      message << "\nNumber of hash arguments=" << _amoebaInputArguments.size() << " <key,value> pairs listed below:\n";
      for( StringStringMapCI ii = _amoebaInputArguments.begin(); ii != _amoebaInputArguments.end(); ii++ ){
         message << "\n   <" << (*ii).first <<  ">=<" << (*ii).second << ">";
      }
*/
      message << "\n\n";
      SimTKOpenMMLog::printMessage( message );
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuImplicitSolvent::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   if( title ){
      message << title;
   }

   message << "\nImplicitSolvent info:";
   message << "\nForce conversion=" << getForceConversionFactor() << " Energy conversion=" << getEnergyConversionFactor();
   message << "\nInclude ACE approximation=";
   if( includeAceApproximation() ){
      message << "Y";
   } else {
      message << "N";
   }

   // get parameters

   message << getImplicitSolventParameters()->getStateString( NULL );

   return message.str();
}

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */

/*
void* CpuImplicitSolvent::operator new( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

//   (void) fprintf( stdout, "\nCpuImplicitSolvent new called -- size=%u", size );
//   (void) fflush( stdout );

   return ptr;
} */

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */

/*
void CpuImplicitSolvent::operator delete( void *ptr ){

   // (void) fprintf( stdout, "\nCpuImplicitSolvent delete called." );
   // (void) fflush( stdout );

   sfree( ptr ); 
} */

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */
/*
void* CpuImplicitSolvent::operator new[]( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

   // (void) fprintf( stdout, "\nCpuImplicitSolvent new[] called -- size=%u", size );
   // (void) fflush( stdout );

   return ptr;
} */

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */
/*
void CpuImplicitSolvent::operator delete[]( void *ptr ){

   // (void) fprintf( stdout, "\nCpuImplicitSolvent delete[] called." );
   // (void) fflush( stdout );

   sfree( ptr ); 
}
*/
