
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
#include "ImplicitSolventParameters.h"
#include <sstream>

/**---------------------------------------------------------------------------------------

   ImplicitSolventParameters:

		Calculates for each atom

			(1) the van der Waal radii
         (2) volume
         (3) fixed terms in ImplicitSolvent equation gPol
         (4) list of atoms that should be excluded in calculating
				 force -- nonbonded atoms (1-2, and 1-3 atoms)

	Implementation:

		Slightly different sequence of calls when running on CPU vs GPU.
		Difference arise because the CPU-side data arrays for the Brook
		streams are allocated by the BrookStreamWrapper objects. These
		arrays are then used by ImplicitSolventParameters when initializing the
		the values (vdwRadii, volume, ...) to be used in the calculation.

		Cpu:
			 ImplicitSolventParameters* implicitSolventParameters = new ImplicitSolventParameters( numberOfAtoms, log );
          implicitSolventParameters->initializeParameters( top );

		Gpu:

			implicitSolventParameters   = new ImplicitSolventParameters( gpu->natoms, log );
			
			// set arrays for cpu using stream data field; 
			// initializeParameters() only allocates space for arrays if they are not set (==NULL)
			// also set flag so that ImplicitSolventParameters destructor does not free arrays 
			
			implicitSolventParameters->setVdwRadii(  getBrookStreamWrapperAtIndex( GpuImplicitSolvent::implicitSolventVdwRadii  )->getData() );
			implicitSolventParameters->setVolume(    getBrookStreamWrapperAtIndex( GpuImplicitSolvent::implicitSolventVolume    )->getData() );
			implicitSolventParameters->setGPolFixed( getBrookStreamWrapperAtIndex( GpuImplicitSolvent::implicitSolventGpolFixed )->getData() );
			implicitSolventParameters->setAtomicRadii( getBrookStreamWrapperAtIndex( GpuImplicitSolvent::implicitSolventAtomicRadii )->getData() );
			
			implicitSolventParameters->setFreeArrays( false );
			
			implicitSolventParameters->initializeParameters( top );
 

   Issues:

		Tinker's atom radii are used. 
      The logic for mapping the Gromacs atom names to Tinker type may be incomplete;
      only tested for generic proteins
		see mapGmxAtomNameToTinkerAtomNumber()

   --------------------------------------------------------------------------------------- */


/**---------------------------------------------------------------------------------------

   ImplicitSolventParameters constructor (Simbios) 

   @param numberOfAtoms       number of atoms

   --------------------------------------------------------------------------------------- */

ImplicitSolventParameters::ImplicitSolventParameters( int numberOfAtoms ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::ImplicitSolventParameters";
   
   // ---------------------------------------------------------------------------------------

   _numberOfAtoms          = numberOfAtoms;

   _ownAtomicRadii         = false;
   _atomicRadii            = NULL;

   // see comments in ~ImplicitSolventParameters for explanation

   _freeArrays             = false;

   _initializeImplicitSolventConstants();

}

/**---------------------------------------------------------------------------------------

   ImplicitSolventParameters destructor (Simbios) 

   --------------------------------------------------------------------------------------- */

ImplicitSolventParameters::~ImplicitSolventParameters( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::~ImplicitSolventParameters";
   
   // ---------------------------------------------------------------------------------------

   if( _ownAtomicRadii ){
      delete[] _atomicRadii;
   }

}

/**---------------------------------------------------------------------------------------

   Get number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::getNumberOfAtoms( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getNumberOfAtoms:";

   // ---------------------------------------------------------------------------------------

   return _numberOfAtoms;
}

/**---------------------------------------------------------------------------------------

   Get solvent dielectric

   @return solvent dielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getSolventDielectric( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getSolventDielectric:";

   // ---------------------------------------------------------------------------------------

   return _solventDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solvent dielectric

   @param solventDielectric solvent dielectric

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setSolventDielectric( RealOpenMM solventDielectric ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setSolventDielectric:";

   // ---------------------------------------------------------------------------------------

   _solventDielectric = solventDielectric;
   _resetPreFactor();
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get solute dielectric

   @return soluteDielectric

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getSoluteDielectric( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getSoluteDielectric:";

   // ---------------------------------------------------------------------------------------

   return _soluteDielectric;
}

/**---------------------------------------------------------------------------------------

   Set solute dielectric

   @param soluteDielectric solute dielectric

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setSoluteDielectric( RealOpenMM soluteDielectric ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setSoluteDielectric:";

   // ---------------------------------------------------------------------------------------

   _soluteDielectric = soluteDielectric;
   _resetPreFactor();
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get electric constant (Simbios) 

   @return electricConstant

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getElectricConstant( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getElectricConstant:";

   // ---------------------------------------------------------------------------------------

   return _electricConstant;
}

/**---------------------------------------------------------------------------------------

   Set electric constant (Simbios) 

   @param electricConstant			electric constant

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setElectricConstant( RealOpenMM electricConstant ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setElectricConstant:";

   // ---------------------------------------------------------------------------------------

   _electricConstant = electricConstant;
   _resetPreFactor();

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get probe radius (Simbios) 

   @return probeRadius

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getProbeRadius( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getProbeRadius:";

   // ---------------------------------------------------------------------------------------

   return _probeRadius;
}

/**---------------------------------------------------------------------------------------

   Set probe radius (Simbios) 

   @param probeRadius	probe radius

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setProbeRadius( RealOpenMM probeRadius ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setProbeRadius:";

   // ---------------------------------------------------------------------------------------

   _probeRadius = probeRadius;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get pi*4*Asolv:  used in ACE approximation for nonpolar term  
         ((RealOpenMM) M_PI)*4.0f*0.0049*1000.0; (Still) 
         ((RealOpenMM) M_PI)*4.0f*0.0054*1000.0; (OBC) 

   @return pi4Asolv

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getPi4Asolv( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getPi4Asolv:";

   // ---------------------------------------------------------------------------------------

   return _pi4Asolv;
}

/**---------------------------------------------------------------------------------------

   Get prefactor

   @return prefactor

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getPreFactor( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getPreFactor:";

   // ---------------------------------------------------------------------------------------

   return _preFactor;
}

/**---------------------------------------------------------------------------------------

   Set pi4Asolv:  used in ACE approximation for nonpolar term  
         ((RealOpenMM) M_PI)*4.0f*0.0049*1000.0; (Still) 
         ((RealOpenMM) M_PI)*4.0f*0.0054*1000.0; (OBC) 

   @param pi4Asolv	probe radius

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setPi4Asolv( RealOpenMM pi4Asolv ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setPi4Asolv:";

   // ---------------------------------------------------------------------------------------

   _pi4Asolv = pi4Asolv;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get conversion factor for kcal/A -> kJ/nm (Simbios) 

   @return kcalA_To_kJNm factor

   --------------------------------------------------------------------------------------- */

RealOpenMM ImplicitSolventParameters::getKcalA_To_kJNm( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getKcalA_To_kJNm:";

   // ---------------------------------------------------------------------------------------

   return _kcalA_To_kJNm;
}

/**---------------------------------------------------------------------------------------

   Set KcalA_To_kJNm

   @param kcalA_To_kJNm	probe radius

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setKcalA_To_kJNm( RealOpenMM kcalA_To_kJNm ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setKcalA_To_kJNm:";

   // ---------------------------------------------------------------------------------------

   _kcalA_To_kJNm = kcalA_To_kJNm;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get AtomicRadii array

   @return array of atomic radii

   --------------------------------------------------------------------------------------- */

RealOpenMM* ImplicitSolventParameters::getAtomicRadii( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   if( _atomicRadii == NULL ){
      ImplicitSolventParameters* localThis = const_cast<ImplicitSolventParameters* const>(this);
      localThis->_atomicRadii   = new RealOpenMM[getNumberOfAtoms()];
      localThis->_ownAtomicRadii = true;
      memset( localThis->_atomicRadii, 0, sizeof( RealOpenMM )*getNumberOfAtoms() );
   }
   return _atomicRadii;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii array of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setAtomicRadii( RealOpenMM* atomicRadii ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   if( _ownAtomicRadii && atomicRadii != _atomicRadii ){
      delete[] _atomicRadii;
      _ownAtomicRadii = false;
   }
   _atomicRadii = atomicRadii;

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii array of atomic radii
   @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                  SimTKOpenMMCommon::MdUnits 

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setAtomicRadii( const RealOpenMMVector& atomicRadii, int units ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nImplicitSolventParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = getNumberOfAtoms();

   if( _ownAtomicRadii ){
      delete[] _atomicRadii;
      _atomicRadii = new RealOpenMM[numberOfAtoms];
   } else if( _atomicRadii == NULL ){
      _atomicRadii    = new RealOpenMM[numberOfAtoms];
      _ownAtomicRadii = true;
   }
   
   if( numberOfAtoms != (int) atomicRadii.size() ){
      std::stringstream message;
      message << methodName;
      message << " number of object atoms=" << numberOfAtoms << " does not match number in input vector=" << atomicRadii.size();
      SimTKOpenMMLog::printWarning( message );
      numberOfAtoms = numberOfAtoms < (int) atomicRadii.size() ? numberOfAtoms : (int) atomicRadii.size();
   }

   // force kcal/A units

   if( units == SimTKOpenMMCommon::MdUnits ){
      RealOpenMM ten = (RealOpenMM) 10.0;
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         _atomicRadii[ii] = ten*atomicRadii[ii];
      }
   } else {
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         _atomicRadii[ii] = atomicRadii[ii];
      }
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii array of atomic radii
   @param units       units flag: SimTKOpenMMCommon::KcalAngUnits or
                                  SimTKOpenMMCommon::MdUnits 

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setAtomicRadii( RealOpenMM* atomicRadii, int units ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nImplicitSolventParameters::setAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms = getNumberOfAtoms();

   if( _ownAtomicRadii ){
      delete[] _atomicRadii;
      _atomicRadii = new RealOpenMM[numberOfAtoms];
   } else if( _atomicRadii == NULL ){
      _atomicRadii    = new RealOpenMM[numberOfAtoms];
      _ownAtomicRadii = true;
   }
   
   // force kcal/A units

   if( units == SimTKOpenMMCommon::MdUnits ){
      RealOpenMM ten = (RealOpenMM) 10.0;
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         _atomicRadii[ii] = ten*atomicRadii[ii];
      }
   } else {
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         _atomicRadii[ii] = atomicRadii[ii];
      }
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set AtomicRadii array

   @param atomicRadii vector of atomic radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setOwnAtomicRadii( int ownAtomicRadii ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName = "\nImplicitSolventParameters::setOwnAtomicRadii:";

   // ---------------------------------------------------------------------------------------

   _ownAtomicRadii = ownAtomicRadii;

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get free array flag -- if set then work arrays are freed when destructor is called

   @return freeArrays flag

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::getFreeArrays( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setFreeArrays:";

   // ---------------------------------------------------------------------------------------

    return _freeArrays;

}

/**---------------------------------------------------------------------------------------

   Set free array flag -- if set then work arrays are freed when destructor is called

   @param freeArrays flag

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ImplicitSolventParameters::setFreeArrays( int freeArrays ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::setFreeArrays:";

   // ---------------------------------------------------------------------------------------

    _freeArrays = freeArrays;

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Initialize ImplicitSolvent Parameters (Simbios) 

   --------------------------------------------------------------------------------------- */

void ImplicitSolventParameters::_initializeImplicitSolventConstants( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::_initializeImplicitSolventConstants:";

   // ---------------------------------------------------------------------------------------

   _soluteDielectric        = (RealOpenMM)    1.0;
   _solventDielectric       = (RealOpenMM)   78.3;
   _kcalA_To_kJNm           = (RealOpenMM)    0.4184;
   _probeRadius             = (RealOpenMM)    1.4;
   _electricConstant        = (RealOpenMM) -166.02691;
   // _pi4Asolv                = (RealOpenMM) PI_M*4.0*0.0049*1000.0;
   //_pi4Asolv                = (RealOpenMM) PI_M*19.6;

   // _pi4Asolv                = (RealOpenMM) PI_M*4.0*0.0054;
   _pi4Asolv                = (RealOpenMM) PI_M*0.0216;
   
   _resetPreFactor();
}

/**---------------------------------------------------------------------------------------

   Reset prefactor (Simbios) 

	called when _electricConstant, _soluteDielectric, or _solventDielectric are modified

   --------------------------------------------------------------------------------------- */

void ImplicitSolventParameters::_resetPreFactor( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::_resetPreFactor:";

   static const RealOpenMM zero = 0.0;
   static const RealOpenMM one  = 1.0;
   static const RealOpenMM two  = 2.0;

   // ---------------------------------------------------------------------------------------

   if( getSoluteDielectric() != zero && getSolventDielectric() != zero ){
      _preFactor = two*getElectricConstant()*( (one/getSoluteDielectric()) - (one/getSolventDielectric()) );
   } else {
      _preFactor = zero;
   }
}

/**---------------------------------------------------------------------------------------
      
   Get state (Simbios)
   
   @param title               title (optional)
      
   @return SimTKOpenMMCommon::DefaultReturn always
      
   --------------------------------------------------------------------------------------- */

std::string ImplicitSolventParameters::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   if( title ){
      message << title;
   }
   std::string tab = getStringTab();

   message << "\nImplicitSolventParameters :";

   message << tab << "Solute dielectric:    " << getSoluteDielectric();
   message << tab << "Solvent dielectric:   " << getSolventDielectric();
   message << tab << "Electric constant:    " << getElectricConstant();
   message << tab << "Probe radius:         " << getProbeRadius();
   message << tab << "Number of atoms:      " << getNumberOfAtoms();
   message << tab << "Free arrays:          " << getFreeArrays();

   return message.str();

}

/**---------------------------------------------------------------------------------------
            
   Return zero value if all parameters set; else return nonzero
         
   @return ready status
            
   --------------------------------------------------------------------------------------- */
 
int ImplicitSolventParameters::isNotReady( void ) const {

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nImplicitSolventParameters::isNotReady";

   // ---------------------------------------------------------------------------------------

   int errors  = 0;
   int warning = 0;
   std::stringstream message;
   message << methodName;

   if( getNumberOfAtoms() <= 0 ){
      errors++;
      message << "\n   number of atoms=" << getNumberOfAtoms() << " is invalid.";
   }

   RealOpenMM* atomicRadii = getAtomicRadii();
   if( atomicRadii == NULL ){
      errors++;
      message << "\n   scaledRadiusFactors is not set";
   }

   // check radii are in correct units

   RealOpenMM average, stdDev, maxValue, minValue;
   int minIndex, maxIndex;
   SimTKOpenMMUtilities::getArrayStatistics( getNumberOfAtoms(), atomicRadii, &average,
                                             &stdDev, &minValue, &minIndex,
                                             &maxValue, &maxIndex );

   if( average < 1.0 || average > 10.0 || minValue < 0.5 ){
      errors++;
      message << "\n   atomic radii appear not to be set correctly -- radii should be in Angstroms";
      message << "\n   average radius=" << average << " min radius=" << minValue << " at atom index=" << minIndex;
   }


   if( getPreFactor() == 0.0 ){
      errors++;
      message << "\n   prefactor is not set.";
   }

   if( getSolventDielectric() <= 0.0 ){
      errors++;
      message << "\n   solvent dielectric=" << getSolventDielectric() << " is invalid.";
   }

   if( getSolventDielectric() >= 1000.0 ){
      warning++;
      message << "\n   Warning: solvent dielectric=" << getSolventDielectric() << " is large.";
   }

   if( getSoluteDielectric() <= 0.0 ){
      errors++;
      message << "\n   solute dielectric=" << getSoluteDielectric() << " is invalid.";
   }

   if( getSoluteDielectric() >= 1000.0 ){
      warning++;
      message << "\n   Warning: solute dielectric=" << getSoluteDielectric() << " is large.";
   }

   if( getElectricConstant() >= 0.0 ){
      errors++;
      message << "\n   electric constant=" << getElectricConstant() << " is invalid.";
   }

   if( getElectricConstant() >= 1000.0 ){
      warning++;
      message << "\n   Warning: electric constant=" << getElectricConstant() << " is large.";
   }

   if( getProbeRadius() <= 0.0 ){
      errors++;
      message << "\n   probe radius=" << getProbeRadius() << " is invalid.";
   }

   if( getProbeRadius() >= 10.0 ){
      warning++;
      message << "\n   Warning: probe radius=" << getProbeRadius() << " is large.";
   }

   if( getPi4Asolv() <= 0.0 ){
      errors++;
      message << "\n   Pi4Asolv=" << getPi4Asolv() << " is invalid.";
   }

   if( getPi4Asolv() >= 1000.0 ){
      warning++;
      message << "\n   Warning: Pi4Asolv=" << getPi4Asolv() << " is large.";
   }

   if( errors || warning ){
      message << std::endl;
      SimTKOpenMMLog::printMessage( message );
   }

   return errors;
}

/**---------------------------------------------------------------------------------------
      
   Get string tab -- centralized
   
   @return tab string
      
   --------------------------------------------------------------------------------------- */

std::string ImplicitSolventParameters::getStringTab( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nImplicitSolventParameters::getStringTab";

   // ---------------------------------------------------------------------------------------

   return std::string( "\n   " );

}
