
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

#ifndef __CpuImplicitSolvent_H__
#define __CpuImplicitSolvent_H__

#include "ImplicitSolventParameters.h"

// ---------------------------------------------------------------------------------------

class CpuImplicitSolvent {

   public:

      // fields in info file

      const static std::string CpuImplicitSolventBaseFileName;
      const static std::string CpuImplicitSolventFileGenerationFrequency;

   private:

      // default info file name

      static std::string _defaultInfoFileName;

      // base file name

      std::string _baseFileName;

      // frequency to output diagnostic files

      int _outputFileFrequency;

      // used for direct calls 

      static CpuImplicitSolvent* _cpuImplicitSolvent;

      // parameters

      ImplicitSolventParameters* _implicitSolventParameters;

      // flag to signal whether ACE approximation
      // is to be included

      int _includeAceApproximation;

      // force index call 

      int _forceCallIndex;

      // work arrays

      RealOpenMM* _bornForce;

      // Born radii and force

      RealOpenMM* _bornRadii;
      RealOpenMM* _tempBornRadii;

      // convert units for energy/force

      RealOpenMM _forceConversionFactor;
      RealOpenMM _energyConversionFactor;

      // Ed, 2007-04-27: Store the energy internally

      RealOpenMM _implicitSolventEnergy; 

      /**---------------------------------------------------------------------------------------
      
         Initialize data members -- potentially more than
         one constructor, so centralize intialization here
      
         --------------------------------------------------------------------------------------- */

      void _initializeDataMembers( void );

   protected:

      /**---------------------------------------------------------------------------------------
      
         Return implicitSolventBornForce, a work array of size _implicitSolventParameters->getNumberOfAtoms()*sizeof( RealOpenMM )
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getBornForce( void );

      /**---------------------------------------------------------------------------------------
      
         Return Born radii temp work array of size=_implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getBornRadiiTemp( void );

      /**---------------------------------------------------------------------------------------
      
         Set energy 

         @param energy new energy
      
         ireturn SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

		int setEnergy( RealOpenMM energy );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuImplicitSolvent( ImplicitSolventParameters* implicitSolventParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuImplicitSolvent( );

      // override of new/delete -- used when run in PS3 framework(?)

      // static void* operator new( size_t size ); 
      // static void  operator delete( void *p );

      // static void* operator new[]( size_t size ); 
      // static void  operator delete[]( void *p );

      /**---------------------------------------------------------------------------------------
      
         Delete static _cpuImplicitSolvent object if set
      
         @return SimTKOpenMMCommon::DefaultReturn if _cpuImplicitSolvent was set; 
                 otherwise return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */

      static int deleteCpuImplicitSolvent( void );

      /**---------------------------------------------------------------------------------------
      
         Set static member _cpuImplicitSolvent
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      static int setCpuImplicitSolvent( CpuImplicitSolvent* cpuImplicitSolvent );

      /**---------------------------------------------------------------------------------------
      
         Get static member cpuImplicitSolvent
      
         @return static member cpuImplicitSolvent
      
         --------------------------------------------------------------------------------------- */
      
      static CpuImplicitSolvent* getCpuImplicitSolvent( void );

      /**---------------------------------------------------------------------------------------
      
         Get number of atoms
      
         @return number of atoms
      
         --------------------------------------------------------------------------------------- */

		int getNumberOfAtoms( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Get energy 
      
         @return energy
      
         --------------------------------------------------------------------------------------- */

		RealOpenMM getEnergy( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Return ImplicitSolventParameters
      
         @return ImplicitSolventParameters
      
         --------------------------------------------------------------------------------------- */
      
      ImplicitSolventParameters* getImplicitSolventParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setImplicitSolventParameters( ImplicitSolventParameters* implicitSolventParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Return flag signalling whether AceApproximation for nonpolar term is to be included
      
         @return flag
      
         --------------------------------------------------------------------------------------- */

      int includeAceApproximation( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set flag indicating whether AceApproximation is to be included
      
         @param includeAceApproximation new includeAceApproximation value
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setIncludeAceApproximation( int includeAceApproximation );

      /**---------------------------------------------------------------------------------------
      
         Return ForceConversionFactor for units
      
         @return ForceConversionFactor
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getForceConversionFactor(  void  ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ForceConversionFactor
      
         @param ForceConversionFactor (units conversion)
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setForceConversionFactor(  RealOpenMM forceConversionFactor  );

      /**---------------------------------------------------------------------------------------
      
         Return EnergyConversionFactor for units
      
         @return EnergyConversionFactor
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM getEnergyConversionFactor( void  ) const;

      /**---------------------------------------------------------------------------------------
      
         Set EnergyConversionFactor
      
         @param EnergyConversionFactor (units conversion)
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      int setEnergyConversionFactor( RealOpenMM energyConversionFactor );

      /**---------------------------------------------------------------------------------------
      
         Return ForceCallIndex -- number of times forces have been calculated
      
         @return ForceCallIndex
      
         --------------------------------------------------------------------------------------- */
      
      int getForceCallIndex( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Increment ForceCallIndex
      
         @return incremented forceCallIndex
      
         --------------------------------------------------------------------------------------- */
      
      int incrementForceCallIndex( void );
      
      /**---------------------------------------------------------------------------------------
      
         Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getBornRadii( void );
      
      /**---------------------------------------------------------------------------------------
      
         Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      const RealOpenMM* getBornRadiiConst( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces (output); if not set on input, then memory is allocated
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int computeImplicitSolventForces( RealOpenMM** atomCoordinates, const RealOpenMM* partialCharges,
                                               RealOpenMM** forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param atomCoordinates   atomic coordinates dimension: [0-numberAtoms-1][0-2]
         @param bornRadii         output array of Born radii
         @param obcChain          output array of OBC chain derivative
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      virtual int computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii, RealOpenMM* obcChain = NULL );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      virtual int computeBornEnergyForces( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                           const RealOpenMM* partialCharges, RealOpenMM** forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get nonpolar solvation force constribution via ACE approximation
      
         @param implicitSolventParameters parameters
         @param bornRadii                 Born radii
         @param energy                    energy (output): value is incremented from input value 
         @param forces                    forces: values are incremented from input values
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      int computeAceNonPolarForce( const ImplicitSolventParameters* implicitSolventParameters,
                                   const RealOpenMM* bornRadii, RealOpenMM* energy, 
                                   RealOpenMM* forces ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Write Born energy and forces (Simbios)
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial atom charges
         @param forces            force array
         @param resultsFileName   output file name
      
         @return SimTKOpenMMCommon::DefaultReturn if file opened; else return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      virtual int writeBornEnergyForces( RealOpenMM** atomCoordinates,
                                         const RealOpenMM* partialCharges, RealOpenMM** forces,
                                         const std::string& resultsFileName ) const;

      /**---------------------------------------------------------------------------------------
            
         Get string w/ state 
         
         @param title               title (optional)
            
         @return string containing state
            
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;
      
      /**---------------------------------------------------------------------------------------

         Get BaseFileName

         @return    baseFileName
      
         --------------------------------------------------------------------------------------- */

      const std::string& getBaseFileName( void ) const;

      /**---------------------------------------------------------------------------------------

         Set BaseFileName

         @param    input baseFileName

         @return   SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setBaseFileName( const std::string& baseFileName );

      /**---------------------------------------------------------------------------------------
      
         Get OutputFileFrequency
      
         @return    outputFileFrequency
      
         --------------------------------------------------------------------------------------- */
     
      int getOutputFileFrequency( void ) const;
     
      /**---------------------------------------------------------------------------------------
      
         Set OutputFileFrequency
      
         @param    input outputFileFrequency
      
         @return   SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
     
      int setOutputFileFrequency( int outputFileFrequency );

      /**---------------------------------------------------------------------------------------
      
         Read info file
      
         @param infoFileName       file name to read
      
         @return CpuImplicitSolventCommon::DefaultReturn if ok; else CpuImplicitSolventCommon::ErrorReturn; if major 
                 problem, program will exit
      
         --------------------------------------------------------------------------------------- */
      
      int readInfoFile( const std::string infoFileName );
      
      /**---------------------------------------------------------------------------------------
      
         Get output file name (helper method) (Simbios)
      
         fileName = (inputFileName || getBaseFileName()) . getForceCallIndex() . suffix
      
         if inputFileName is NULL, then use getBaseFileName() as base file name
      
         if getForceCallIndex() > 0, insert '.getForceCallIndex()' 
      
         append suffix, if not NULL
      
         @param inputFileName                inputFileName
         @param suffix                       file suffix
      
         @return file name
      
         --------------------------------------------------------------------------------------- */

      std::string getOutputFileName( const std::string* inputFileName, const std::string* suffix ) const;

      /**---------------------------------------------------------------------------------------
      
         Write Tinker xyz file (Simbios)
      
         @param numberOfAtoms      number of atoms
         @param atomCoordinates    atom coordinates
         @param atomNames          atom names
         @param header             header
         @param xyzFileName        output file name
         @param bondsArray         bond array -- used to print 1-2 bonds
      
         @return 0 unless error detected
      
         Currently no attempt is made to get the atom name/type to accurately 
         reflect the Tinker names/types. Rather method is used to output atoms
         in Gromacs order and then reorder those in a corresponding xyz file
         w/ the correct atom names/types so that they match the Gromacs order
         This makes it easier to compare results between Gromacs and Tinker
      
         --------------------------------------------------------------------------------------- */
      
      /*
      static int writeXyzFile( int numberOfAtoms, const RealOpenMM** atomCoordinates, 
                               const char** atomNames,
                               const std::string& header, const std::string& xyzFileName,
                               const implicitSolventBonds** bondsArray ); */
      
};

// ---------------------------------------------------------------------------------------

#endif // __CpuImplicitSolvent_H__
