
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

#ifndef __CpuObc_H__
#define __CpuObc_H__

#include "ObcParameters.h"
#include "CpuImplicitSolvent.h"

// ---------------------------------------------------------------------------------------

class CpuObc : public CpuImplicitSolvent {

   private:

      // GBSA/OBC parameters

      ObcParameters* _obcParameters;

      // arrays containing OBC chain derivative 

      RealOpenMM* _obcChain;
      RealOpenMM* _obcChainTemp;

      // initialize data members (more than
      // one constructor, so centralize intialization here)

      void _initializeObcDataMembers( void );

   public:

      /**---------------------------------------------------------------------------------------
      
         Constructor
      
         @param implicitSolventParameters    ImplicitSolventParameters reference
      
         @return CpuImplicitSolvent object
      
         --------------------------------------------------------------------------------------- */

       CpuObc( ImplicitSolventParameters* obcParameters );

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~CpuObc( );

      /**---------------------------------------------------------------------------------------
      
         Return ObcParameters
      
         @return ObcParameters
      
         --------------------------------------------------------------------------------------- */

      ObcParameters* getObcParameters( void ) const;

      /**---------------------------------------------------------------------------------------
      
         Set ImplicitSolventParameters
      
         @param ImplicitSolventParameters
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */

      int setObcParameters( ObcParameters* obcParameters );
 
      /**---------------------------------------------------------------------------------------
      
         Return OBC chain derivative: size = _implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getObcChain( void );
      RealOpenMM* getObcChainConst( void ) const;
      
      /**---------------------------------------------------------------------------------------
      
         Return OBC chain temp work array of size=_implicitSolventParameters->getNumberOfAtoms()
         On first call, memory for array is allocated if not set
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      RealOpenMM* getObcChainTemp( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born radii based on OBC 

         @param atomCoordinates   atomic coordinates
         @param bornRadii         output array of Born radii
         @param obcChain          output array of OBC chain derivative factors; if NULL,
                                  then ignored
      
         @return array of Born radii
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornRadii( RealOpenMM** atomCoordinates, RealOpenMM* bornRadii,
                            RealOpenMM* obcChain = NULL );
      
      /**---------------------------------------------------------------------------------------
      
         Get Born energy and forces based on OBC 
      
         @param bornRadii         Born radii
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial charges
         @param forces            forces
      
         @return force array
      
         --------------------------------------------------------------------------------------- */
      
      int computeBornEnergyForces( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                   const RealOpenMM* partialCharges, RealOpenMM** forces );
      
      int computeBornEnergyForcesPrint( RealOpenMM* bornRadii, RealOpenMM** atomCoordinates,
                                        const RealOpenMM* partialCharges, RealOpenMM** forces );
      
      /**---------------------------------------------------------------------------------------
      
         Get state 
      
         title             title (optional)
      
         @return state string
      
         --------------------------------------------------------------------------------------- */
      
      std::string getStateString( const char* title ) const;

      /**---------------------------------------------------------------------------------------
      
         Write Born energy and forces (Simbios)
      
         @param atomCoordinates   atomic coordinates
         @param partialCharges    partial atom charges
         @param forces            force array
         @param resultsFileName   output file name
      
         @return SimTKOpenMMCommon::DefaultReturn if file opened; else return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
          
      int writeBornEnergyForces( RealOpenMM** atomCoordinates,
                                 const RealOpenMM* partialCharges, RealOpenMM** forces,
                                 const std::string& resultsFileName ) const;

      /**---------------------------------------------------------------------------------------
      
         Write  results from first loop
      
         @param atomCoordinates     atomic coordinates
         @param RealOpenMM forces         forces
         @param outputFileName      output file name
      
         @return SimTKOpenMMCommon::DefaultReturn unless
                 file cannot be opened
                 in which case return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int writeForceLoop1( int numberOfAtoms, RealOpenMM** forces, const RealOpenMM* bornForce,
                                  const std::string& outputFileName );
      
      /**---------------------------------------------------------------------------------------
      
         Write results
      
         @param numberOfAtoms       number of atoms
         @param chunkSizes          vector of chunk sizes for realRealOpenMMVector
         @param realRealOpenMMVector      vector of RealOpenMM**
         @param realVector          vector of RealOpenMM*
         @param outputFileName      output file name
      
         @return SimTKOpenMMCommon::DefaultReturn unless
                 file cannot be opened
                 in which case return SimTKOpenMMCommon::ErrorReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int writeForceLoop( int numberOfAtoms, const IntVector& chunkSizes,
                                 const RealOpenMMPtrPtrVector& realRealOpenMMVector, 
                                 const RealOpenMMPtrVector& realVector,
                                 const std::string& outputFileName );
      
};

// ---------------------------------------------------------------------------------------

#endif // __CpuObc_H__
