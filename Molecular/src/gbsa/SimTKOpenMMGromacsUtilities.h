
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

#ifndef __SimTKOpenMMGromacsUtilities_H_
#define __SimTKOpenMMGromacsUtilities_H_

// ---------------------------------------------------------------------------------------

// class of shared, static utility methods

#include "SimTKOpenMMCommon.h" 

#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>

// ---------------------------------------------------------------------------------------

// Gromac's types

#include "typedefs.h" 

// ---------------------------------------------------------------------------------------

#define ATOM_ID_STRING_TAB 12
#define MAX_DEBUG_FIELDS   20

/**---------------------------------------------------------------------------------------

   Class of static methods to be shared
   Most methods are standalone 'utility' methods

--------------------------------------------------------------------------------------- */

class SimTKOpenMMGromacsUtilities {

   public:

      // dummy constructor/destructor

       SimTKOpenMMGromacsUtilities(){};
      ~SimTKOpenMMGromacsUtilities(){};

      /**---------------------------------------------------------------------------------------
      
         Find distances**2 from a given atom (Simbios)
      
         @param atomCoordinates     atom coordinates
         @param atomIndex           atom index to find distances from
         @param numberOfAtoms       number of atoms
         @param distances           array of distances squared on @return; array size must be at least
                                    numberOfAtoms
      
         @return distances
      
         --------------------------------------------------------------------------------------- */
      
      static int getDistanceSquaredFromSpecifiedAtom( const rvec* atomCoordinates, int atomIndex,
                                                      int numberOfAtoms, float* distances );
      
      /**---------------------------------------------------------------------------------------
      
         Find distances**2 from a given point (Simbios)
      
         @param atomCoordinates     atom coordinates
         @param point               point to find distances from
         @param numberOfAtoms       number of atoms
         @param distances           array of distances squared on @return; array size must be at least
                                    numberOfAtoms
      
         @return distances
      
         --------------------------------------------------------------------------------------- */
      
      static int getDistanceSquaredFromSpecifiedPoint( const rvec* atomCoordinates, float* point, 
                                                       int numberOfAtoms, float* distances );

      /**---------------------------------------------------------------------------------------
      
         Get atom name from top data struct
      
         @param atomIndex           atom index
         @param outputAtomName      output atom name
         @param top                 GMX t_topology struct
      
         @return SimTKOpenMMCommon::DefaultReturn

         --------------------------------------------------------------------------------------- */
      
      static int getAtomNameGivenAtomIndex( int atomIndex, char* outputAtomName, const t_topology* top );

      /**---------------------------------------------------------------------------------------
      
         Get residue name from top data struct given atom index
      
         @param atomIndex           atom index
         @param top                 GMX t_topology struct
         @param outputResidueName   output residue name
         @param outputResidueIndex  if not null, then *outputResidueIndex is residue index
      
         @return SimTKOpenMMCommon::DefaultReturn

         --------------------------------------------------------------------------------------- */
      
      static int getResidueNameGivenAtomIndex( int atomIndex, const t_topology* top,
                                               char* outputResidueName, int* outputResidueIndex );

      /**---------------------------------------------------------------------------------------
      
         Get atom name from top data struct
      
         @param atomIndex           atom index
         @param top                 GMX t_topology struct
         @param buffer              output buffer (enough space should have been reserved) 
         @param maxAtoms            max number of atoms for this run (may change -- used mainly
                                    to keep from reallocating cache array)
         @param tab                 tab spacing

         @return SimTKOpenMMCommon::DefaultReturn

         --------------------------------------------------------------------------------------- */
      
      static int getAtomIdStringGivenAtomIndex( int atomIndex, const t_topology* top,
                                                int sizeOfBuffer, char* buffer,
                                                int maxAtoms, unsigned int tab );

      /**---------------------------------------------------------------------------------------
      
         Get (1-2) bonds (Simbios) 
      
         @param maxAtoms            max number of atoms
         @param IntSetVector        vector of integer sets
         @param top                 Gromacs t_topolgy struct
      
         @return 0 if no errors or
         return x, where x is the number of errors encountered
      
         covalentBonds[i] = set of atom indices that are covalent partners
      
         --------------------------------------------------------------------------------------- */
      
      static int getCovalentBondIndices( int maxAtoms, IntSetVector& covalentBonds,
                                         const t_topology* top );
      
      /**---------------------------------------------------------------------------------------
      
         Get SETTLE stretch (1-2) bonds (Simbios) 
      
         @param maxAtoms            max number of atoms
         @param IntSetVector        vector of integer sets
         @param top                 Gromacs t_topolgy struct
      
         @return SimTKOpenMMCommon::DefaultReturn if no errors or
         return x, where x is the number of errors encountered
      
         --------------------------------------------------------------------------------------- */
      
      static int getSettleCovalentBondIndices( int maxAtoms, IntSetVector& covalentBonds,
                                               const t_topology* top );

      /**---------------------------------------------------------------------------------------
      
         Write Tinker xyz file (Simbios)
      
         @param numberOfAtoms        number of atoms
         @param atomCoordinates      atom coordinates
         @param header               header
         @param xyzFileName          output file name
         @param top                  Gromacs topology struct
      
         @return 0 unless error detected
      
         Currently no attempt is made to get the atom name/type to accurately 
         reflect the Tinker names/types. Rather method is used to output atoms
         in Gromacs order and then reorder those in a corresponding xyz file
         w/ the correct atom names/types so that they match the Gromacs order
         This makes it easier to compare results between Gromacs and Tinker
      
         --------------------------------------------------------------------------------------- */
     
      static int writeTinkerXyzFile( int numberOfAtoms, const rvec* atomCoordinates, 
                                     const std::string& header, const std::string& xyzFileName,
                                     const t_topology* top );
      
      /**---------------------------------------------------------------------------------------
      
         Get Tinker biotypes (Simbios) 
      
         @param numberOfAtoms        number of atoms
         @param top                  Gromacs topology struct
         @param tinkerAtomNames      Tinker atom names upon return
         @param tinkerResidueNames   Tinker residue names upon return
         @param tinkerBiotypes       Tinker biotypes upon return
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */ 
      
      static int getTinkerBiotypes( int numberOfAtoms, const t_topology* top, StringVector& tinkerAtomNames,
                                    StringVector& tinkerResidueNames, IntVector& tinkerBiotypes );
      
      /**---------------------------------------------------------------------------------------
      
         Get atomic numbers
      
         @param top          Gromac's topology struct
         @param atomicNumber output atomic numbers
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getAtomicNumbers( const t_topology* top, IntVector& atomicNumber );
      
      /**---------------------------------------------------------------------------------------
      
         Get OBC scale factors
      
         @param top          Gromac's topology struct
         @param scaleFactors output atomic numbers
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getObcScaleFactors( const t_topology* top, RealOpenMM* scaleFactors );
      
      /**---------------------------------------------------------------------------------------
      
         Get OBC scale factors
      
         @param top          Gromac's topology struct
         @param scaleFactors output atomic numbers
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getObcScaleFactors( const t_topology* top, RealOpenMMVector& scaleFactors );
      
      /**---------------------------------------------------------------------------------------
         
         Get Tinker biotype index given residue and atom names (Simbios)
         
         @param tinkerResidueName         Tinker residue name
         @param tinkerAtomName            Tinker atom name
      
         @return biotype if able to map names; otherwise return -1
      
         --------------------------------------------------------------------------------------- */
      
      static int getBiotypeGivenResidueAtomNames( const std::string& tinkerResidueName,
                                                  const std::string& tinkerAtomName );
      
      /**---------------------------------------------------------------------------------------
      
         Get Tinker biotype residue name given Gromacs residue name (Simbios)
      
         @param gromacsResidueName Gromac's residue name.
      
         @return AmoebaCommon::AmoebaNotSet if name not found in mapping
      
         --------------------------------------------------------------------------------------- */
      
      static std::string getTinkerBiotypeResidueNameGivenGromacsResidueName( const std::string gromacsResidueName );
      
      
      /**---------------------------------------------------------------------------------------
      
         Get Tinker-Gromacs residue name map (Simbios)
      
         The hash map is implemented as singleton
      
         @return StringMap mapping Gromac's residue names to Tinker biotype residue names
      
         --------------------------------------------------------------------------------------- */
      
      static StringMap* getTinkerGromacsResidueNameMap( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get residue name map ( Gromacs -> Tinker ) string (Simbios)
      
         @return string containing contents of residue map
      
         --------------------------------------------------------------------------------------- */
      
      static std::string getTinkerGromacsResidueNameMapString( void );
            
      /**---------------------------------------------------------------------------------------
      
         Get Tinker residue/atom name -> biotype map (Simbios)
      
         The hash map is implemented as singleton
      
         @return StringIntMap mapping Tinker residue_atomName -> Tinker biotype
      
         --------------------------------------------------------------------------------------- */
      
      static StringIntMap* getTinkerResidueAtomNameBiotypeMap( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get solvent radii from parameter file (Simbios) 
            
         @param numberOfAtoms       number of atoms
         @param parameterFileName   parameter file name
         @param top                 Gromacs topology data struct
         @param radii               array store Macromodel radii for each atom
         @param scaleFactor         scale factor
            
         @return SimTKOpenMMCommon::DefaultReturn unless paramter file not opened
                 in which case return SimTKOpenMMCommon::ErrorReturn
            
         --------------------------------------------------------------------------------------- */
      
      static int getMacroModelAtomicRadii( int numberOfAtoms, const std::string parameterFileName,
                                           const t_topology* top, RealOpenMM* radii, RealOpenMM scaleFactor = 1.0 );

      /**---------------------------------------------------------------------------------------
      
         Get solvent radii from parameter file (Simbios) 
            
         @param numberOfAtoms       number of atoms
         @param parameterFileName   parameter file name
         @param top                 Gromacs topology data struct
         @param radii               vector to store Macromodel radii for each atom
         @param scaleFactor         scale factor
            
         @return SimTKOpenMMCommon::DefaultReturn unless paramter file not opened
                 in which case return SimTKOpenMMCommon::ErrorReturn
            
         --------------------------------------------------------------------------------------- */
      
      static int getMacroModelAtomicRadii( int numberOfAtoms, const std::string parameterFileName,
                                           const t_topology* top, RealOpenMMVector& radii,
                                           RealOpenMM scaleFactor = 1.0 );
      
      /**---------------------------------------------------------------------------------------
      
         Get string containing atom types
            
         @param top                 Gromacs topology data struct
            
         @return string
            
         --------------------------------------------------------------------------------------- */
      
      static std::string getAtomTypesString( const t_topology* top );
      
      /**---------------------------------------------------------------------------------------
      
         Copy contents of Gromacs rvec array to RealOpenMM array
      
         If realArray == NULL on input, then memory is allocated -- callee is responsible for
         freeing memory:
      
         SimTKOpenMMUtilities::Xfree( "realArrayBlock",  __FILE__, __LINE__, realArray[0] );
         SimTKOpenMMUtilities::Xfree( "realArray",       __FILE__, __LINE__, realArray );
      
         @param numberOfEntries      number of entries in array
         @param gromacsArray         Gromac's array
         @param realArray            RealOpenMM** array (allocated if NULL on input)
         @param scaleFactor          scale factor
      
         @return realArray
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM** copyRvecArrayToRealOpenMMArray( int numberOfEntries, rvec* gromacsArray,
                                                          RealOpenMM** realArray,
                                                          RealOpenMM scaleFactor = 1.0 );
      
      /**---------------------------------------------------------------------------------------
      
         Get vectors of atom names, residue indices, and residue names from 
         Gromacs data structs (Simbios)
      
         @param top            Gromacs t_topolgy struct.
         @param residueNames   output residue names
         @param residueIndices output residue indices
         @param atomNames      output atom names
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getAtomResidueNames( const t_topology* top, StringVector& residueNames,
                                      IntVector& residueIndices, StringVector& atomNames );
      
   // ---------------------------------------------------------------------------------------

};
   
// ---------------------------------------------------------------------------------------

#endif // __SimTKOpenMMGromacsUtilities_H__
