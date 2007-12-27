
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

#ifdef NOTDEF

// class of shared, static utility methods

#include "SimTKOpenMMGromacsUtilities.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"

// fabs(), ...

#include <math.h>
#include <stdlib.h>

#define UseGromacsMalloc 1
#ifdef UseGromacsMalloc
extern "C" {
#include "smalloc.h" 
}
#endif

/* ---------------------------------------------------------------------------------------

   Find distances**2 from a given atom (Simbios)

   @param atomCoordinates     atom coordinates
   @param atomIndex           atom index to find distances from
   @param numberOfAtoms       number of atoms
   @param distances           array of distances squared on return; array size must be at least
                              numberOfAtoms

   @return distances

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getDistanceSquaredFromSpecifiedAtom( const rvec* atomCoordinates,
                                                                      int atomIndex,
                                                                      int numberOfAtoms, float* distances ){

   // ---------------------------------------------------------------------------------------

   float atomXyz[3];
   // static const char* methodName    = "\nSimTKOpenMMGromacsUtilities::getDistanceSquaredFromSpecifiedAtom";

   // ---------------------------------------------------------------------------------------

   for( int jj = 0; jj < 3; jj++ ){
      atomXyz[jj] = atomCoordinates[atomIndex][jj];
   }
      
   return getDistanceSquaredFromSpecifiedPoint( atomCoordinates, atomXyz,
                                                numberOfAtoms, distances );
}

/* ---------------------------------------------------------------------------------------

   Find distances**2 from a given point (Simbios)

   @param atomCoordinates     atom coordinates
   @param point               point to find distances from
   @param numberOfAtoms       number of atoms
   @param distances           array of distances squared on return; array size must be at least \n
                              numberOfAtoms

   @return SimTKOpenMMCommon::DefaultReturn;

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getDistanceSquaredFromSpecifiedPoint( const rvec* atomCoordinates,
                                                                       float* point,
                                                                       int numberOfAtoms,
                                                                       float* distances ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nSimTKOpenMMGromacsUtilities::getDistanceSquaredFromSpecifiedPoint";

   // ---------------------------------------------------------------------------------------

   memset( distances, 0, sizeof( float )*numberOfAtoms );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      for( int jj = 0; jj < 3; jj++ ){
         float diff = (point[jj] - atomCoordinates[ii][jj]);
         distances[ii] += diff*diff;
      }
   }
      
   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------

   Get atom name from top data struct

   @param atomIndex           atom index
   @param outputAtomName      output atom name
   @param top                 GMX t_topology struct

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getAtomNameGivenAtomIndex( int atomIndex, char* outputAtomName,
                                                            const t_topology* top ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nSimTKOpenMMGromacsUtilities::getAtomNameGivenAtomIndex";

   // ---------------------------------------------------------------------------------------

   char*** atomNames       = top->atoms.atomname;
   const char* atomName    =  (*(atomNames[atomIndex])) == NULL || strlen( (*(atomNames[atomIndex])) ) < 1 ||
                              strlen( (*(atomNames[atomIndex])) ) > 100 ? "NA" : (*(atomNames[atomIndex]));

#ifdef WIN32
   (void) strcpy_s( outputAtomName, sizeof( atomName ), atomName );
#else
   (void) strcpy( outputAtomName, atomName );
#endif

   return SimTKOpenMMCommon::DefaultReturn;

}

/* ---------------------------------------------------------------------------------------

   Get residue name from top data struct given atom index

   @param atomIndex           atom index
   @param top                 GMX t_topology struct
   @param outputResidueName   output residue name (assume enough memory has been allocated)
   @param outputResidueIndex  if not null, then *outputResidueIndex is residue index

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getResidueNameGivenAtomIndex( int atomIndex, const t_topology* top,
                                                               char* outputResidueName,
                                                               int* outputResidueIndex ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMGromacsUtilities::getResidueNameGivenAtomIndex";

   // ---------------------------------------------------------------------------------------

   char*** residueNames = top->atoms.resname;
   int residueIndex     = top->atoms.atom[atomIndex].resnr;

   if( outputResidueIndex != NULL ){
      *outputResidueIndex = residueIndex;
   }

   const char* residueName =  (*(residueNames[residueIndex])) == NULL          ||
                              strlen( (*(residueNames[residueIndex])) ) < 1 ||
                              strlen( (*(residueNames[residueIndex])) ) > 100 ? "NA" : (*(residueNames[residueIndex]));

#ifdef WIN32
   (void) strcpy_s( outputResidueName, sizeof( residueName ), residueName );
#else
   (void) strcpy( outputResidueName, residueName );
#endif

   return SimTKOpenMMCommon::DefaultReturn;

   // ---------------------------------------------------------------------------------------
}

/* ---------------------------------------------------------------------------------------

   Get atom name from top data struct

   @param atomIndex           atom index
   @param top                 GMX t_topology struct
   @param buffer              output buffer (enough space should have been reserved) 
   @param maxAtoms            max number of atoms for this run (may change -- used mainly
                              to keep from reallocating cache array)
   @param tab                 tab spacing

   @return SimTKOpenMMCommon::DefaultReturn

   Atom info is cached in atomIdStrings to speed things up.
   The cache memory may be freed by calling the method w/
   atomIndex == -1

   @see freeArrayOfStrings()

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex( int atomIndex, const t_topology* top,
                                                                int sizeOfBuffer, char* buffer, 
                                                                int maxAtoms, unsigned int tab ){

   // ---------------------------------------------------------------------------------------

   static int maxAtomIndex         = -1;
   static char** atomIdStrings     = NULL;

   char atomName[32];
   char residueName[32];

   int residueIndex;

   // static const char* methodName = "\nSimTKOpenMMGromacsUtilities::getAtomIdStringGivenAtomIndex";

   // ---------------------------------------------------------------------------------------

   // free cache memory if allocated

   if( atomIndex == -1 ){
      // (void) fprintf( stdout, "SimTKOpenMMGromacsUtilities: called getAtomIdStringGivenAtomIndex to delete cached strings %d.", maxAtomIndex );
      if( maxAtomIndex > 0 && atomIdStrings != NULL ){
         SimTKOpenMMUtilities::freeArrayOfStrings( maxAtomIndex, atomIdStrings );
         atomIdStrings = NULL;
         maxAtomIndex  = -1;
      }         
      return SimTKOpenMMCommon::DefaultReturn;
   }         

   // allocate cache memory?

   if( maxAtoms > 0 && (maxAtomIndex == -1 || maxAtoms > maxAtomIndex) ){
      if(  maxAtoms > maxAtomIndex && maxAtomIndex > 0 ){
         SimTKOpenMMUtilities::freeArrayOfStrings( maxAtomIndex, atomIdStrings );
      }
      maxAtomIndex = maxAtoms + 1;

//    atomIdStrings = (char**) malloc( maxAtomIndex*sizeof( char* ) ); 
      atomIdStrings = (char**) SimTKOpenMMUtilities::Xmalloc( "atomIdStrings", __FILE__, __LINE__, maxAtomIndex*sizeof( char* ) );
      memset( atomIdStrings, 0, maxAtomIndex*sizeof( char* ) );

   }

   // if id is cached, return it

   if( atomIndex < maxAtomIndex && atomIdStrings[atomIndex] != NULL ){ 
#ifdef WIN32
      (void) strcpy_s( buffer, sizeof( atomIdStrings[atomIndex]), atomIdStrings[atomIndex] );
#else
      (void) strcpy( buffer, atomIdStrings[atomIndex] );
#endif
      return SimTKOpenMMCommon::DefaultReturn;
   }

   // not cached -- assemble info

   getAtomNameGivenAtomIndex( atomIndex, atomName, top );
   getResidueNameGivenAtomIndex( atomIndex, top, residueName, &residueIndex );

#ifdef WIN32
   (void) sprintf_s( buffer, sizeOfBuffer, "%s_%d %s", residueName, residueIndex, atomName );
#else
   (void) sprintf( buffer, "%s_%d %s", residueName, residueIndex, atomName );
#endif

   // tab string

   if( tab > 0 && strlen( buffer ) < tab ){
      SimTKOpenMMUtilities::tabStringInPlace( buffer, tab );
   }

   // cache info if atomIdStrings array is allocated

   if( atomIndex < maxAtomIndex && atomIdStrings ){

      if( atomIdStrings[atomIndex] != NULL ){
         SimTKOpenMMUtilities::Xfree( "atomIdStrings", __FILE__, __LINE__,  atomIdStrings[atomIndex] );
      }
      unsigned int bufferSz  = (unsigned int) sizeof( char );
      bufferSz              *= (unsigned int) (strlen( buffer ) + 1);
      atomIdStrings[atomIndex] = (char*) SimTKOpenMMUtilities::Xmalloc( "atomIdStrings[atomIndex]", __FILE__, __LINE__, bufferSz );
#ifdef WIN32
      (void) strcpy_s( atomIdStrings[atomIndex], sizeof( buffer ), buffer );
#else
      (void) strcpy( atomIdStrings[atomIndex], buffer );
#endif
   }

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get (1-2) bonds (Simbios) 

   @param maxAtoms            max number of atoms
   @param IntSetVector        vector of integer sets
   @param top                 Gromacs t_topolgy struct

   @return 0 if no errors or
   return x, where x is the number of errors encountered

   Upon return covalentBonds[i] = set of atom indices that are covalent partners

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getCovalentBondIndices( int maxAtoms, IntSetVector& covalentBonds,
                                                         const t_topology* top ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nSimTKOpenMMGromacsUtilities::getCovalentBonds";

   // ---------------------------------------------------------------------------------------

   // indices

   int idefArrayIndex;
   int offset          = 3;
   int atomIndexOffset = 1;

   // load 1-2 bonds

   // t_iparams* params   = top->idef.iparams;
   int errors          = 0;

   covalentBonds.resize( maxAtoms + 1 );

   for( int jj = 0; jj < 2; jj++ ){
 
      idefArrayIndex      = jj ? F_SHAKE : F_BONDS;
      t_iatom* atoms      = top->idef.il[idefArrayIndex].iatoms;

      for( int ii = 0; ii < top->idef.il[idefArrayIndex].nr; ii += offset ){
   
         // int type  = (int) atoms[ ii ];
         int atomI = (int) atoms[ ii + atomIndexOffset ];
         int atomJ = (int) atoms[ ii + atomIndexOffset + 1 ];
   
         // validate indices
   
         if( atomI >= maxAtoms || atomI < 0 ){
            std::stringstream message;
            message << methodName;
            message << " atom index=" << atomI << " (Gromacs index=" << ii << ") too large: max=" << maxAtoms;
            SimTKOpenMMLog::printMessage( message );
            errors++;
            atomI = -1;
         }
   
         if( atomJ >= maxAtoms || atomJ < 0 ){
            std::stringstream message;
            message << methodName;
            message << " atom index=" << atomJ << " (Gromacs index=" << ii << ") too large: max=" << maxAtoms;
            SimTKOpenMMLog::printMessage( message );
            errors++;
            atomJ = -1;
         }
   
         // RealOpenMM bondLength = params[type].harmonic.rA; 
   
         if( atomI >= 0 && atomJ >= 0 ){
            covalentBonds[atomI].insert( atomJ );
            covalentBonds[atomJ].insert( atomI );
         }
   
      }
   }

   // waters?

   errors += SimTKOpenMMGromacsUtilities::getSettleCovalentBondIndices( maxAtoms, covalentBonds, top );

   return errors;
}

/**---------------------------------------------------------------------------------------

   Add SETTLE stretch (1-2) bonds (Simbios) 

   @param maxAtoms            max number of atoms
   @param IntSetVector        vector of integer sets
   @param top                 Gromacs t_topolgy struct

   @return SimTKOpenMMCommon::DefaultReturn if no errors or
   return x, where x is the number of errors encountered

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getSettleCovalentBondIndices( int maxAtoms,
                                                               IntSetVector& covalentBonds,
                                                               const t_topology* top ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "SimTKOpenMMGromacsUtilities::getSettleCovalentBondIndices";

   // ---------------------------------------------------------------------------------------

   // indices

   int idefArrayIndex  = F_SETTLE;
   int offset          = 2;
   int atomIndexOffset = 1;

   // load bonds via SETTLE parameter

   t_iatom* atoms      = top->idef.il[idefArrayIndex].iatoms;
   t_iparams* params   = top->idef.iparams;
   int errors          = 0;

   covalentBonds.resize( maxAtoms + 1 );

   for( int ii = 0; ii < top->idef.il[idefArrayIndex].nr; ii += offset ){

      //int type            = (int) atoms[ ii ];
      int atomI           = (int) atoms[ ii + atomIndexOffset ];
      //RealOpenMM OHBondLength  = params[type].harmonic.rA; 
      //RealOpenMM HHBondLength  = params[type].harmonic.krA;
      //RealOpenMM bondAngle     = 1.0f - ( (HHBondLength*HHBondLength)/ (2.0f*OHBondLength*OHBondLength));
      //      bondAngle     = acosf( bondAngle );

      // validate indices

      if( (atomI + 2) >= maxAtoms || atomI < 0 ){
         std::stringstream message;
         message << methodName;
         message << " atom index=" << atomI << " (Gromacs index=" << ii << ") too large: max=" << maxAtoms;
         SimTKOpenMMLog::printMessage( message );
         errors++;
         atomI = -1;
      }

      // add 2 O-H bonds

      if( atomI >= 0 ){

         covalentBonds[atomI].insert( atomI + 1 );
         covalentBonds[atomI+1].insert( atomI );

         covalentBonds[atomI].insert( atomI + 2 );
         covalentBonds[atomI+2].insert( atomI );
      }

   }

   return errors;
}

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

int SimTKOpenMMGromacsUtilities::writeTinkerXyzFile( int numberOfAtoms, const rvec* atomCoordinates, 
                                                     const std::string& header, const std::string& xyzFileName,
                                                     const t_topology* top ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName        = "\nSimTKOpenMMGromacsUtilities::writeXyzFile";

   // ---------------------------------------------------------------------------------------

   // get covalent bonds

   IntSetVector covalentBonds;
   SimTKOpenMMGromacsUtilities::getCovalentBondIndices( numberOfAtoms, covalentBonds, top );

   // get Tinker biotypes

   StringVector tinkerAtomNames;
   StringVector tinkerResidueNames;
   IntVector tinkerBiotypes;
   SimTKOpenMMGromacsUtilities::getTinkerBiotypes( numberOfAtoms, top, tinkerAtomNames, tinkerResidueNames, tinkerBiotypes );

   // open file

   FILE* xyzFile = NULL;
#ifdef WIN32
   fopen_s( &xyzFile, xyzFileName.c_str(), "w" );
#else
   xyzFile = fopen( xyzFileName.c_str(), "w" );
#endif

   if( xyzFile != NULL ){
      std::stringstream message;
      message << methodName;
      message << " Opened file=<" << xyzFileName.c_str() << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName;
      message << "  could not open file=<" << xyzFileName.c_str() << "> -- abort output.";
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;
   }

   // first line

   (void) fprintf( xyzFile, "%d %s\n", numberOfAtoms, header.c_str() );

/*
  1232  CHROMOSOMAL PROTEIN                     02-JAN-87   1UBQ 
     1  N3    27.340000   24.430000    2.614000   472     2     5     6     7    
     2  CT    26.266000   25.413000    2.842000   473     1     3     8     9    
     3  C     26.913000   26.639000    3.531000   474     2     4    20   
*/

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // scale coordinates by 10 for Angstrom -> nanometer conversion

      (void) fprintf( xyzFile, "%6d  %-4s %16.9f %16.9f %16.9f %6d ", atomI+1, tinkerAtomNames[atomI].c_str(),
                      10.0f*atomCoordinates[atomI][0], 10.0f*atomCoordinates[atomI][1], 10.0f*atomCoordinates[atomI][2],
                      tinkerBiotypes[atomI] );

      // include 1-2 bonds

      IntSet atomCovalentBonds = covalentBonds[atomI]; 
      for( IntSetCI kk = atomCovalentBonds.begin(); kk != atomCovalentBonds.end(); kk++ ){
         (void) fprintf( xyzFile, "%6d ", (*kk+1) );
      }
      (void) fprintf( xyzFile, "\n" );
   }
   (void) fflush( xyzFile );
   (void) fclose( xyzFile );

   // diagnostics

   std::stringstream message;
   message << methodName;
   message << " closed file=<" << xyzFileName.c_str() << ">.";
   SimTKOpenMMLog::printMessage( message );

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get Tinker biotypes (Simbios)

   @param numberOfAtoms        number of atoms
   @param top                  Gromacs topology struct
   @param tinkerAtomNames      Tinker atom names upon return
   @param tinkerResidueNames   Tinker residue names upon return
   @param tinkerBiotypes       Tinker biotypes upon return

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getTinkerBiotypes( int numberOfAtoms, const t_topology* top,
                                                    StringVector& tinkerAtomNames,
                                                    StringVector& tinkerResidueNames,
                                                    IntVector& tinkerBiotypes ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName        = "\nSimTKOpenMMGromacsUtilities::getTinkerBiotypes";

   static const unsigned int bufferSz   = 128;
   char atomName[128];

   // ---------------------------------------------------------------------------------------

   tinkerAtomNames.resize( numberOfAtoms );
   tinkerResidueNames.resize( numberOfAtoms );
   tinkerBiotypes.resize( numberOfAtoms );

   // Gromacs array of atom and residue names

   char*** atomNames     = top->atoms.atomname;

   // char*** atomType      = top->atoms.atomtype;

   int numberOfResidues  = top->atoms.nres;
   char*** residueNames  = top->atoms.resname;

   // for each atom, set names, indices, ...

   int residueI          = 0;    
   int firstAtomNameL    = -1;   
   char* firstAtomName   = NULL; 
   int firstWaterAtom    = -1;   

   std::string residueName;
   std::string tinkerResidueName;
   
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // check if hitting new residue
      // save first atom in list to allow for check
      // when another atom is hit

#ifdef WIN32
       (void) strcpy_s( atomName, bufferSz, *(atomNames[atomI]) );
#else
       (void) strcpy( atomName, *(atomNames[atomI]) );
#endif

      if( !firstAtomName ){

         firstAtomName     = *(atomNames[atomI]);
         firstAtomNameL    = (int) strlen( *(atomNames[atomI]) );
         residueName       = std::string( *(residueNames[residueI]) );
         tinkerResidueName = SimTKOpenMMGromacsUtilities::getTinkerBiotypeResidueNameGivenGromacsResidueName( residueName );

      } else if( ((int) strlen( *(atomNames[atomI]) ) == firstAtomNameL && !strcmp( atomName, firstAtomName )) ||
                 !strcmp( atomName, "OW" ) ){
         residueI++;
         residueName       = std::string( *(residueNames[residueI]) );
         residueI          = (residueI >= numberOfResidues) ? numberOfResidues - 1 : residueI;
         if( firstWaterAtom < 0 && !strcmp( atomName, "OW" ) ){
            firstWaterAtom = atomI;
         }
         tinkerResidueName = SimTKOpenMMGromacsUtilities::getTinkerBiotypeResidueNameGivenGromacsResidueName( residueName );
      }

/*
if( residueI >= 0 ){
   std::stringstream message;
   message << "\n" << residueI << ". " << *(residueNames[residueI]) << " Tinker=<" << tinkerResidueName << "> <" << *(atomNames[atomI]) << ">";
   SimTKOpenMMLog::printMessage( message );
} */

      // OC1, OC2 -> OXT at C-terminus

      if( atomName[0] == 'O' && atomName[1] == 'C' ){

         atomName[1] = 'X';
         atomName[2] = 'T';
         atomName[3] = '\0';

      // H1, H2, H3 -> HN at N-terminus

      } else if( residueI == 0 && atomName[0] == 'H' && (atomName[1] == '1' || atomName[1] == '2' || atomName[1] == '3') ){
         atomName[1] = 'N';

      // solvent atoms

      } else if( atomName[1] == 'W' && ( atomName[0] == 'O' || atomName[0] == 'H' ) ){
         atomName[1] = '\0';
      }

      int bioType =SimTKOpenMMGromacsUtilities::getBiotypeGivenResidueAtomNames( tinkerResidueName, std::string( atomName ) );

      // handle cases like "NLEU" where type is in "LEU"

      if( bioType < 0 && residueName.length() == 4 ){ 
         std::string removeFirstChar        = residueName.substr( 1, 4 ); 
         std::string localTinkerResidueName = SimTKOpenMMGromacsUtilities::getTinkerBiotypeResidueNameGivenGromacsResidueName( removeFirstChar );
         bioType =SimTKOpenMMGromacsUtilities::getBiotypeGivenResidueAtomNames( localTinkerResidueName, std::string( atomName ) );
      }

      tinkerAtomNames[atomI]    = atomName;
      tinkerResidueNames[atomI] = residueName;
      tinkerBiotypes[atomI]     = bioType;

   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get Tinker biotype index given residue and atom names (Simbios)

   @param tinkerResidueName         Tinker residue name
   @param tinkerAtomName				Tinker atom name

   @return biotype if able to map names; otherwise return -1

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getBiotypeGivenResidueAtomNames( const std::string& tinkerResidueName,
                                                                  const std::string& tinkerAtomName ){

   // ---------------------------------------------------------------------------------------

   char localAtomName[256];

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getBiotypeGivenResidueAtomNames";

   // ---------------------------------------------------------------------------------------

   // build residue/atom name

   std::string fullName = tinkerResidueName;
   fullName.append( "_" );
   fullName.append( tinkerAtomName );

   StringIntMap* nameMap = SimTKOpenMMGromacsUtilities::getTinkerResidueAtomNameBiotypeMap( );

   StringIntMapCI reference = nameMap->find( fullName );
   if(  reference != nameMap->end( ) ){
      return (*reference).second;
   }

   // if atom name ends with two digits -> remove last digit
   // HD12 -> HD1

#ifdef WIN32
   (void) sprintf_s( localAtomName, 256, "%s", tinkerAtomName.c_str() );
#else
   (void) sprintf( localAtomName, "%s", tinkerAtomName.c_str() );
#endif
   unsigned int lastChar = (unsigned int) (strlen( localAtomName ) - 1);
   if( lastChar > 2 && isdigit( localAtomName[lastChar] ) && isdigit( localAtomName[lastChar-1] ) ){
      localAtomName[lastChar] = '\0';
      lastChar--;
   }

   // 'H' -> 'HN' (skipping H in water)

   if( localAtomName[0] == 'H' && lastChar == 0 && strcmp( tinkerResidueName.c_str(), "AMOEBA_Water" ) ){
      localAtomName[1] = 'N';
      localAtomName[2] = '\0';
   }

   // Solvent: OW -> O && HW1/HW2 -> H

   if( lastChar > 0 && localAtomName[1] == 'W' ){
      localAtomName[2] = '\0';
   }

   std::stringstream messageX;

   std::string modifiedAtomName( localAtomName );
   fullName = tinkerResidueName;
   fullName.append( "_" );
   fullName.append( modifiedAtomName );

//messageX << " " << modifiedAtomName;

   // check if _tinkerResidueBiotypeParameterMap contains Tinker atom name
   // if last name ends in a digit and name was not found, try removing trailing digit

// (void) fprintf( log, "%s check tinkerBiotypeParameterMap (residue=%s) Tinker atom type=<%s>; input=<%s>",
//                methodName.c_str(), tinkerResidueName->c_str(), modifiedAtomName.c_str(), tinkerAtomName->c_str() );
// (void) fflush( log );

   reference = nameMap->find( fullName );
   if( reference == nameMap->end( ) ){

      if( isdigit( localAtomName[lastChar] ) ){
         localAtomName[lastChar] = '\0';
         modifiedAtomName = std::string( localAtomName );

         fullName = tinkerResidueName;
         fullName.append( "_" );
         fullName.append( modifiedAtomName );

//messageX << " " << modifiedAtomName;
         reference = nameMap->find( fullName);
      }
   }

   if( reference == nameMap->end( ) ){

      // no message for atoms w/ residue names of length 4
      // program will try agian w/ N/C prefix removed from name (e.g., NLEU -> LEU )
      // printing message will confuse user 

      if( tinkerResidueName.size() < 4 ){
         std::stringstream message;
         message <<  methodName.c_str() << " Missing Tinker residue=<" <<
                     tinkerResidueName.c_str() << "> <" <<  tinkerAtomName.c_str() << ">.";
//message << " " << messageX.str();
         SimTKOpenMMLog::printMessage( message );
      }
      return -1;
   } else {
      return (*reference).second;
   }
}

/**---------------------------------------------------------------------------------------

   Get Tinker biotype residue name given Gromacs residue name (Simbios)

   @param gromacsResidueName Gromac's residue name.

   @return AmoebaCommon::AmoebaNotSet if name not found in mapping

   --------------------------------------------------------------------------------------- */

std::string SimTKOpenMMGromacsUtilities::getTinkerBiotypeResidueNameGivenGromacsResidueName(
            const std::string gromacsResidueName ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getTinkerBiotypeResidueNameGivenGromacsResidueName";

   // ---------------------------------------------------------------------------------------

   StringMap* residueNameMap = getTinkerGromacsResidueNameMap();
   if( residueNameMap != NULL ){
      StringMapI mapIterator   = residueNameMap->find( gromacsResidueName );
      if( mapIterator != residueNameMap->end() ){
         return (*mapIterator).second;
      }
   }

   return SimTKOpenMMCommon::NotSet;
}

/**---------------------------------------------------------------------------------------

   Get atomic numbers

   @param top          Gromac's topology struct
   @param atomicNumber output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getAtomicNumbers( const t_topology* top, IntVector& atomicNumber ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getAtomicNumbers";

   // ---------------------------------------------------------------------------------------

   atomicNumber.resize( top->atoms.nr );
   for( int atomI = 0; atomI < top->atoms.nr; atomI++ ){

      int element = 0;
      RealOpenMM mass = (RealOpenMM) top->atoms.atom[atomI].m;

      if( mass < (RealOpenMM) 1.2 && mass >= (RealOpenMM) 1.0 ){          // hydrogen
         element = 1; 
      } else if( mass > (RealOpenMM) 11.8 && mass < (RealOpenMM) 12.2 ){  // carbon
         element = 6; 
      } else if( mass > (RealOpenMM) 14.0 && mass < (RealOpenMM) 15.0 ){  // nitrogen
         element = 7; 
      } else if( mass > (RealOpenMM) 15.5 && mass < (RealOpenMM) 16.5 ){  // oxygen
         element = 8; 
      } else if( mass > (RealOpenMM) 31.5 && mass < (RealOpenMM) 32.5 ){  // sulphur
         element = 16;
      } else if( mass > (RealOpenMM) 29.5 && mass < (RealOpenMM) 30.5 ){  // phosphorus
         element = 15;
      } else {
         std::stringstream message;
         char*** atomNames = top->atoms.atomname;
         message << methodName;
         message << " Warning: mass for atom=<" << (*(atomNames[atomI])) << "> mass=" << mass << "> not recognized.";
         SimTKOpenMMLog::printMessage( message );
      }

      atomicNumber[atomI] = element;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get OBC scale factors

   @param top          Gromac's topology struct
   @param scaleFactors output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getObcScaleFactors( const t_topology* top, RealOpenMM* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getObcScaleFactors";

   // ---------------------------------------------------------------------------------------

   RealOpenMMVector scaleFactorV;
   int status = SimTKOpenMMGromacsUtilities::getObcScaleFactors( top, scaleFactorV );
   for( int ii = 0; ii < top->atoms.nr; ii++ ){
      scaleFactors[ii] = scaleFactorV[ii];
   }
   return status;
}

/**---------------------------------------------------------------------------------------

   Get OBC scale factors

   @param top          Gromac's topology struct
   @param scaleFactors output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getObcScaleFactors( const t_topology* top, 
                                                     RealOpenMMVector& scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getObcScaleFactors";

   // ---------------------------------------------------------------------------------------

   scaleFactors.resize( top->atoms.nr );

   for( int atomI = 0; atomI < top->atoms.nr; atomI++ ){

      RealOpenMM scaleFactor = (RealOpenMM) 0.8;
      RealOpenMM mass        = (RealOpenMM) top->atoms.atom[atomI].m;

      if ( mass < (RealOpenMM) 1.2 && mass >= (RealOpenMM) 1.0 ){        // hydrogen
         scaleFactor  = (RealOpenMM) 0.85; 
      } else if( mass > (RealOpenMM) 11.8 && mass < (RealOpenMM) 12.2 ){ // carbon
         scaleFactor  = (RealOpenMM) 0.72; 
      } else if( mass > (RealOpenMM) 14.0 && mass < (RealOpenMM) 15.0 ){ // nitrogen
         scaleFactor  = (RealOpenMM) 0.79;
      } else if( mass > (RealOpenMM) 15.5 && mass < (RealOpenMM) 16.5 ){ // oxygen
         scaleFactor  = (RealOpenMM) 0.85; 
      } else if( mass > (RealOpenMM) 31.5 && mass < (RealOpenMM) 32.5 ){ // sulphur
         scaleFactor  = (RealOpenMM) 0.96;
      } else if( mass > (RealOpenMM) 29.5 && mass < (RealOpenMM) 30.5 ){ // phosphorus
         scaleFactor  = (RealOpenMM) 0.86;
      } else {
         std::stringstream message;
         char*** atomNames = top->atoms.atomname;
         message << methodName;
         message << " Warning: mass for atom=<" << (*(atomNames[atomI])) << "> mass=" << mass << "> not recognized.";
         SimTKOpenMMLog::printMessage( message );
      }

      scaleFactors[atomI] = scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get Tinker-Gromacs residue name map (Simbios)

   The hash map is implemented as singleton

   @return StringMap mapping Gromac's residue names to Tinker biotype residue names

   --------------------------------------------------------------------------------------- */

StringMap* SimTKOpenMMGromacsUtilities::getTinkerGromacsResidueNameMap( void ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getTinkerBiotypeResidueNameGivenGromacsResidueName";

   // CYS2 == CYX ?

   static char* nameMap[][3] = {
                                  { "ALA",  "Alanine",              NULL       },
                                  { "ARG",  "Arginine",             NULL       },
                                  { "ASN",  "Asparagine",           NULL       },
                                  { "ASP",  "Aspartic_Acid",        NULL       },
                                  { "CYS",  "Cysteine_(-SH)",       NULL       },
                                  { "CYX",  "Cystine_(-SS-)",       NULL       },
                                  { "CYS2", "Cystine_(-SS-)",       NULL       },
                                  { "GLU",  "Glutamic_Acid",        NULL       },
                                  { "GLN",  "Glutamine",            NULL       },
                                  { "GLY",  "Glycine",              NULL       },
                                  { "HIP",  "Histidine_(+)",        NULL       },
                                  { "HID",  "Histidine_(HD)",       NULL       },
                                  { "HIE",  "Histidine_(HE)",       NULL       },
                                  { "ILE",  "Isoleucine",           NULL       },
                                  { "LEU",  "Leucine",              NULL       },
                                  { "NLE",  "Lysine",               "LYS"      },
                                  { "LYS",  "Lysine",               NULL       },
                                  { "LYN",  "Lysine",               "LYS"      },
                                  { "LYSH", "Lysine",               "LYS"      },
                                  { "LYP",  "Lysine",               "LYS"      },
                                  { "MET",  "Methionine",           NULL       },
                                  { "AIB",  "MethylAlanine_(AIB)",  NULL       },
                                  { "PHE",  "Phenylalanine",        NULL       },
                                  { "PRO",  "Proline",              NULL       },
                                  { "SER",  "Serine",               NULL       },
                                  { "THR",  "Threonine",            NULL       },
                                  { "TRP",  "Tryptophan",           NULL       },
                                  { "TYR",  "Tyrosine",             NULL       },
                                  { "VAL",  "Valine",               NULL       },
                                  { "SOL",  "AMOEBA_Water",         NULL       },
                                  {  NULL,  NULL,                   NULL       }
                               };


   static const int bufferSz = 128;
   char auxiliaryName[2][bufferSz];

   static StringMap* residueNameMap = NULL;

   // ---------------------------------------------------------------------------------------

   // residueNameMap is singleton

   if( residueNameMap == NULL ){

      residueNameMap = new StringMap();

      // unsigned int initialMapSize = sizeof( nameMap )/(2*sizeof( char* ));

      int index = 0;
      while( nameMap[index][0] != NULL ){

         (*residueNameMap)[nameMap[index][0]] = nameMap[index][1];

         // add C & N terminal residues
         // e.g., NMET -> N-Terminal_MET

         if( nameMap[index][2] != NULL ){
#ifdef WIN32
            (void) sprintf_s( auxiliaryName[1], bufferSz, "N-Terminal_%s", nameMap[index][2] );
#else
            (void) sprintf( auxiliaryName[1], "N-Terminal_%s", nameMap[index][2] );
#endif
         } else {

#ifdef WIN32
            (void) sprintf_s( auxiliaryName[1], bufferSz, "N-Terminal_%s", nameMap[index][0] );
#else
            (void) sprintf( auxiliaryName[1], "N-Terminal_%s", nameMap[index][0] );
#endif
         }

#ifdef WIN32
         (void) sprintf_s( auxiliaryName[0], bufferSz, "N%s", nameMap[index][0] );
#else
         (void) sprintf( auxiliaryName[0], "N%s", nameMap[index][0] );
#endif

         (*residueNameMap)[auxiliaryName[0]] = auxiliaryName[1];

#ifdef WIN32
         (void) sprintf_s( auxiliaryName[1], bufferSz, "C-Terminal_%s", nameMap[index][0] );
         (void) sprintf_s( auxiliaryName[0], bufferSz, "C%s", nameMap[index][0] );
#else
         (void) sprintf( auxiliaryName[1], "C-Terminal_%s", nameMap[index][0] );
         (void) sprintf( auxiliaryName[0], "C%s", nameMap[index][0] );
#endif

         (*residueNameMap)[auxiliaryName[0]] = auxiliaryName[1];

         index++;
      }

   }

   return residueNameMap;
}

/**---------------------------------------------------------------------------------------

   Get residue name map ( Gromacs -> Tinker ) string (Simbios)

   @return string containing contents of residue map

   --------------------------------------------------------------------------------------- */

std::string SimTKOpenMMGromacsUtilities::getTinkerGromacsResidueNameMapString( void ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getTinkerGromacsResidueNameMapString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   StringMap* residueNameMap = SimTKOpenMMGromacsUtilities::getTinkerGromacsResidueNameMap();

   message << std::endl;
   message << "ResidueNameMap" << std::endl;
   int count=1;
   for( StringMapCI ii = residueNameMap->begin(); ii != residueNameMap->end(); ii++ ){
      message << count++ << " <";
      message << (*ii).first.c_str()  << "> <";
      message << (*ii).second.c_str() << ">";
      message << std::endl;
   }

   return message.str();
}

/**---------------------------------------------------------------------------------------

   Get Tinker residue/atom name -> biotype map (Simbios)

   The hash map is implemented as singleton

   @return StringIntMap mapping Tinker residue_atomName -> Tinker biotype

   --------------------------------------------------------------------------------------- */

StringIntMap* SimTKOpenMMGromacsUtilities::getTinkerResidueAtomNameBiotypeMap( void ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getTinkerResidueAtomNameBiotypeMap";

   // entries for map generated by TinkerBiotypeParameterSet::getResidueAtomMapBiotypeString()
   // in TinkerParameters::_readFile()

   static char* nameMap[][2] = {
                     { "AMOEBA_Water_H",               "203" },
                     { "AMOEBA_Water_O",               "202" },
                     { "Acetyl_N-Terminus_C",          "182" },
                     { "Acetyl_N-Terminus_CH3",        "180" },
                     { "Acetyl_N-Terminus_H",          "181" },
                     { "Acetyl_N-Terminus_O",          "183" },
                     { "Alanine_C",                    "9" },
                     { "Alanine_CA",                   "8" },
                     { "Alanine_CB",                   "13" },
                     { "Alanine_HA",                   "12" },
                     { "Alanine_HB",                   "14" },
                     { "Alanine_HN",                   "10" },
                     { "Alanine_N",                    "7" },
                     { "Alanine_O",                    "11" },
                     { "Amide_C-Terminus_HN",          "185" },
                     { "Amide_C-Terminus_N",           "184" },
                     { "Arginine_C",                   "9" },
                     { "Arginine_CA",                  "8" },
                     { "Arginine_CB",                  "169" },
                     { "Arginine_CD",                  "173" },
                     { "Arginine_CG",                  "171" },
                     { "Arginine_CZ",                  "177" },
                     { "Arginine_HA",                  "12" },
                     { "Arginine_HB",                  "170" },
                     { "Arginine_HD",                  "174" },
                     { "Arginine_HE",                  "176" },
                     { "Arginine_HG",                  "172" },
                     { "Arginine_HH",                  "179" },
                     { "Arginine_HN",                  "10" },
                     { "Arginine_N",                   "7" },
                     { "Arginine_NE",                  "175" },
                     { "Arginine_NH",                  "178" },
                     { "Arginine_O",                   "11" },
                     { "Asparagine_C",                 "9" },
                     { "Asparagine_CA",                "8" },
                     { "Asparagine_CB",                "138" },
                     { "Asparagine_CG",                "140" },
                     { "Asparagine_HA",                "12" },
                     { "Asparagine_HB",                "139" },
                     { "Asparagine_HD2",               "143" },
                     { "Asparagine_HN",                "10" },
                     { "Asparagine_N",                 "7" },
                     { "Asparagine_ND2",               "142" },
                     { "Asparagine_O",                 "11" },
                     { "Asparagine_OD1",               "141" },
                     { "Aspartic_Acid_C",              "9" },
                     { "Aspartic_Acid_CA",             "8" },
                     { "Aspartic_Acid_CB",             "128" },
                     { "Aspartic_Acid_CG",             "130" },
                     { "Aspartic_Acid_HA",             "12" },
                     { "Aspartic_Acid_HB",             "129" },
                     { "Aspartic_Acid_HN",             "10" },
                     { "Aspartic_Acid_N",              "7" },
                     { "Aspartic_Acid_O",              "11" },
                     { "Aspartic_Acid_OD",             "131" },
                     { "C-Terminal_AIB_C",             "999" },
                     { "C-Terminal_AIB_CA",            "999" },
                     { "C-Terminal_AIB_HN",            "999" },
                     { "C-Terminal_AIB_N",             "999" },
                     { "C-Terminal_AIB_OXT",           "999" },
                     { "C-Terminal_ALA_C",             "192" },
                     { "C-Terminal_ALA_CA",            "8" },
                     { "C-Terminal_ALA_HA",            "12" },
                     { "C-Terminal_ALA_HN",            "10" },
                     { "C-Terminal_ALA_N",             "7" },
                     { "C-Terminal_ALA_OXT",           "193" },
                     { "C-Terminal_ARG_C",             "192" },
                     { "C-Terminal_ARG_CA",            "8" },
                     { "C-Terminal_ARG_HA",            "12" },
                     { "C-Terminal_ARG_HN",            "10" },
                     { "C-Terminal_ARG_N",             "7" },
                     { "C-Terminal_ARG_OXT",           "193" },
                     { "C-Terminal_ASN_C",             "192" },
                     { "C-Terminal_ASN_CA",            "8" },
                     { "C-Terminal_ASN_HA",            "12" },
                     { "C-Terminal_ASN_HN",            "10" },
                     { "C-Terminal_ASN_N",             "7" },
                     { "C-Terminal_ASN_OXT",           "193" },
                     { "C-Terminal_ASP_C",             "192" },
                     { "C-Terminal_ASP_CA",            "8" },
                     { "C-Terminal_ASP_HA",            "12" },
                     { "C-Terminal_ASP_HN",            "10" },
                     { "C-Terminal_ASP_N",             "7" },
                     { "C-Terminal_ASP_OXT",           "193" },
                     { "C-Terminal_CYS_(-SH)_C",       "192" },
                     { "C-Terminal_CYS_(-SH)_CA",      "44" },
                     { "C-Terminal_CYS_(-SH)_HA",      "12" },
                     { "C-Terminal_CYS_(-SH)_HN",      "10" },
                     { "C-Terminal_CYS_(-SH)_N",       "7" },
                     { "C-Terminal_CYS_(-SH)_OXT",     "193" },
                     { "C-Terminal_CYS_(-SS)_C",       "192" },
                     { "C-Terminal_CYS_(-SS)_CA",      "44" },
                     { "C-Terminal_CYS_(-SS)_HA",      "12" },
                     { "C-Terminal_CYS_(-SS)_HN",      "10" },
                     { "C-Terminal_CYS_(-SS)_N",       "7" },
                     { "C-Terminal_CYS_(-SS)_OXT",     "193" },
                     { "C-Terminal_GLN_C",             "192" },
                     { "C-Terminal_GLN_CA",            "8" },
                     { "C-Terminal_GLN_HA",            "12" },
                     { "C-Terminal_GLN_HN",            "10" },
                     { "C-Terminal_GLN_N",             "7" },
                     { "C-Terminal_GLN_OXT",           "193" },
                     { "C-Terminal_GLU_C",             "192" },
                     { "C-Terminal_GLU_CA",            "8" },
                     { "C-Terminal_GLU_HA",            "12" },
                     { "C-Terminal_GLU_HN",            "10" },
                     { "C-Terminal_GLU_N",             "7" },
                     { "C-Terminal_GLU_OXT",           "193" },
                     { "C-Terminal_GLY_C",             "192" },
                     { "C-Terminal_GLY_CA",            "2" },
                     { "C-Terminal_GLY_HA",            "6" },
                     { "C-Terminal_GLY_HN",            "4" },
                     { "C-Terminal_GLY_N",             "1" },
                     { "C-Terminal_GLY_OXT",           "193" },
                     { "C-Terminal_HIS_(+)_C",         "192" },
                     { "C-Terminal_HIS_(+)_CA",        "8" },
                     { "C-Terminal_HIS_(+)_HA",        "12" },
                     { "C-Terminal_HIS_(+)_HN",        "10" },
                     { "C-Terminal_HIS_(+)_N",         "7" },
                     { "C-Terminal_HIS_(+)_OXT",       "193" },
                     { "C-Terminal_HIS_(HD)_C",        "192" },
                     { "C-Terminal_HIS_(HD)_CA",       "8" },
                     { "C-Terminal_HIS_(HD)_HA",       "12" },
                     { "C-Terminal_HIS_(HD)_HN",       "10" },
                     { "C-Terminal_HIS_(HD)_N",        "7" },
                     { "C-Terminal_HIS_(HD)_OXT",      "193" },
                     { "C-Terminal_HIS_(HE)_C",        "192" },
                     { "C-Terminal_HIS_(HE)_CA",       "8" },
                     { "C-Terminal_HIS_(HE)_HA",       "12" },
                     { "C-Terminal_HIS_(HE)_HN",       "10" },
                     { "C-Terminal_HIS_(HE)_N",        "7" },
                     { "C-Terminal_HIS_(HE)_OXT",      "193" },
                     { "C-Terminal_ILE_C",             "192" },
                     { "C-Terminal_ILE_CA",            "8" },
                     { "C-Terminal_ILE_HA",            "12" },
                     { "C-Terminal_ILE_HN",            "10" },
                     { "C-Terminal_ILE_N",             "7" },
                     { "C-Terminal_ILE_OXT",           "193" },
                     { "C-Terminal_LEU_C",             "192" },
                     { "C-Terminal_LEU_CA",            "8" },
                     { "C-Terminal_LEU_HA",            "12" },
                     { "C-Terminal_LEU_HN",            "10" },
                     { "C-Terminal_LEU_N",             "7" },
                     { "C-Terminal_LEU_OXT",           "193" },
                     { "C-Terminal_LYS_C",             "192" },
                     { "C-Terminal_LYS_CA",            "8" },
                     { "C-Terminal_LYS_HA",            "12" },
                     { "C-Terminal_LYS_HN",            "10" },
                     { "C-Terminal_LYS_N",             "7" },
                     { "C-Terminal_LYS_OXT",           "193" },
                     { "C-Terminal_MET_C",             "192" },
                     { "C-Terminal_MET_CA",            "8" },
                     { "C-Terminal_MET_HA",            "12" },
                     { "C-Terminal_MET_HN",            "10" },
                     { "C-Terminal_MET_N",             "7" },
                     { "C-Terminal_MET_OXT",           "193" },
                     { "C-Terminal_ORN_C",             "999" },
                     { "C-Terminal_ORN_CA",            "999" },
                     { "C-Terminal_ORN_HA",            "999" },
                     { "C-Terminal_ORN_HN",            "999" },
                     { "C-Terminal_ORN_N",             "999" },
                     { "C-Terminal_ORN_OXT",           "999" },
                     { "C-Terminal_PHE_C",             "192" },
                     { "C-Terminal_PHE_CA",            "8" },
                     { "C-Terminal_PHE_HA",            "12" },
                     { "C-Terminal_PHE_HN",            "10" },
                     { "C-Terminal_PHE_N",             "7" },
                     { "C-Terminal_PHE_OXT",           "193" },
                     { "C-Terminal_PRO_C",             "192" },
                     { "C-Terminal_PRO_CA",            "51" },
                     { "C-Terminal_PRO_HA",            "54" },
                     { "C-Terminal_PRO_N",             "50" },
                     { "C-Terminal_PRO_OXT",           "193" },
                     { "C-Terminal_SER_C",             "192" },
                     { "C-Terminal_SER_CA",            "33" },
                     { "C-Terminal_SER_HA",            "12" },
                     { "C-Terminal_SER_HN",            "10" },
                     { "C-Terminal_SER_N",             "7" },
                     { "C-Terminal_SER_OXT",           "193" },
                     { "C-Terminal_THR_C",             "192" },
                     { "C-Terminal_THR_CA",            "33" },
                     { "C-Terminal_THR_HA",            "12" },
                     { "C-Terminal_THR_HN",            "10" },
                     { "C-Terminal_THR_N",             "7" },
                     { "C-Terminal_THR_OXT",           "193" },
                     { "C-Terminal_TRP_C",             "192" },
                     { "C-Terminal_TRP_CA",            "8" },
                     { "C-Terminal_TRP_HA",            "12" },
                     { "C-Terminal_TRP_HN",            "10" },
                     { "C-Terminal_TRP_N",             "7" },
                     { "C-Terminal_TRP_OXT",           "193" },
                     { "C-Terminal_TYR_C",             "192" },
                     { "C-Terminal_TYR_CA",            "8" },
                     { "C-Terminal_TYR_HA",            "12" },
                     { "C-Terminal_TYR_HN",            "10" },
                     { "C-Terminal_TYR_N",             "7" },
                     { "C-Terminal_TYR_OXT",           "193" },
                     { "C-Terminal_VAL_C",             "192" },
                     { "C-Terminal_VAL_CA",            "8" },
                     { "C-Terminal_VAL_HA",            "12" },
                     { "C-Terminal_VAL_HN",            "10" },
                     { "C-Terminal_VAL_N",             "7" },
                     { "C-Terminal_VAL_OXT",           "193" },
                     { "Cysteine_(-SH)_C",             "9" },
                     { "Cysteine_(-SH)_CA",            "44" },
                     { "Cysteine_(-SH)_CB",            "45" },
                     { "Cysteine_(-SH)_HA",            "12" },
                     { "Cysteine_(-SH)_HB",            "46" },
                     { "Cysteine_(-SH)_HG",            "48" },
                     { "Cysteine_(-SH)_HN",            "10" },
                     { "Cysteine_(-SH)_N",             "7" },
                     { "Cysteine_(-SH)_O",             "11" },
                     { "Cysteine_(-SH)_SG",            "47" },
                     { "Cystine_(-SS-)_C",             "9" },
                     { "Cystine_(-SS-)_CA",            "44" },
                     { "Cystine_(-SS-)_CB",            "45" },
                     { "Cystine_(-SS-)_HA",            "12" },
                     { "Cystine_(-SS-)_HB",            "46" },
                     { "Cystine_(-SS-)_HN",            "10" },
                     { "Cystine_(-SS-)_N",             "7" },
                     { "Cystine_(-SS-)_O",             "11" },
                     { "Cystine_(-SS-)_SG",            "49" },
                     { "Formyl_N-Terminus_C",          "999" },
                     { "Formyl_N-Terminus_H",          "999" },
                     { "Formyl_N-Terminus_O",          "999" },
                     { "Glutamic_Acid_C",              "9" },
                     { "Glutamic_Acid_CA",             "8" },
                     { "Glutamic_Acid_CB",             "144" },
                     { "Glutamic_Acid_CD",             "148" },
                     { "Glutamic_Acid_CG",             "146" },
                     { "Glutamic_Acid_HA",             "12" },
                     { "Glutamic_Acid_HB",             "145" },
                     { "Glutamic_Acid_HG",             "147" },
                     { "Glutamic_Acid_HN",             "10" },
                     { "Glutamic_Acid_N",              "7" },
                     { "Glutamic_Acid_O",              "11" },
                     { "Glutamic_Acid_OE",             "149" },
                     { "Glutamine_C",                  "9" },
                     { "Glutamine_CA",                 "8" },
                     { "Glutamine_CB",                "150" },
                     { "Glutamine_CD",                "140" },
                     { "Glutamine_CG",             "138" },
                     { "Glutamine_HA",             "12" },
                     { "Glutamine_HB",             "151" },
                     { "Glutamine_HE2",             "143" },
                     { "Glutamine_HG",             "139" },
                     { "Glutamine_HN",             "10" },
                     { "Glutamine_N",             "7" },
                     { "Glutamine_NE2",             "142" },
                     { "Glutamine_O",             "11" },
                     { "Glutamine_OE1",             "141" },
                     { "Glycine_C",             "3" },
                     { "Glycine_CA",             "2" },
                     { "Glycine_HA",             "6" },
                     { "Glycine_HN",             "4" },
                     { "Glycine_N",             "1" },
                     { "Glycine_O",             "5" },
                     { "Histidine_(+)_C",             "9" },
                     { "Histidine_(+)_CA",             "8" },
                     { "Histidine_(+)_CB",             "97" },
                     { "Histidine_(+)_CD2",             "102" },
                     { "Histidine_(+)_CE1",             "104" },
                     { "Histidine_(+)_CG",             "99" },
                     { "Histidine_(+)_HA",             "12" },
                     { "Histidine_(+)_HB",             "98" },
                     { "Histidine_(+)_HD1",             "101" },
                     { "Histidine_(+)_HD2",             "103" },
                     { "Histidine_(+)_HE1",             "105" },
                     { "Histidine_(+)_HE2",             "107" },
                     { "Histidine_(+)_HN",             "10" },
                     { "Histidine_(+)_N",             "7" },
                     { "Histidine_(+)_ND1",             "100" },
                     { "Histidine_(+)_NE2",             "106" },
                     { "Histidine_(+)_O",             "11" },
                     { "Histidine_(HD)_C",             "9" },
                     { "Histidine_(HD)_CA",             "8" },
                     { "Histidine_(HD)_CB",             "108" },
                     { "Histidine_(HD)_CD2",             "113" },
                     { "Histidine_(HD)_CE1",             "115" },
                     { "Histidine_(HD)_CG",             "110" },
                     { "Histidine_(HD)_HA",             "12" },
                     { "Histidine_(HD)_HB",             "109" },
                     { "Histidine_(HD)_HD1",             "112" },
                     { "Histidine_(HD)_HD2",             "114" },
                     { "Histidine_(HD)_HE1",             "116" },
                     { "Histidine_(HD)_HN",             "10" },
                     { "Histidine_(HD)_N",             "7" },
                     { "Histidine_(HD)_ND1",             "111" },
                     { "Histidine_(HD)_NE2",             "117" },
                     { "Histidine_(HD)_O",             "11" },
                     { "Histidine_(HE)_C",             "9" },
                     { "Histidine_(HE)_CA",             "8" },
                     { "Histidine_(HE)_CB",             "118" },
                     { "Histidine_(HE)_CD2",             "122" },
                     { "Histidine_(HE)_CE1",             "124" },
                     { "Histidine_(HE)_CG",             "120" },
                     { "Histidine_(HE)_HA",             "12" },
                     { "Histidine_(HE)_HB",             "119" },
                     { "Histidine_(HE)_HD2",             "123" },
                     { "Histidine_(HE)_HE1",             "125" },
                     { "Histidine_(HE)_HE2",             "127" },
                     { "Histidine_(HE)_HN",             "10" },
                     { "Histidine_(HE)_N",             "7" },
                     { "Histidine_(HE)_ND1",             "121" },
                     { "Histidine_(HE)_NE2",             "126" },
                     { "Histidine_(HE)_O",             "11" },
                     { "Isoleucine_C",             "9" },
                     { "Isoleucine_CA",             "8" },
                     { "Isoleucine_CB",             "25" },
                     { "Isoleucine_CD",             "31" },
                     { "Isoleucine_CG1",             "29" },
                     { "Isoleucine_CG2",             "27" },
                     { "Isoleucine_HA",             "12" },
                     { "Isoleucine_HB",             "26" },
                     { "Isoleucine_HD",             "32" },
                     { "Isoleucine_HG1",             "30" },
                     { "Isoleucine_HG2",             "28" },
                     { "Isoleucine_HN",             "10" },
                     { "Isoleucine_N",             "7" },
                     { "Isoleucine_O",             "11" },
                     { "Leucine_C",             "9" },
                     { "Leucine_CA",             "8" },
                     { "Leucine_CB",             "19" },
                     { "Leucine_CD1",             "23" },
                     { "Leucine_CD2",             "23" },
                     { "Leucine_CG",             "21" },
                     { "Leucine_HA",             "12" },
                     { "Leucine_HB",             "20" },
                     { "Leucine_HD1",             "24" },
                     { "Leucine_HD2",             "24" },
                     { "Leucine_HG",             "22" },
                     { "Leucine_HN",             "10" },
                     { "Leucine_N",             "7" },
                     { "Leucine_O",             "11" },
                     { "Lysine_C",             "9" },
                     { "Lysine_CA",             "8" },
                     { "Lysine_CB",             "159" },
                     { "Lysine_CD",             "163" },
                     { "Lysine_CE",             "165" },
                     { "Lysine_CG",             "161" },
                     { "Lysine_HA",             "12" },
                     { "Lysine_HB",             "160" },
                     { "Lysine_HD",             "164" },
                     { "Lysine_HE",             "166" },
                     { "Lysine_HG",             "162" },
                     { "Lysine_HN",             "10" },
                     { "Lysine_HZ",             "168" },
                     { "Lysine_N",             "7" },
                     { "Lysine_NZ",             "167" },
                     { "Lysine_O",             "11" },
                     { "Methionine_C",             "9" },
                     { "Methionine_CA",             "8" },
                     { "Methionine_CB",             "152" },
                     { "Methionine_CE",             "157" },
                     { "Methionine_CG",             "154" },
                     { "Methionine_HA",             "12" },
                     { "Methionine_HB",             "153" },
                     { "Methionine_HE",             "158" },
                     { "Methionine_HG",             "155" },
                     { "Methionine_HN",             "10" },
                     { "Methionine_N",             "7" },
                     { "Methionine_O",             "11" },
                     { "Methionine_SD",             "156" },
                     { "MethylAlanine_(AIB)_C",             "999" },
                     { "MethylAlanine_(AIB)_CA",             "999" },
                     { "MethylAlanine_(AIB)_CB",             "999" },
                     { "MethylAlanine_(AIB)_HB",             "999" },
                     { "MethylAlanine_(AIB)_HN",             "999" },
                     { "MethylAlanine_(AIB)_N",             "999" },
                     { "MethylAlanine_(AIB)_O",             "999" },
                     { "N-MeAmide_C-Terminus_CH3",             "188" },
                     { "N-MeAmide_C-Terminus_H",             "189" },
                     { "N-MeAmide_C-Terminus_HN",             "187" },
                     { "N-MeAmide_C-Terminus_N",             "186" },
                     { "N-Terminal_AIB_C",             "999" },
                     { "N-Terminal_AIB_CA",             "999" },
                     { "N-Terminal_AIB_HN",             "999" },
                     { "N-Terminal_AIB_N",             "999" },
                     { "N-Terminal_AIB_O",             "999" },
                     { "N-Terminal_ALA_C",             "9" },
                     { "N-Terminal_ALA_CA",             "8" },
                     { "N-Terminal_ALA_HA",             "12" },
                     { "N-Terminal_ALA_HN",             "191" },
                     { "N-Terminal_ALA_N",             "190" },
                     { "N-Terminal_ALA_O",             "11" },
                     { "N-Terminal_ARG_C",             "9" },
                     { "N-Terminal_ARG_CA",             "8" },
                     { "N-Terminal_ARG_HA",             "12" },
                     { "N-Terminal_ARG_HN",             "191" },
                     { "N-Terminal_ARG_N",             "190" },
                     { "N-Terminal_ARG_O",             "11" },
                     { "N-Terminal_ASN_C",             "9" },
                     { "N-Terminal_ASN_CA",             "8" },
                     { "N-Terminal_ASN_HA",             "12" },
                     { "N-Terminal_ASN_HN",             "191" },
                     { "N-Terminal_ASN_N",             "190" },
                     { "N-Terminal_ASN_O",             "11" },
                     { "N-Terminal_ASP_C",             "9" },
                     { "N-Terminal_ASP_CA",             "8" },
                     { "N-Terminal_ASP_HA",             "12" },
                     { "N-Terminal_ASP_HN",             "191" },
                     { "N-Terminal_ASP_N",             "190" },
                     { "N-Terminal_ASP_O",             "11" },
                     { "N-Terminal_CYS_(-SH)_C",             "9" },
                     { "N-Terminal_CYS_(-SH)_CA",             "44" },
                     { "N-Terminal_CYS_(-SH)_HA",             "12" },
                     { "N-Terminal_CYS_(-SH)_HN",             "191" },
                     { "N-Terminal_CYS_(-SH)_N",             "190" },
                     { "N-Terminal_CYS_(-SH)_O",             "11" },
                     { "N-Terminal_CYS_(-SS)_C",             "9" },
                     { "N-Terminal_CYS_(-SS)_CA",             "44" },
                     { "N-Terminal_CYS_(-SS)_HA",             "12" },
                     { "N-Terminal_CYS_(-SS)_HN",             "191" },
                     { "N-Terminal_CYS_(-SS)_N",             "190" },
                     { "N-Terminal_CYS_(-SS)_O",             "11" },
                     { "N-Terminal_GLN_C",             "9" },
                     { "N-Terminal_GLN_CA",             "8" },
                     { "N-Terminal_GLN_HA",             "12" },
                     { "N-Terminal_GLN_HN",             "191" },
                     { "N-Terminal_GLN_N",             "190" },
                     { "N-Terminal_GLN_O",             "11" },
                     { "N-Terminal_GLU_C",             "9" },
                     { "N-Terminal_GLU_CA",             "8" },
                     { "N-Terminal_GLU_HA",             "12" },
                     { "N-Terminal_GLU_HN",             "191" },
                     { "N-Terminal_GLU_N",             "190" },
                     { "N-Terminal_GLU_O",             "11" },
                     { "N-Terminal_GLY_C",             "3" },
                     { "N-Terminal_GLY_CA",             "2" },
                     { "N-Terminal_GLY_HA",             "6" },
                     { "N-Terminal_GLY_HN",             "191" },
                     { "N-Terminal_GLY_N",             "190" },
                     { "N-Terminal_GLY_O",             "5" },
                     { "N-Terminal_HIS_(+)_C",             "9" },
                     { "N-Terminal_HIS_(+)_CA",             "8" },
                     { "N-Terminal_HIS_(+)_HA",             "12" },
                     { "N-Terminal_HIS_(+)_HN",             "191" },
                     { "N-Terminal_HIS_(+)_N",             "190" },
                     { "N-Terminal_HIS_(+)_O",             "11" },
                     { "N-Terminal_HIS_(HD)_C",             "9" },
                     { "N-Terminal_HIS_(HD)_CA",             "8" },
                     { "N-Terminal_HIS_(HD)_HA",             "12" },
                     { "N-Terminal_HIS_(HD)_HN",             "191" },
                     { "N-Terminal_HIS_(HD)_N",             "190" },
                     { "N-Terminal_HIS_(HD)_O",             "11" },
                     { "N-Terminal_HIS_(HE)_C",             "9" },
                     { "N-Terminal_HIS_(HE)_CA",             "8" },
                     { "N-Terminal_HIS_(HE)_HA",             "12" },
                     { "N-Terminal_HIS_(HE)_HN",             "191" },
                     { "N-Terminal_HIS_(HE)_N",             "190" },
                     { "N-Terminal_HIS_(HE)_O",             "11" },
                     { "N-Terminal_ILE_C",             "9" },
                     { "N-Terminal_ILE_CA",             "8" },
                     { "N-Terminal_ILE_HA",             "12" },
                     { "N-Terminal_ILE_HN",             "191" },
                     { "N-Terminal_ILE_N",             "190" },
                     { "N-Terminal_ILE_O",             "11" },
                     { "N-Terminal_LEU_C",             "9" },
                     { "N-Terminal_LEU_CA",             "8" },
                     { "N-Terminal_LEU_HA",             "12" },
                     { "N-Terminal_LEU_HN",             "191" },
                     { "N-Terminal_LEU_N",             "190" },
                     { "N-Terminal_LEU_O",             "11" },
                     { "N-Terminal_LYS_C",             "9" },
                     { "N-Terminal_LYS_CA",             "8" },
                     { "N-Terminal_LYS_HA",             "12" },
                     { "N-Terminal_LYS_HN",             "191" },
                     { "N-Terminal_LYS_N",             "190" },
                     { "N-Terminal_LYS_O",             "11" },
                     { "N-Terminal_MET_C",             "9" },
                     { "N-Terminal_MET_CA",             "8" },
                     { "N-Terminal_MET_HA",             "12" },
                     { "N-Terminal_MET_HN",             "191" },
                     { "N-Terminal_MET_N",             "190" },
                     { "N-Terminal_MET_O",             "11" },
                     { "N-Terminal_ORN_C",             "999" },
                     { "N-Terminal_ORN_CA",             "999" },
                     { "N-Terminal_ORN_HA",             "999" },
                     { "N-Terminal_ORN_HN",             "999" },
                     { "N-Terminal_ORN_N",             "999" },
                     { "N-Terminal_ORN_O",             "999" },
                     { "N-Terminal_PHE_C",             "9" },
                     { "N-Terminal_PHE_CA",             "8" },
                     { "N-Terminal_PHE_HA",             "12" },
                     { "N-Terminal_PHE_HN",             "191" },
                     { "N-Terminal_PHE_N",             "190" },
                     { "N-Terminal_PHE_O",             "11" },
                     { "N-Terminal_PRO_C",             "197" },
                     { "N-Terminal_PRO_CA",             "196" },
                     { "N-Terminal_PRO_CD",             "200" },
                     { "N-Terminal_PRO_HA",             "199" },
                     { "N-Terminal_PRO_HD",             "201" },
                     { "N-Terminal_PRO_HN",             "195" },
                     { "N-Terminal_PRO_N",             "194" },
                     { "N-Terminal_PRO_O",             "198" },
                     { "N-Terminal_SER_C",             "9" },
                     { "N-Terminal_SER_CA",             "33" },
                     { "N-Terminal_SER_HA",             "12" },
                     { "N-Terminal_SER_HN",             "191" },
                     { "N-Terminal_SER_N",             "190" },
                     { "N-Terminal_SER_O",             "11" },
                     { "N-Terminal_THR_C",             "9" },
                     { "N-Terminal_THR_CA",             "33" },
                     { "N-Terminal_THR_HA",             "12" },
                     { "N-Terminal_THR_HN",             "191" },
                     { "N-Terminal_THR_N",             "190" },
                     { "N-Terminal_THR_O",             "11" },
                     { "N-Terminal_TRP_C",             "9" },
                     { "N-Terminal_TRP_CA",             "8" },
                     { "N-Terminal_TRP_HA",             "12" },
                     { "N-Terminal_TRP_HN",             "191" },
                     { "N-Terminal_TRP_N",             "190" },
                     { "N-Terminal_TRP_O",             "11" },
                     { "N-Terminal_TYR_C",             "9" },
                     { "N-Terminal_TYR_CA",             "8" },
                     { "N-Terminal_TYR_HA",             "12" },
                     { "N-Terminal_TYR_HN",             "191" },
                     { "N-Terminal_TYR_N",             "190" },
                     { "N-Terminal_TYR_O",             "11" },
                     { "N-Terminal_VAL_C",             "9" },
                     { "N-Terminal_VAL_CA",             "8" },
                     { "N-Terminal_VAL_HA",             "12" },
                     { "N-Terminal_VAL_HN",             "191" },
                     { "N-Terminal_VAL_N",             "190" },
                     { "N-Terminal_VAL_O",             "11" },
                     { "Ornithine_C",             "999" },
                     { "Ornithine_CA",             "999" },
                     { "Ornithine_CB",             "999" },
                     { "Ornithine_CD",             "999" },
                     { "Ornithine_CG",             "999" },
                     { "Ornithine_HA",             "999" },
                     { "Ornithine_HB",             "999" },
                     { "Ornithine_HD",             "999" },
                     { "Ornithine_HE",             "999" },
                     { "Ornithine_HG",             "999" },
                     { "Ornithine_HN",             "999" },
                     { "Ornithine_N",             "999" },
                     { "Ornithine_NE",             "999" },
                     { "Ornithine_O",             "999" },
                     { "Phenylalanine_C",             "9" },
                     { "Phenylalanine_CA",             "8" },
                     { "Phenylalanine_CB",             "61" },
                     { "Phenylalanine_CD",             "64" },
                     { "Phenylalanine_CE",             "66" },
                     { "Phenylalanine_CG",             "63" },
                     { "Phenylalanine_CZ",             "68" },
                     { "Phenylalanine_HA",             "12" },
                     { "Phenylalanine_HB",             "62" },
                     { "Phenylalanine_HD",             "65" },
                     { "Phenylalanine_HE",             "67" },
                     { "Phenylalanine_HN",             "10" },
                     { "Phenylalanine_HZ",             "69" },
                     { "Phenylalanine_N",             "7" },
                     { "Phenylalanine_O",             "11" },
                     { "Proline_C",             "52" },
                     { "Proline_CA",             "51" },
                     { "Proline_CB",             "55" },
                     { "Proline_CD",             "59" },
                     { "Proline_CG",             "57" },
                     { "Proline_HA",             "54" },
                     { "Proline_HB",             "56" },
                     { "Proline_HD",             "60" },
                     { "Proline_HG",             "58" },
                     { "Proline_N",             "50" },
                     { "Proline_O",             "53" },
                     { "Pyroglutamic_Acid_C",             "999" },
                     { "Pyroglutamic_Acid_CA",             "999" },
                     { "Pyroglutamic_Acid_CB",             "999" },
                     { "Pyroglutamic_Acid_CD",             "999" },
                     { "Pyroglutamic_Acid_CG",             "999" },
                     { "Pyroglutamic_Acid_HA",             "999" },
                     { "Pyroglutamic_Acid_HB",             "999" },
                     { "Pyroglutamic_Acid_HG",             "999" },
                     { "Pyroglutamic_Acid_HN",             "999" },
                     { "Pyroglutamic_Acid_N",             "999" },
                     { "Pyroglutamic_Acid_O",             "999" },
                     { "Pyroglutamic_Acid_OE",             "999" },
                     { "Serine_C",             "9" },
                     { "Serine_CA",             "33" },
                     { "Serine_CB",             "34" },
                     { "Serine_HA",             "12" },
                     { "Serine_HB",             "35" },
                     { "Serine_HG",             "37" },
                     { "Serine_HN",             "10" },
                     { "Serine_N",             "7" },
                     { "Serine_O",             "11" },
                     { "Serine_OG",             "36" },
                     { "Threonine_C",             "9" },
                     { "Threonine_CA",             "33" },
                     { "Threonine_CB",             "38" },
                     { "Threonine_CG2",             "40" },
                     { "Threonine_HA",             "12" },
                     { "Threonine_HB",             "39" },
                     { "Threonine_HG1",             "43" },
                     { "Threonine_HG2",             "41" },
                     { "Threonine_HN",             "10" },
                     { "Threonine_N",             "7" },
                     { "Threonine_O",             "11" },
                     { "Threonine_OG1",             "42" },
                     { "Tryptophan_C",             "9" },
                     { "Tryptophan_CA",             "8" },
                     { "Tryptophan_CB",             "80" },
                     { "Tryptophan_CD1",             "83" },
                     { "Tryptophan_CD2",             "85" },
                     { "Tryptophan_CE2",             "88" },
                     { "Tryptophan_CE3",             "89" },
                     { "Tryptophan_CG",             "82" },
                     { "Tryptophan_CH2",             "95" },
                     { "Tryptophan_CZ2",             "91" },
                     { "Tryptophan_CZ3",             "93" },
                     { "Tryptophan_HA",             "12" },
                     { "Tryptophan_HB",             "81" },
                     { "Tryptophan_HD1",             "84" },
                     { "Tryptophan_HE1",             "87" },
                     { "Tryptophan_HE3",             "90" },
                     { "Tryptophan_HH2",             "96" },
                     { "Tryptophan_HN",             "10" },
                     { "Tryptophan_HZ2",             "92" },
                     { "Tryptophan_HZ3",             "94" },
                     { "Tryptophan_N",             "7" },
                     { "Tryptophan_NE1",             "86" },
                     { "Tryptophan_O",             "11" },
                     { "Tyrosine_C",             "9" },
                     { "Tyrosine_CA",             "8" },
                     { "Tyrosine_CB",             "70" },
                     { "Tyrosine_CD",             "73" },
                     { "Tyrosine_CE",             "75" },
                     { "Tyrosine_CG",             "72" },
                     { "Tyrosine_CZ",             "77" },
                     { "Tyrosine_HA",             "12" },
                     { "Tyrosine_HB",             "71" },
                     { "Tyrosine_HD",             "74" },
                     { "Tyrosine_HE",             "76" },
                     { "Tyrosine_HH",             "79" },
                     { "Tyrosine_HN",             "10" },
                     { "Tyrosine_N",             "7" },
                     { "Tyrosine_O",             "11" },
                     { "Tyrosine_OH",             "78" },
                     { "Valine_C",             "9" },
                     { "Valine_CA",             "8" },
                     { "Valine_CB",             "15" },
                     { "Valine_CG1",             "17" },
                     { "Valine_CG2",             "17" },
                     { "Valine_HA",             "12" },
                     { "Valine_HB",             "16" },
                     { "Valine_HG1",             "18" },
                     { "Valine_HG2",             "18" },
                     { "Valine_HN",             "10" },
                     { "Valine_N",             "7" },
                     { "Valine_O",             "11" },
                     {  NULL,  NULL     }
                               };


   static StringIntMap* residueAtomNameBiotypeMap = NULL;

   // ---------------------------------------------------------------------------------------

   // residueNameMap is singleton

   if( residueAtomNameBiotypeMap == NULL ){

      residueAtomNameBiotypeMap = new StringIntMap();

      // unsigned int initialMapSize = sizeof( nameMap )/(2*sizeof( char* ));

      int index = 0;
      while( nameMap[index][0] != NULL ){
         (*residueAtomNameBiotypeMap)[nameMap[index][0]] = atoi( nameMap[index][1] );
         index++;
      }

   }

   return residueAtomNameBiotypeMap;
}

/**---------------------------------------------------------------------------------------

   Get solvent radii from parameter file (Simbios) 
      
   @param numberOfAtoms       number of atoms
   @param parameterFileName   parameter file name
   @param top                 Gromacs topology data struct
   @param radii               array to store Macromodel radii for each atom
   @param scaleFactor         scale factor
      
   @return SimTKOpenMMCommon::DefaultReturn unless paramter file not opened
           in which case return SimTKOpenMMCommon::ErrorReturn
      
   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getMacroModelAtomicRadii( int numberOfAtoms,
                                                           const std::string parameterFileName,
                                                           const t_topology* top, RealOpenMM* radii,
                                                           RealOpenMM scaleFactor ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getMacroModelAtomicRadii";

   // ---------------------------------------------------------------------------------------

 
   RealOpenMMVector radiiV;
   int status = SimTKOpenMMGromacsUtilities::getMacroModelAtomicRadii( numberOfAtoms, parameterFileName,
                                                                       top, radiiV, scaleFactor );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      radii[ii] = radiiV[ii];
   }

   return status;

}

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

int SimTKOpenMMGromacsUtilities::getMacroModelAtomicRadii( int numberOfAtoms, const std::string parameterFileName,
                                                           const t_topology* top, RealOpenMMVector& radii,
                                                           RealOpenMM scaleFactor ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getMacroModelAtomicRadii";

   // ---------------------------------------------------------------------------------------

   // read in file contents -- return error flag if file not opened

   StringVector fileContents;
   int returnStatus = SimTKOpenMMUtilities::readFileIntoStringVector( parameterFileName, fileContents ); 
   std::stringstream message;

   if( returnStatus != SimTKOpenMMCommon::DefaultReturn ){
      message << " agb parameter file not opened=<" << parameterFileName <<">";
      return SimTKOpenMMCommon::ErrorReturn;
   } else {
      message << methodName << " read " << fileContents.size() << " lines from agb parameter file=<" << parameterFileName << ">.";
   }
   SimTKOpenMMLog::printMessage( message );

   // parse file: parameterMap->[atomType] = radius

   StringRealOpenMMMap parameterMap;
   std::string delimiter = " ";
   StringVector tokenVector;

   std::stringstream warningMessage;
   warningMessage << methodName << " file=<" << parameterFileName << "> ";
   int printWarning = 0;

   for( StringVectorCI ii = fileContents.begin(); ii != fileContents.end(); ii++ ){

      std::string line = *ii;
      SimTKOpenMMUtilities::tokenizeString( line, tokenVector, delimiter, 1 );
      if( tokenVector[0] != "#" && tokenVector[0] != "@" && tokenVector.size() >= 2 ){
         parameterMap[tokenVector[0]] = (RealOpenMM) atof( tokenVector[1].c_str() );
      } else {
         printWarning = 1;
         warningMessage << "\n   agb parameter file line=<" << line << "> is being skipped.";
      }
   }
   if( printWarning ){
      SimTKOpenMMLog::printMessage( warningMessage );
   }

   char*** atomNames = top->atoms.atomname;

   // get radii

   radii.resize( numberOfAtoms );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      char* atomTypeC = *(top->atoms.atomtype[ii]);   
      std::string atomType  = atomTypeC;
      StringRealOpenMMMapCI entry = parameterMap.find( atomType ); 
      if( entry != parameterMap.end() ){

         radii[ii]   = scaleFactor*((*entry).second);

         //std::stringstream message;
         //message << "\ntype found for atom=<" << (*(atomNames[ii])) << "> type=<" << atomType << "> r=" << radii[ii];
         //SimTKOpenMMLog::printMessage( message );

      } else {
         radii[ii] = (RealOpenMM) 0.000001;
         std::stringstream message;
         message << methodName;
         message << " no type found for atom=<" << (*(atomNames[ii])) << "> type=<" << atomType << ">";
         SimTKOpenMMLog::printMessage( message );
      }
   }
   
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get string containing atom types
      
   @param top                 Gromacs topology data struct
      
   @return string
      
   --------------------------------------------------------------------------------------- */

std::string SimTKOpenMMGromacsUtilities::getAtomTypesString( const t_topology* top ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getAtomTypesString";

   // ---------------------------------------------------------------------------------------

   char*** atomNames = top->atoms.atomname;

   std::stringstream message;

   for( int ii = 0; ii < top->atoms.nr; ii++ ){
      char* atomTypeC = *(top->atoms.atomtype[ii]);   
      message << "\n" << (ii+1) << " " << (*(atomNames[ii])) << " " << atomTypeC;
   }
   
   return message.str();
}

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

RealOpenMM** SimTKOpenMMGromacsUtilities::copyRvecArrayToRealOpenMMArray( int numberOfEntries,
                                                                          rvec* gromacsArray,
                                                                          RealOpenMM** realArray,
                                                                          RealOpenMM scaleFactor ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::copyRvecArrayToRealOpenMMArray";

   // ---------------------------------------------------------------------------------------

   if( realArray == NULL ){
      realArray = SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( numberOfEntries, 3, NULL, false, 0.0 ); 
   }

   for( int ii = 0; ii < numberOfEntries; ii++ ){
      realArray[ii][0] = scaleFactor*gromacsArray[ii][0];
      realArray[ii][1] = scaleFactor*gromacsArray[ii][1];
      realArray[ii][2] = scaleFactor*gromacsArray[ii][2];
   }

   return realArray;
}

/**---------------------------------------------------------------------------------------

   Get vectors of atom names, residue indices, and residue names from 
   Gromacs data structs (Simbios)

   @param top            Gromacs t_topolgy struct.
   @param residueNames   output residue names
   @param residueIndices output residue indices
   @param atomNames      output atom names

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMGromacsUtilities::getAtomResidueNames( const t_topology* top,
                                                      StringVector& residueNamesVector,
                                                      IntVector& residueIndicesVector,
                                                      StringVector& atomNamesVector ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMGromacsUtilities::getAtomResidueNames";

   // ---------------------------------------------------------------------------------------

   int numberOfAtoms     = top->atoms.nr;

   residueNamesVector.resize( numberOfAtoms );
   residueIndicesVector.resize( numberOfAtoms );
   atomNamesVector.resize( numberOfAtoms );

   // Gromacs arrays of atom and residue names

   char*** atomNames     = top->atoms.atomname;
   char*** residueNames  = top->atoms.resname;

   // for each atom, set names, indices, ...

   int residueI          = 0;
   char* firstAtomName   = NULL;

   std::string residueName;
   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // check if hitting new residue
      // save first atom in list to allow for check
      // when another atom is hit

      if( !firstAtomName ){
         firstAtomName  = *(atomNames[atomI]);
         residueName    = std::string( *(residueNames[residueI]) );
      } else if( !strcmp( *(atomNames[atomI]), firstAtomName ) ){
         residueI++;
         residueName    = std::string( *(residueNames[residueI]) );
      }

      atomNamesVector.push_back( *(atomNames[atomI]) );
      residueNamesVector.push_back( residueName );
      residueIndicesVector.push_back( residueI + 1 );

   }

   return SimTKOpenMMCommon::DefaultReturn;

}
#endif