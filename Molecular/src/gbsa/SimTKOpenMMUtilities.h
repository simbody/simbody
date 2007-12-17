
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

#ifndef __SimTKOpenMMUtilities_H_
#define __SimTKOpenMMUtilities_H_

// class of shared, static utility methods

#include "SimTKOpenMMCommon.h" 

#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sstream>

// template is used to check if a string is integer, real, ...

template <class T>
bool checkString( T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&) ){
   std::istringstream iss(s); 
   return !(iss >> f >> t).fail();
}

/**---------------------------------------------------------------------------------------

   Class of static methods to be shared
   Most methods are standalone 'utility' methods

--------------------------------------------------------------------------------------- */

class SimTKOpenMMUtilities {

   public:

      // file flag enums

      enum FileFlags { OpenDebugFile, WriteDebugFile, CloseDebugFile };

      // dummy constructor/destructor

       SimTKOpenMMUtilities(){};
      ~SimTKOpenMMUtilities(){};

      /**---------------------------------------------------------------------------------------
      
         Find distances**2 from a given atom (Simbios)
      
         @param atomCoordinates     atom coordinates
         @param atomIndex           atom index to find distances from
         @param numberOfAtoms       number of atoms
         @param distances           array of distances squared on @return; array size must be at least
                                    numberOfAtoms
         @param log                 if set, then print error messages to log file
      
         @return distances
      
         --------------------------------------------------------------------------------------- */
      
      static int getDistanceSquaredFromSpecifiedAtom( RealOpenMM** atomCoordinates, int atomIndex,
                                                      int numberOfAtoms, RealOpenMM* distances,
                                                      FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Find distances**2 from a given point (Simbios)
      
         @param atomCoordinates     atom coordinates
         @param point               point to find distances from
         @param numberOfAtoms       number of atoms
         @param distances           array of distances squared on @return; array size must be at least
                                    numberOfAtoms
         @param log                 if set, then print error messages to log file
      
         @return distances
      
         --------------------------------------------------------------------------------------- */
      
      static int getDistanceSquaredFromSpecifiedPoint( RealOpenMM** atomCoordinates, RealOpenMM* point, 
                                                       int numberOfAtoms, RealOpenMM* distances,
                                                       FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Helper method to allocate RealOpenMM arrays (Simbios)
      
         @param bufferIndex         buffer index
         @param allocatedSz         array of allocated sizes
         @param bufferArray         array of allocated RealOpenMM arrays
         @param requestedSize       requested size
         @param dataAction          action flag: -1 = free memory \n
                                                  1 = zero memory
      
         @return SimTKOpenMMCommon::DefaultReturn 
      
         --------------------------------------------------------------------------------------- */
      
      static int allocateRealOpenMMBufferArray( int bufferIndex, int* allocatedSz,
                                                RealOpenMM** bufferArray,
                                                int requestedSize, int dataAction );
      
      /**---------------------------------------------------------------------------------------
      
         Print atom coordinates, ...
      
         @param numberAtoms         numberAtoms
         @param atomCoordinates     atomCoordinates (may be NULL)
         @param numberOf1Darrays     number of 1-d arrays (may be 0)
         @param oneDArrays          1-d arrays 
         @param idString            id string to be printed if set
         @param log                 print messages to log file
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      static int printCoordinateAnd1DArrays( int numberAtoms, RealOpenMM** atomCoordinates,
                                             int numberOf1Darrays, RealOpenMM** oneDArrays,
                                             const char* idString, FILE* log );

      /**---------------------------------------------------------------------------------------
      
         Free array of strings
      
         @param arraySz             atom index
         @param arrayOfStrings      array of strings
      
         @return SimTKOpenMMCommon::DefaultReturn

         --------------------------------------------------------------------------------------- */
      
      static int freeArrayOfStrings( int arraySz, char** arrayOfStrings );

      /**---------------------------------------------------------------------------------------
      
         Tab string in place
      
         @param string              string to tab; assume string is of at least length=tab + 1
         @param tab                 tab length
      
         --------------------------------------------------------------------------------------- */
      
      static int tabStringInPlace( char* string, int tab  );
      
      /**---------------------------------------------------------------------------------------
      
         Write debug fields (Simbios)
      
         @param numberOfFields       number of fields to print
         @param fields               fields
         @param numberOfStringFields number of string fields to print
         @param stringFields         string fields
         @param comment              comment (optinal -- ignored if NULL)
         @param debugFileName        output debug file name
         @param action               0 open file and @return w/o printing
                                     1 open file and print
                                     2 close file (no print)
         @param debugFile            debug file reference
         @param log                  if set, then print error messages to log file
      
         @return debugFile unless file is closed
      
         stringFields printed after RealOpenMM fields
      
         --------------------------------------------------------------------------------------- */
      
      static FILE* writeDebugFile( int numberOfFields, const RealOpenMM* fields,
                                   int numberOfStringFields, const  StringVector& stringFields,
                                   const char* comment, const char* debugFileName, int action,
                                   FILE* debugFile, FILE* log );
      
      /**---------------------------------------------------------------------------------------
      
         Allocate 2D RealOpenMM array (Simbios)
      
         array[i][j]
      
         @param iSize                i-dimension
         @param jSize                j-dimension
         @param array2D              array (if null on entry allocated)
         @param initialize           if true, then initialize array
         @param initialValue         intitial value
         @param idString             id string
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      static RealOpenMM** allocateTwoDRealOpenMMArray( int iSize, int jSize,
                                                       RealOpenMM** array2D, int initialize,
                                                       RealOpenMM initialValue,
                                                       const std::string& idString = std::string( "2DArray" ) );
      
      /* ---------------------------------------------------------------------------------------
      
         Free 2D RealOpenMM array (Simbios)
      
         array[i][j]
      
         @param array2D              array (if null on entry allocated)
         @param idString             id string
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int freeTwoDRealOpenMMArray( RealOpenMM** array2D,
                                          const std::string& idString = std::string( "2DArray" ) );
      
      /**---------------------------------------------------------------------------------------
      
         Initialize 2D RealOpenMM array (Simbios)
      
         array[i][j]
      
         @param iSize                i-dimension
         @param jSize                j-dimension
         @param array2D              array (if null on entry allocated)
         @param initialValue         intitial value
      
         @return array
      
         --------------------------------------------------------------------------------------- */
      
      static int initialize2DRealOpenMMArray( int iSize, int jSize,
                                              RealOpenMM** array2D, RealOpenMM initialValue );

      /**---------------------------------------------------------------------------------------
      
         Malloc memory of size bytesToAllocate and zero
      
         @param bytesToAllocate      bytes to allocate
      
         @return ptr to allocated memory; NULL if bytesToAllocate <= 0
      
         --------------------------------------------------------------------------------------- */
      
      static char* allocateAndZero( unsigned int bytesToAllocate );

      /**---------------------------------------------------------------------------------------
      
         Normalize 3-vector -- helper method
      
         @param vector              vector to normalize
      
         --------------------------------------------------------------------------------------- */
            
      static void normalizeVector3( RealOpenMM* vector );
      
      /**---------------------------------------------------------------------------------------
      
         Remove 3-vector -- helper method
      
         @param vectorToRemove      vector to remove
         @param vector              vector to from which 'vectorToRemove' is to be removed \n
                                    vector is normalized after the component is subtracted out
      
         --------------------------------------------------------------------------------------- */
      
      static void removeVector3( RealOpenMM* vectorToRemove, RealOpenMM* vector );
      
      /**---------------------------------------------------------------------------------------
      
         Compute cross product of two 3-vectors and place in 3rd vector  -- helper method
      
         @param vectorZ = vectorX x vectorY
      
         @param vectorX             x-vector
         @param vectorY             y-vector
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void crossProductVector3( RealOpenMM* vectorX, RealOpenMM* vectorY, RealOpenMM* vectorZ );

      /**---------------------------------------------------------------------------------------
      
         Compute matrix product of 3x3 matrix and 3-vector and place in 3rd vector  -- helper method
      
         @param vectorZ = matrixX . vectorY
      
         @param matrixX             matrixX
         @param vectorY             y-vector
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void matrixProductVector3( RealOpenMM* matrixX, RealOpenMM* vectorY, RealOpenMM* vectorZ );

      /**---------------------------------------------------------------------------------------
      
         Compute cross product between two 3x3 matrices
      
         @param vectorZ = matrixX . matrixY
      
         @param matrixX             matrixX
         @param matrixY             matrixY
         @param vectorZ             z-vector
      
         @return vector is vectorZ
      
         --------------------------------------------------------------------------------------- */
      
      static void matrixCrossProductMatrix3( RealOpenMM* matrixX, RealOpenMM* matrixY, RealOpenMM* vectorZ );

      /* ---------------------------------------------------------------------------------------
      
         Centralized malloc/new
      
         @param name                ptr name
         @param fileName            file name
         @param line                file line no.
         @param file line           size in bytes to be allocated
      
         @return ptr to allocated object
      
         --------------------------------------------------------------------------------------- */
               
      static void* Xmalloc( const char* name, char* fileName, int line, unsigned int size );
      
      /* ---------------------------------------------------------------------------------------
      
         Centralized free/delete
      
         @param name                ptr name
         @param fileName            file name
         @param line                file line no.
         @param ptr                 ptr to be freed
      
         --------------------------------------------------------------------------------------- */
      
      static void Xfree( const char* name, char* fileName, int line, void* ptr );
      
      /* ---------------------------------------------------------------------------------------
      
         Format array of reals
      
         @param message             input string stream
         @param realArray				array of RealOpenMMs
         @param numberOfFields      number of fields (optional - defaults to 3)
      
         @return SimTKOpenMMCommon::DefaultReturn

         --------------------------------------------------------------------------------------- */
      
      static int formatRealStringStream( std::stringstream& message, const RealOpenMM* realArray, 
                                         int numberOfFields = 3, RealOpenMM factor = (RealOpenMM) 1.0f );
      
      /**---------------------------------------------------------------------------------------
      
         Tokenize a string (static method) (Simbios)
      
         @param lineBuffer           string to tokenize
         @param tokenArray           upon return vectory of tokens
         @param delimiter            token delimter
      
         @return number of args
      
         --------------------------------------------------------------------------------------- */
      
      static int tokenizeString( char* lineBuffer, StringVector& tokenArray, const std::string delimiter = "\t\n " );
      
      /**---------------------------------------------------------------------------------------

         Tokenize a string (static method) (Simbios)


         @param line                 string to tokenize
         @param tokenVector          upon return vector of tokens
         @param delimiter            token delimter
         @param clearTokenVector     if true, clear tokenVector
      
         @return SimTKOpenMMCommon::DefaultReturn
         
         --------------------------------------------------------------------------------------- */
      
      static int tokenizeString( const std::string& line, StringVector& tokenVector,
                                 const std::string& delimiter, int clearTokenVector );
      
      /**---------------------------------------------------------------------------------------
      
         Local version of strncasecmp (missing in Windows) (static method) (Simbios)
      
         @param string1                 first string
         @param string2                 second string
         @param matchLength             match length
         
         @return SimTKOpenMMCommon::DefaultReturn
         
         --------------------------------------------------------------------------------------- */
      
      static int localStrncasecmp( const char *string1, const char *string2, int matchLength );

      /**---------------------------------------------------------------------------------------
      
         Check that string is valid integer
      
         @param stringToCheck string to check
      
         @return true if string is a valid integer
      
         --------------------------------------------------------------------------------------- */
      
      static bool isValidInteger( std::string stringToCheck );

      /**---------------------------------------------------------------------------------------
      
         Check that string is valid RealOpenMM
      
         @param stringToCheck string to check
      
         @return true if string is a valid RealOpenMM
      
         --------------------------------------------------------------------------------------- */
      
      static bool isValidRealOpenMM( std::string stringToCheck );

      /**---------------------------------------------------------------------------------------
      
         Read file into string vector (Simbios) 
      
         @param fileName       file name
         @param fileContents   string vector containing file contents upon return
                               one string per line
      
         @return SimTKOpenMMCommon::DefaultReturn unless file could not be opened  
      
         --------------------------------------------------------------------------------------- */
      
      static int readFileIntoStringVector( const std::string& fileName, StringVector& fileContents );

      /**---------------------------------------------------------------------------------------
         
         Replacement of sorts for strtok() (static method) (Simbios)
         Used to parse parameter file lines
      
         Should be moved to Utilities file
      
         @param lineBuffer           string to tokenize
         @param delimiter            token delimter
      
         @return number of args; if return value equals maxTokens, then more tokens than allocated
         
         --------------------------------------------------------------------------------------- */
      
      static char* strsep( char** lineBuffer, const char* delimiter );
      
      /**---------------------------------------------------------------------------------------
      
         Write file (helper method) (Simbios)
      
         @param lineVector                   line entries for file
         @param inputFileName                inputFileName
      
         @return SimTKOpenMMCommon::ErrorReturn if error -- else SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int writeFile( const StringVector& lineVector, const std::string& fileName );

      /**---------------------------------------------------------------------------------------
      
         Get statistics on an array
      
         @param numberOfEntries        number of entries in array
         @param array                  array
         @param average                average of array on output
         @param stdDev                 std dev of array on output
         @param minValue               min value in array on output
         @param minIndex               index of min value in array on output
         @param maxValue               max value in array on output
         @param maxIndex               index of max value in array on output
      
         if numberOfEntries <= 0, return 0 for all RealOpenMM values and -1 for index values
      
         @return SimTKOpenMMCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */
      
      static int getArrayStatistics( int numberOfEntries, const RealOpenMM* array,
                                     RealOpenMM* average, RealOpenMM* stdDev,
                                     RealOpenMM* minValue, int* minIndex,
                                     RealOpenMM* maxValue, int* maxIndex );
      
};
   
// ---------------------------------------------------------------------------------------

#endif // __SimTKOpenMMUtilities_H__
