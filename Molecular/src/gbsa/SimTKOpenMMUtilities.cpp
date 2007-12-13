
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

// class of shared, static utility methods

#include "SimTKOpenMMUtilities.h"
#include "SimTKOpenMMLog.h"

// fabs(), ...

#include <math.h>

/* ---------------------------------------------------------------------------------------

   Find distances**2 from a given atom (Simbios)

   @param atomCoordinates     atom coordinates
   @param atomIndex           atom index to find distances from
   @param numberOfAtoms       number of atoms
   @param distances           array of distances squared on return; array size must be at least
                              numberOfAtoms
   @param log                 if set, then print error messages to log file

   @return distances

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::getDistanceSquaredFromSpecifiedAtom( RealOpenMM** atomCoordinates, int atomIndex,
                                                               int numberOfAtoms, RealOpenMM* distances,
                                                               FILE* log ){

   // ---------------------------------------------------------------------------------------

   RealOpenMM atomXyz[3];
   // static const char* methodName    = "\nSimTKOpenMMUtilities::getDistanceSquaredFromSpecifiedAtom";

   // ---------------------------------------------------------------------------------------

   for( int jj = 0; jj < 3; jj++ ){
      atomXyz[jj] = atomCoordinates[atomIndex][jj];
   }
      
   return getDistanceSquaredFromSpecifiedPoint( atomCoordinates, atomXyz,
                                                numberOfAtoms, distances, log );
}

/* ---------------------------------------------------------------------------------------

   Find distances**2 from a given point (Simbios)

   @param atomCoordinates     atom coordinates
   @param point               point to find distances from
   @param numberOfAtoms       number of atoms
   @param distances           array of distances squared on return; array size must be at least \n
                              numberOfAtoms
   @param log                 if set, then print error messages to log file

   @return distances

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::getDistanceSquaredFromSpecifiedPoint( RealOpenMM** atomCoordinates,
                                                                RealOpenMM* point,
                                                                int numberOfAtoms,
                                                                RealOpenMM* distances, FILE* log ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nSimTKOpenMMUtilities::getDistanceSquaredFromSpecifiedPoint";

   // ---------------------------------------------------------------------------------------

   memset( distances, 0, sizeof( RealOpenMM )*numberOfAtoms );
   for( int ii = 0; ii < numberOfAtoms; ii++ ){
      for( int jj = 0; jj < 3; jj++ ){
         RealOpenMM diff = (point[jj] - atomCoordinates[ii][jj]);
         distances[ii] += diff*diff;
      }
   }
      
   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------

   Helper method to allocate RealOpenMM arrays (Simbios)

   @param bufferIndex         buffer index
   @param allocatedSz         array of allocated sizes
   @param bufferArray         array of allocated RealOpenMM arrays
   @param requestedSize       requested size
   @param dataAction          action flag: -1 = free memory \n
                                            1 = zero memory

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::allocateRealOpenMMBufferArray( int bufferIndex, int* allocatedSz,
                                                         RealOpenMM** bufferArray,
                                                         int requestedSize, int dataAction ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName    = "\nSimTKOpenMMUtilities::allocateRealOpenMMBufferArray";

   // ---------------------------------------------------------------------------------------

   // clear data

   if( dataAction == -1 ){
      if( allocatedSz[bufferIndex] && bufferArray[bufferIndex] ){
         SimTKOpenMMUtilities::Xfree( "bufferArray", __FILE__, __LINE__, bufferArray[bufferIndex] );
         allocatedSz[bufferIndex] = 0;
         bufferArray[bufferIndex] = NULL;
      }
   }
 
   // return if requested size is less than allcated

   if( allocatedSz[bufferIndex] > requestedSize ){
      return SimTKOpenMMCommon::DefaultReturn;
   }

   // free space if currently allocated

   if(  allocatedSz[bufferIndex] && bufferArray[bufferIndex] ){
      SimTKOpenMMUtilities::Xfree( "bufferArray", __FILE__, __LINE__, bufferArray[bufferIndex] );
   }

   // allocate

   // bufferArray[bufferIndex] = (RealOpenMM*) malloc( requestedSize*sizeof( RealOpenMM ) );
   bufferArray[bufferIndex] = (RealOpenMM*) SimTKOpenMMUtilities::Xmalloc( "bufferArray", __FILE__, __LINE__, requestedSize*sizeof( RealOpenMM ) );
   allocatedSz[bufferIndex] = requestedSize;

   // zero?

   if( dataAction == 1 ){
      memset( bufferArray[bufferIndex], 0, requestedSize*sizeof( RealOpenMM ) );
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------

   Print atom coordinates, ...

   @param numberAtoms         numberAtoms
   @param atomCoordinates     atomCoordinates (may be NULL)
   @param numberOf1Darrays    number of 1-d arrays (may be 0)
   @param oneDArrays          1-d arrays 
   @param idString            id string to be printed if set
   @param log                 print messages to log file

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::printCoordinateAnd1DArrays( int numberAtoms, RealOpenMM** atomCoordinates,
                                                      int numberOf1Darrays,
                                                      RealOpenMM** oneDArrays, 
                                                      const char* idString, FILE* log ){

   // ---------------------------------------------------------------------------------------

   if( log == NULL ){
      return SimTKOpenMMCommon::DefaultReturn;
   }

   if( idString ){
      (void) fprintf( log, "\n%s", idString );
   }
   for( int ii = 0; ii < numberAtoms; ii++ ){
      if( atomCoordinates != NULL ){
         (void) fprintf( log, "\n%d %12.4f %12.4f %12.4f",
                         ii + 1, atomCoordinates[ii][0], atomCoordinates[ii][1], atomCoordinates[ii][2] );
      } else {
         (void) fprintf( log, "\n" );
      }
      for( int jj = 0; jj < numberOf1Darrays; jj++ ){
         (void) fprintf( log, " %12.4f", oneDArrays[jj][ii] );
      }
   }

   (void) fflush(log);

   // ---------------------------------------------------------------------------------------
	
	return SimTKOpenMMCommon::DefaultReturn;	
}

/* ---------------------------------------------------------------------------------------

   Free array of strings

   @param arraySz             atom index
   @param arrayOfStrings      array of strings

   @return SimTKOpenMMCommon::DefaultReturn

   Used for freeing an array of strings

   @see getAtomIdStringGivenAtomIndex()

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::freeArrayOfStrings( int arraySz, char** arrayOfStrings ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::freeArrayOfStrings";

   // ---------------------------------------------------------------------------------------

   // free memory if allocated

   for( int ii = 0; ii < arraySz; ii++ ){
      if( arrayOfStrings[ii] ){
         SimTKOpenMMUtilities::Xfree( "atomIdStrings", __FILE__, __LINE__,  arrayOfStrings[ii] );
      }
   }
   SimTKOpenMMUtilities::Xfree( "arrayOfStrings", __FILE__, __LINE__,  arrayOfStrings );

   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------

   Tab string in place

   @param string   string to tab; assume string is of at least length=tab + 1
   @param tab      tab spacing length

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::tabStringInPlace( char* string, int tab ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::tabStringInPlace";

   // ---------------------------------------------------------------------------------------

   for( int ii = (int) strlen( string ); ii < tab; ii++ ){
      string[ii] = ' ';
   }
   string[tab] = '\0';

   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------

   Write debug fields (Simbios)

   @param numberOfFields       number of fields to print
   @param fields               fields
   @param numberOfStringFields number of string fields to print
   @param stringFields         string fields
   @param comment              comment (optinal -- ignored if NULL)
   @param debugFileName        output debug file name
   @param action               0 open file and return w/o printing \n
                               1 open file and print \n
                               2 close file (no print) \n
   @param log                  if set, then print error messages to log file

   @return debugFile unless file coud not be opened (or other errors )
   or file is closed -- for these cases return NULL

   stringFields printed after RealOpenMM fields

   --------------------------------------------------------------------------------------- */

FILE* SimTKOpenMMUtilities::writeDebugFile( int numberOfFields, const RealOpenMM* fields,
                                            int numberOfStringFields,
                                            const StringVector& stringFields,
                                            const char* comment, const char* debugFileName, int action,
                                            FILE* debugFile, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nSimTKOpenMMUtilities::writeDebugFile";

   // ---------------------------------------------------------------------------------------

   if( debugFileName != NULL && (action <= WriteDebugFile || debugFile == NULL) ){

   // open file

#ifdef WIN32
      fopen_s( &debugFile, debugFileName, "w" );
#else
      debugFile = fopen( debugFileName, "w" );
#endif

      if( debugFile != NULL ){
         if( log != NULL ){
            (void) fprintf( log, "%s opened file=<%s>.", methodName, debugFileName );
            (void) fflush( log );
         }
      } else {
         if( log != NULL ){
            (void) fprintf( log, "%s could not open file=<%s> -- abort output.",
                            methodName, debugFileName );
            (void) fflush( log );
         }
         return NULL;
      }

      if( action == OpenDebugFile ){
         return debugFile;
      }

   } else if( action == CloseDebugFile ){

   // close file

      if( debugFile ){
         if( log != NULL ){
            (void) fprintf( log, "%s closed debug file=<%s>.", 
                            methodName, debugFileName == NULL ? "??" : debugFileName );
            (void) fflush( log );
         }
         (void) fclose( debugFile );
      }
      return NULL;
   }   

   if( comment != NULL ){
      (void) fprintf( debugFile, "%s", comment );
   }

   if( numberOfFields > 0 || numberOfStringFields > 0 ){
      (void) fprintf( debugFile, "\n" );
      if( numberOfFields > 0 || fields != NULL ){
         for( int ii = 0; ii < numberOfFields; ii++ ){
            (void) fprintf( debugFile, "%.5e ", fields[ii] );
         }
      }
      if( numberOfStringFields > 0 ){
         for( StringVectorCI ii = stringFields.begin(); ii != stringFields.end() ; ii++ ){
            (void) fprintf( debugFile, "%s ", (*ii).c_str()  );
         }
      }
   }

   (void) fflush( debugFile );

   return debugFile;
}

/* ---------------------------------------------------------------------------------------

   Allocate 2D RealOpenMM array (Simbios)

   array[i][j]

   @param iSize                i-dimension
   @param jSize                j-dimension
   @param array2D              array (if null on entry allocated)
   @param initialize           if true, then initialize array
   @param initialValue         intitial value
   @param idString             id string

   @return allocated array

   --------------------------------------------------------------------------------------- */

RealOpenMM** SimTKOpenMMUtilities::allocateTwoDRealOpenMMArray( int iSize, int jSize,
                                                                RealOpenMM** array2D, 
                                                                int initialize, RealOpenMM initialValue,
                                                                const std::string& idString ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::allocate2DRealOpenMMArray";

   // ---------------------------------------------------------------------------------------

   if( array2D == NULL ){

      array2D     = (RealOpenMM**) SimTKOpenMMUtilities::Xmalloc( idString.c_str(), __FILE__, __LINE__, iSize*sizeof( RealOpenMM* ) );
      std::string blockString = idString;
      blockString.append( "Block" );

      RealOpenMM* block = (RealOpenMM*)  SimTKOpenMMUtilities::Xmalloc( blockString.c_str(),   __FILE__, __LINE__, jSize*iSize*sizeof( RealOpenMM ) );

      for( int ii = 0; ii < iSize; ii++ ){
         array2D[ii]  = block;
         block       += jSize;
      }    
   }

   if( initialize ){
      initialize2DRealOpenMMArray( iSize, jSize, array2D, initialValue );
   }

   return array2D;
}

/* ---------------------------------------------------------------------------------------

   Free 2D RealOpenMM array (Simbios)

   array[i][j]

   @param array2D              array (if null on entry allocated)
   @param idString             id string

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::freeTwoDRealOpenMMArray( RealOpenMM** array2D, const std::string& idString ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::allocate2DRealOpenMMArray";

   // ---------------------------------------------------------------------------------------

   if( array2D != NULL ){

      std::string blockString = idString;
      blockString.append( "Block" );

      SimTKOpenMMUtilities::Xfree( blockString.c_str(), __FILE__, __LINE__, array2D[0] );
      SimTKOpenMMUtilities::Xfree( idString.c_str(), __FILE__, __LINE__, array2D );
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------

   Initialize 2D RealOpenMM array (Simbios)

   array[i][j]

   @param iSize                i-dimension
   @param jSize                j-dimension
   @param array2D              array (if null on entry allocated)
   @param initialValue         intitial value

   @return array

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::initialize2DRealOpenMMArray( int iSize, int jSize,
                                                       RealOpenMM** array2D, RealOpenMM initialValue ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::initialize2DRealOpenMMArray";

   // ---------------------------------------------------------------------------------------

   bool useMemset;
   bool useMemsetSingleBlock;

   if( initialValue == 0.0f ){
      useMemset = true;
      if( jSize > 1 && (array2D[0] + jSize) == array2D[1] ){  
         useMemsetSingleBlock = true;
      } else {
         useMemsetSingleBlock = false;
      }

   } else {
      useMemset = false;
   }

   if( useMemset ){
      if( useMemsetSingleBlock ){
         memset( array2D[0], 0, iSize*jSize*sizeof( RealOpenMM ) );
      } else {
         for( int ii = 0; ii < iSize; ii++ ){
            memset( array2D[ii], 0, jSize*sizeof( RealOpenMM ) );
         }
      }
   } else {
      for( int ii = 0; ii < iSize; ii++ ){
         for( int jj = 0; jj < jSize; jj++ ){
            array2D[ii][jj] = initialValue;
         }
      }
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/* ---------------------------------------------------------------------------------------
      
   Malloc memory of size bytesToAllocate and zero
      
   @param bytesToAllocate      bytes to allocate
      
   @return ptr to allocated memory; NULL if bytesToAllocate <= 0
      
   --------------------------------------------------------------------------------------- */
          
char* SimTKOpenMMUtilities::allocateAndZero( unsigned int bytesToAllocate ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::allocateAndZero";

   // ---------------------------------------------------------------------------------------

   if( bytesToAllocate <= 0 ){
      return NULL;
   }

   char* ptrToMemory = (char*) SimTKOpenMMUtilities::Xmalloc( "ptrToMemory", __FILE__, __LINE__, bytesToAllocate*sizeof( char ) );
   memset( ptrToMemory, 0, bytesToAllocate );

   return ptrToMemory;
}

/* ---------------------------------------------------------------------------------------

   Normalize 3-vector -- helper method

   @param vector 3-vector to normalize

   --------------------------------------------------------------------------------------- */
     
void SimTKOpenMMUtilities::normalizeVector3( RealOpenMM* vector ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::normalizeVector3";

   // ---------------------------------------------------------------------------------------

   RealOpenMM sum   = vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
   sum              = sum > 0.0 ? (RealOpenMM) (1.0/SQRT( sum )) : (RealOpenMM) 0.0;

   vector[0]       *= sum;
   vector[1]       *= sum;
   vector[2]       *= sum;

   return;
}

/* ---------------------------------------------------------------------------------------

   Remove 3-vector -- helper method

   @param vectorToRemove      vector to remove
   @param vector              vector to from which 'vectorToRemove' is to be removed
                              vector is normalized after the component is subtracted out

   --------------------------------------------------------------------------------------- */
     
void SimTKOpenMMUtilities::removeVector3( RealOpenMM* vectorToRemove, RealOpenMM* vector ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::removeVector3";

   // ---------------------------------------------------------------------------------------

   RealOpenMM dot   = vectorToRemove[0]*vector[0] + vectorToRemove[1]*vector[1] + vectorToRemove[2]*vector[2];

   vector[0]       -= dot*vectorToRemove[0];
   vector[1]       -= dot*vectorToRemove[1];
   vector[2]       -= dot*vectorToRemove[2];

   normalizeVector3( vector );
}

/* ---------------------------------------------------------------------------------------

   Compute cross product of two 3-vectors and place in 3rd vector  -- helper method

   vectorZ = vectorX x vectorY

   @param vectorX             x-vector
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
     
void SimTKOpenMMUtilities::crossProductVector3( RealOpenMM* vectorX, RealOpenMM* vectorY, RealOpenMM* vectorZ ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::crossProductVector3";

   // ---------------------------------------------------------------------------------------

   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];

   return;
}

/**---------------------------------------------------------------------------------------

   Compute matrix product of 3x3 matrix and 3-vector and place in 3rd vector  -- helper method

   vectorZ = matrixX . vectorY

   @param matrixX             matrixX
   @param vectorY             y-vector
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
    
void SimTKOpenMMUtilities::matrixProductVector3( RealOpenMM* matrixX, RealOpenMM* vectorY,
                                                 RealOpenMM* vectorZ ){
     
   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::matrixProductVector3";

   // ---------------------------------------------------------------------------------------

   vectorZ[0]  = matrixX[0]*vectorY[0] + matrixX[3]*vectorY[1] + matrixX[6]*vectorY[2];
   vectorZ[1]  = matrixX[1]*vectorY[0] + matrixX[4]*vectorY[1] + matrixX[7]*vectorY[2];
   vectorZ[2]  = matrixX[2]*vectorY[0] + matrixX[5]*vectorY[1] + matrixX[8]*vectorY[2];

   return;
}

/**---------------------------------------------------------------------------------------

   Compute cross product between two 3x3 matrices

   @param vectorZ = matrixX x matrixY

   @param matrixX             matrixX
   @param matrixY             matrixY
   @param vectorZ             z-vector

   @return vector is vectorZ

   --------------------------------------------------------------------------------------- */
    
void SimTKOpenMMUtilities::matrixCrossProductMatrix3( RealOpenMM* matrixX, RealOpenMM* matrixY,
                                                      RealOpenMM* vectorZ ){

   // ---------------------------------------------------------------------------------------

   // static const int indices[3][2] = { { 3, 6 }, { 6, 0 }, { 0 , 3 } };
   // static const char* methodName = "\nSimTKOpenMMUtilities::matrixCrossProductMatrix3";
   RealOpenMM* xPtr[3];
   RealOpenMM* yPtr[3];

   // ---------------------------------------------------------------------------------------

   xPtr[0] = matrixX;
   xPtr[1] = matrixX + 3;
   xPtr[2] = matrixX + 6;

   yPtr[0] = matrixY;
   yPtr[1] = matrixY + 3;
   yPtr[2] = matrixY + 6;

   vectorZ[0] = DOT3( xPtr[1], yPtr[2] ) - DOT3( xPtr[2], yPtr[1] );
   vectorZ[1] = DOT3( xPtr[2], yPtr[0] ) - DOT3( xPtr[0], yPtr[2] );
   vectorZ[2] = DOT3( xPtr[0], yPtr[1] ) - DOT3( xPtr[1], yPtr[0] );

   return;
}

/* ---------------------------------------------------------------------------------------

   Centralized malloc/new

   @param name                ptr name
   @param fileName            file name
   @param line                file line no.
   @param file line           size in bytes to be allocated

   @return ptr to allocated object

   --------------------------------------------------------------------------------------- */
     
void* SimTKOpenMMUtilities::Xmalloc( const char* name, char* fileName, int line, unsigned int size ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::Xmalloc";

   // ---------------------------------------------------------------------------------------

#ifdef UseGromacsMalloc
   return save_malloc( name, fileName, line, size );
#else
   return (void*) new char[size];
#endif

}

/* ---------------------------------------------------------------------------------------

   Centralized free/delete

   @param name                ptr name
   @param fileName            file name
   @param line                file line no.
   @param ptr                 ptr to be freed

   --------------------------------------------------------------------------------------- */
     
void SimTKOpenMMUtilities::Xfree( const char* name, char* fileName, int line, void* ptr ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::Xfree";

   // ---------------------------------------------------------------------------------------

#ifdef UseGromacsMalloc
   return save_free( name, fileName, line, ptr );
#else
   delete (char*) ptr;
   return;
#endif
}

/* ---------------------------------------------------------------------------------------
      
   Format array of reals
      
   @param message             input string stream
   @param realArray           array of RealOpenMMs
   @param numberOfFields      number of fields (optional - defaults to 3)
   @param factor					scale factor
      
   @return SimTKOpenMMCommon::DefaultReturn

--------------------------------------------------------------------------------------- */
          
int SimTKOpenMMUtilities::formatRealStringStream( std::stringstream& message, const RealOpenMM* realArray,
   	       	                                   int numberOfFields, RealOpenMM factor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nSimTKOpenMMUtilities::Xfree";

   // ---------------------------------------------------------------------------------------

   for( int ii = 0; ii < numberOfFields; ii++ ){ 
      message << factor*realArray[ii] << " ";
   }
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param lineBuffer           string to tokenize
   @param tokenArray           upon return vector of tokens
   @param delimiter            token delimter

   @return number of args

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::tokenizeString( char* lineBuffer, StringVector& tokenArray,
                                          const std::string delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   // ---------------------------------------------------------------------------------------

// (void) fprintf( stdout, "\nIn SimTKOpenMMUtilities::tokenizeString <%s>", lineBuffer );
// (void) fflush( stdout );

   char *ptr_c = NULL;

   for( ; (ptr_c = SimTKOpenMMUtilities::strsep( &lineBuffer, delimiter.c_str() )) != NULL; ){
      if( *ptr_c ){
         tokenArray.push_back( std::string( ptr_c ) );
      }
   }

   return (int) tokenArray.size();
}

/**---------------------------------------------------------------------------------------

   Local version of strncasecmp (missing in Windows) (static method) (Simbios)

   @param string1                 first string
   @param string2                 second string
   @param matchLength             match length

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::localStrncasecmp( const char *string1, const char *string2, int matchLength ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::localStrncasecmp"

   char ch1,ch2;

   // ---------------------------------------------------------------------------------------

   if( matchLength == 0 ){ 
      return SimTKOpenMMCommon::DefaultReturn;
   }

   do {   
      ch1 = toupper(*(string1++));
      ch2 = toupper(*(string2++));
      if( ch1 != ch2 )return (ch1-ch2);
      matchLength--;
   } while( ch1 && matchLength ); 

   return SimTKOpenMMCommon::DefaultReturn;  

}

/**---------------------------------------------------------------------------------------

   Check that string is valid integer

   @param stringToCheck string to check

   @return true if string is a valid integer

   --------------------------------------------------------------------------------------- */

bool SimTKOpenMMUtilities::isValidInteger( std::string stringToCheck ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::isValidInteger";

   // ---------------------------------------------------------------------------------------

   int ii;

   return checkString<int>(ii, stringToCheck, std::dec );
}

/**---------------------------------------------------------------------------------------

   Check that string is valid RealOpenMM

   @param stringToCheck string to check

   @return true if string is a valid RealOpenMM

   --------------------------------------------------------------------------------------- */

bool SimTKOpenMMUtilities::isValidRealOpenMM( std::string stringToCheck ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::isValidRealOpenMM";

   // ---------------------------------------------------------------------------------------

   RealOpenMM ii;

   return checkString<RealOpenMM>(ii, stringToCheck, std::dec );
}

/**---------------------------------------------------------------------------------------

   Read file into string vector (Simbios) 

   @param fileName       file name
   @param fileContents   string vector containing file contents upon return
                         one string per line

   @return SimTKOpenMMCommon::DefaultReturn unless file could not be opened  

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::readFileIntoStringVector( const std::string& fileName,
                                                    StringVector& fileContents ){

   // ---------------------------------------------------------------------------------------

   static int bufferSize = 2048;
   static char buffer[2048];
   static const std::string methodName = "\nSimTKOpenMMUtilities::readFileIntoStringVector";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* file = NULL;
#ifdef WIN32
   fopen_s( &file, fileName.c_str(), "r" );
#else
   file = fopen( fileName.c_str(), "r" );
#endif

   if( file != NULL ){
      std::stringstream message;
      message << methodName.c_str() << " opened file=<" <<  fileName.c_str() << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message << methodName.c_str() << " could not open file=<" <<  fileName.c_str() << ">.";
(void) fprintf( stderr, "\n%s\n", message.str().c_str() ); 
(void) fflush( stderr );
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;  
   }

   // loop over all lines in file, checking for end-of-file

   int lineNumber = 0;

   while( !feof( file ) ){ 

      // read next line

      int bufferLen;
      lineNumber++;
      if( fgets( buffer, bufferSize, file ) != NULL ){
         bufferLen = (int) strlen( buffer );
         if( bufferLen > 0 ){
            buffer[bufferLen-1] = '\0';
         }
         // (void) fprintf( log, "%s", buffer );
         // (void) fflush( log );
         fileContents.push_back( buffer );
      }
   }

   // done

   (void) fclose( file );

   if( file ){
      //std::stringstream message;
      //message << methodName.c_str() << " read " << lineNumber << " lines from file=<" << fileName->c_str() << ">.";
      // AmoebaLog::printMessage( message );
   }

   return SimTKOpenMMCommon::DefaultReturn;  
}

/**---------------------------------------------------------------------------------------

   Tokenize a string (static method) (Simbios)

   @param line                 string to tokenize
   @param tokenVector          upon return vector of tokens
   @param delimiter            token delimter
   @param clearTokenVector     if true, clear tokenVector

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::tokenizeString( const std::string& line, StringVector& tokenVector,
                                          const std::string& delimiter, int clearTokenVector ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::tokenizeString";

   static int bufferSz       = 2048;
   static char* lineBuffer   = NULL;

   // ---------------------------------------------------------------------------------------

   char *ptr_c;

   // clear tokeen vector
   if( clearTokenVector ){
      tokenVector.clear();
   }   

   // allocate space for line buffer and copy via sprintf()

   if( lineBuffer == NULL || bufferSz < (int) line.size() ){
      if( lineBuffer != NULL ){
         free( lineBuffer );
      }
      if( bufferSz < (int) line.size() ){
         bufferSz = (int) (2*line.size());
      } 
      lineBuffer = (char*) malloc( bufferSz*sizeof( char) );
   }

#ifdef WIN32
   (void) sprintf_s( lineBuffer, bufferSz, "%s", line.c_str() );
#else
   (void) sprintf( lineBuffer, "%s", line.c_str() );
#endif

   // parse

   while( (ptr_c = SimTKOpenMMUtilities::strsep( &lineBuffer, delimiter.c_str() )) != NULL ){
      if( *ptr_c ){
         tokenVector.push_back( std::string( ptr_c ) );
      }   
   }

   return SimTKOpenMMCommon::DefaultReturn;  
}

/**---------------------------------------------------------------------------------------

   Replacement of sorts for strtok() (static method) (Simbios)
   Used to parse parameter file lines

   Should be moved to Utilities file

   @param lineBuffer           string to tokenize
   @param delimiter            token delimter

   @return number of args; if return value equals maxTokens, then more tokens than allocated

   --------------------------------------------------------------------------------------- */

char* SimTKOpenMMUtilities::strsep( char** lineBuffer, const char* delimiter ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nTinkerParameterSet::strsep"

   char *s;
   const char *spanp;
   int c, sc;
   char *tok;

   // ---------------------------------------------------------------------------------------

   s = *lineBuffer;
   if( s == NULL ){
      return (NULL);
   }

   for( tok = s;; ){
      c = *s++;
      spanp = delimiter;
      do {
         if( (sc = *spanp++) == c ){
            if( c == 0 ){
               s = NULL;
            } else {
               s[-1] = 0;
            }
            *lineBuffer = s;
            return( tok );
         }
      } while( sc != 0 );
   }
}

/**---------------------------------------------------------------------------------------

   Return lower case copy of string

   @param string                  string

   @return lower cased string

   --------------------------------------------------------------------------------------- */

//int SimTKOpenMMUtilities::toLower( std::string& string ){

   // ---------------------------------------------------------------------------------------

   // static const std::string methodName = "\nSimTKOpenMMUtilities::toLower"

   // ---------------------------------------------------------------------------------------

	// transform string to lower case
	    
	// std::transform( string.begin(), string.end(), string.begin(), (int(*)(int)) std::tolower);
	//return SimTKOpenMMCommon::DefaultReturn;
		
//}

/**---------------------------------------------------------------------------------------

   Write file (helper method) (Simbios)

   @param lineVector                   line entries for file
   @param inputFileName                inputFileName

   @return SimTKOpenMMCommon::ErrorReturn if error -- else SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMUtilities::writeFile( const StringVector& lineVector, const std::string& fileName ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMUtilities::writeFile";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* file = NULL;
#ifdef WIN32
   fopen_s( &file, fileName.c_str(), "w" );
#else
   file = fopen( fileName.c_str(), "w" );
#endif

   if( file != NULL ){
      std::stringstream message;
      message <<  methodName.c_str() << " opened file=<" <<  fileName.c_str() << ">.";
      SimTKOpenMMLog::printMessage( message );
   } else {
      std::stringstream message;
      message <<  methodName.c_str() << " could not open file=<" <<  fileName.c_str() << ">.";
      SimTKOpenMMLog::printMessage( message );
      return SimTKOpenMMCommon::ErrorReturn;
   }

   // ---------------------------------------------------------------------------------------

   // loop over lines

   int bodyLines = (int) lineVector.size();
   for( StringVectorCI ii = lineVector.begin(); ii != lineVector.end(); ii++ ){
      if( (*ii).length() > 1 ){
         (void) fprintf( file, "%s\n", (*ii).c_str() );
      } else {
         bodyLines--;
      }
   }
   (void) fflush( file );
   (void) fclose( file );

   std::stringstream message;
   message << " :: closed file=<" <<  fileName.c_str() << "> number of body lines in file=" << bodyLines;
   SimTKOpenMMLog::printMessage( message );

   return SimTKOpenMMCommon::DefaultReturn;
}

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

int SimTKOpenMMUtilities::getArrayStatistics( int numberOfEntries, const RealOpenMM* array, RealOpenMM* average,
                                              RealOpenMM* stdDev, RealOpenMM* minValue, int* minIndex,
                                              RealOpenMM* maxValue, int* maxIndex ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\nSimTKOpenMMUtilities::getArrayStatistics";
   static const RealOpenMM zero        = 1.0;
   static const RealOpenMM one         = 1.0;

   // ---------------------------------------------------------------------------------------

   if( numberOfEntries <= 0 ){
      *average  = *stdDev   = *minValue = *maxValue = zero;
      *maxIndex = *minIndex = -1;
      return  SimTKOpenMMCommon::DefaultReturn;
   }
   
   *average  = array[0];
   *stdDev   = (*average)*(*average);
   *minValue = array[0];
   *maxValue = array[0];
   *maxIndex = *minIndex = 0;
   for( int ii = 1; ii < numberOfEntries; ii++ ){

      *average += array[ii];
      *stdDev  += array[ii]*array[ii];

      if( array[ii] < *minValue ){
         *minValue  = array[ii];
         *minIndex = ii;
      }
      if( array[ii] > *maxValue ){
         *maxValue = array[ii];
         *maxIndex = ii;
      }
   }
  
   RealOpenMM numberOfEntriesR  = (RealOpenMM) numberOfEntries;
   *average              /= numberOfEntriesR;
   *stdDev                = *stdDev - numberOfEntriesR*(*average)*(*average);

   if( numberOfEntriesR > one ){
      *stdDev = SQRT( (*stdDev)/(numberOfEntriesR - one) );
   } 
      
   return  SimTKOpenMMCommon::DefaultReturn;
}
