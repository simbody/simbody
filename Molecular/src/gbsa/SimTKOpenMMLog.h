
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

#ifndef __SimTKOpenMMLog_H__
#define __SimTKOpenMMLog_H__

#include <stdio.h>
#include <sstream>
#include "SimTKOpenMMCommon.h"

/** ---------------------------------------------------------------------------------------

   SimTKOpenMMLog class used for logging

   --------------------------------------------------------------------------------------- */

class SimTKOpenMMLog {

   public:

      // log levels

      enum LogLevels { LogOff, LogLowLevel, LogHighLevel };

   private:

      // file to write to

      FILE* _logFile;

      // log level

      LogLevels _logLevel;

      // global reference

      static SimTKOpenMMLog* _simTKOpenMMLog;

   public:

      /**---------------------------------------------------------------------------------------
      
         SimTKOpenMMLog constructor (Simbios)
      
         @param logFile file reference for logging
      
         --------------------------------------------------------------------------------------- */
      
      SimTKOpenMMLog( FILE* logFile = NULL );

      /**---------------------------------------------------------------------------------------
      
         SimTKOpenMMLog destructor (Simbios)
      
         --------------------------------------------------------------------------------------- */
      
      ~SimTKOpenMMLog( );

      /**---------------------------------------------------------------------------------------
      
         SimTKOpenMMLog log message to log (Simbios)
      
         @param message         message to log
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
     int logMessage( const std::stringstream& message ) const;

      /**---------------------------------------------------------------------------------------

         Get LogFile

         @return    logFile
      
         --------------------------------------------------------------------------------------- */

      FILE* getLogFile( void ) const;

      /**---------------------------------------------------------------------------------------

         Set LogFile

         @param    input logFile

         @return   AmoebaCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setLogFile( FILE* logFile );

      /**---------------------------------------------------------------------------------------

         Set LogLevel

         @param    input logLevel

         @return   AmoebaCommon::DefaultReturn
      
         --------------------------------------------------------------------------------------- */

      int setLogLevel( SimTKOpenMMLog::LogLevels logLevel );

      /**---------------------------------------------------------------------------------------
      
         Set global simTKLog (Simbios)
      
         @param logFile         file to log to
      
         @return new SimTKOpenMMLog
      
         --------------------------------------------------------------------------------------- */
      
      static SimTKOpenMMLog* setSimTKOpenMMLog( FILE* logFile = NULL );
      
      /**---------------------------------------------------------------------------------------
      
         Get global simTKLog -- static method (Simbios)
      
         @return static member
      
         --------------------------------------------------------------------------------------- */
      
      static SimTKOpenMMLog* getSimTKOpenMMLog( void );
      
      /**---------------------------------------------------------------------------------------
      
         Get global simTKLog (Simbios)
      
         @return FILE reference
      
         --------------------------------------------------------------------------------------- */
      
      static FILE* getSimTKOpenMMLogFile( void );
      
      /**---------------------------------------------------------------------------------------
      
         Staitc method to print message (Simbios)
      
         @param message         message to log
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      static int printMessage( const std::stringstream& message );
      
      /**---------------------------------------------------------------------------------------
      
         Staitc method to print warning message (Simbios)
      
         @param message         message to log
      
         @return 0
      
         --------------------------------------------------------------------------------------- */
      
      static int printWarning( const std::stringstream& message );
      
      /**---------------------------------------------------------------------------------------
      
         Static method to print error message and exist program (Simbios)
      
         @param message         message to log
      
         @return 0
      
         --------------------------------------------------------------------------------------- */

     static int printError( const std::stringstream& message );

};

#endif //__SimTKOpenMMLog_H__
