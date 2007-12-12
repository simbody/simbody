
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

#include "SimTKOpenMMLog.h"

// static settings

SimTKOpenMMLog* SimTKOpenMMLog::_simTKOpenMMLog = NULL;

/**---------------------------------------------------------------------------------------

   SimTKOpenMMLog constructor (Simbios)

   @param logFile file reference for logging

   --------------------------------------------------------------------------------------- */

SimTKOpenMMLog::SimTKOpenMMLog( FILE* logFile ) : _logFile( logFile ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::SimTKOpenMMLog";

// ---------------------------------------------------------------------------------------

   _logLevel  = LogLowLevel;

}

/**---------------------------------------------------------------------------------------

   SimTKOpenMMLog destructor (Simbios)

   --------------------------------------------------------------------------------------- */

SimTKOpenMMLog::~SimTKOpenMMLog( ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::~SimTKOpenMMLog";

// ---------------------------------------------------------------------------------------

}

/**---------------------------------------------------------------------------------------

   Get LogFile

   @return    logFile

   --------------------------------------------------------------------------------------- */

FILE* SimTKOpenMMLog::getLogFile( void ) const {

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::getLogFile";

// ---------------------------------------------------------------------------------------

   return _logFile;

}

/**---------------------------------------------------------------------------------------

   Set LogFile

   @param    input logFile

   @return   SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMLog::setLogFile( FILE* logFile ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::setLogFile";

// ---------------------------------------------------------------------------------------

   _logFile = logFile;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Set LogLevel

   @param    input logLevel

   @return   SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMLog::setLogLevel( SimTKOpenMMLog::LogLevels logLevel ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::setLogLevel";

// ---------------------------------------------------------------------------------------

   _logLevel = logLevel;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   SimTKOpenMMLog log message to log (Simbios)

   @param message         message to log

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMLog::logMessage( const std::stringstream& message ) const {

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::logMessage";

// ---------------------------------------------------------------------------------------

   if( _logFile ){

//      (void) fprintf( stderr, "%s", message.str().c_str() );
//      (void) fflush( stderr );

      (void) fprintf( _logFile, "%s", message.str().c_str() );
      (void) fflush( _logFile );

   }
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Set global simTKOpenMMLog (Simbios)

   @param logFile         file to log to

   @return new SimTKOpenMMLog

   --------------------------------------------------------------------------------------- */

SimTKOpenMMLog* SimTKOpenMMLog::setSimTKOpenMMLog( FILE* logFile ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::setSimTKOpenMMLog";

// ---------------------------------------------------------------------------------------

   // allow for multiple logs

/*
   if( _simTKOpenMMLog ){
      delete _simTKOpenMMLog;
   }
*/
   _simTKOpenMMLog = new SimTKOpenMMLog( logFile );

   return _simTKOpenMMLog;
}

/**---------------------------------------------------------------------------------------

   Get global simTKOpenMMLog -- static method (Simbios)

   @return static member

   --------------------------------------------------------------------------------------- */

SimTKOpenMMLog* SimTKOpenMMLog::getSimTKOpenMMLog( void ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::getSimTKOpenMMLog";

// ---------------------------------------------------------------------------------------

   if( !_simTKOpenMMLog ){
      _simTKOpenMMLog = new SimTKOpenMMLog( );
   }
   return _simTKOpenMMLog;
}

/**---------------------------------------------------------------------------------------

   Get global simTKOpenMMLog (Simbios)

   @return FILE reference

   --------------------------------------------------------------------------------------- */

FILE* SimTKOpenMMLog::getSimTKOpenMMLogFile( void ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::getSimTKOpenMMLogFile";

// ---------------------------------------------------------------------------------------

   SimTKOpenMMLog* simTKOpenMMLog = getSimTKOpenMMLog( );
   return simTKOpenMMLog->getLogFile();
}

/**---------------------------------------------------------------------------------------

   Static method to print message (Simbios)

   @param message         message to log

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMLog::printMessage( const std::stringstream& message ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::printMessage";

// ---------------------------------------------------------------------------------------

   if( _simTKOpenMMLog ){
      _simTKOpenMMLog->logMessage( message );
   } else {
      (void) fprintf( stderr, "%s", message.str().c_str() );
      (void) fflush( stderr );
   }
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Static method to print warning message (Simbios)
   If global _simTKOpenMMLog is not set, then print to stderr

   @param message         message to log

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMLog::printWarning( const std::stringstream& message ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::printWarning";

// ---------------------------------------------------------------------------------------

   if( _simTKOpenMMLog ){
      std::stringstream messageWithHeader;
      messageWithHeader << "Warning: " << message.str();
      _simTKOpenMMLog->logMessage( messageWithHeader );
   } else {
      (void) fprintf( stderr, "Warning: %s", message.str().c_str() );
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Static method to print error message and exit program (Simbios)
   If global _simTKOpenMMLog is not set, then print to stderr

   @param message   message to log

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int SimTKOpenMMLog::printError( const std::stringstream& message ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nSimTKOpenMMLog::printError";

// ---------------------------------------------------------------------------------------

   if( _simTKOpenMMLog ){
      std::stringstream messageWithHeader;
      messageWithHeader << "Error: " << message.str();
      _simTKOpenMMLog->logMessage( messageWithHeader );
   } else {
      (void) fprintf( stderr, "Error: %s", message.str().c_str() );
      (void) fflush( stderr );
   }
   exit(-1);

   return SimTKOpenMMCommon::DefaultReturn;
}
