//-----------------------------------------------------------------------------
// File:     RotationTest.cpp
// Class:    None
// Parent:   None
// Purpose:  Test 3x3 Rotation class relating two right-handed orthogonal bases
//-----------------------------------------------------------------------------

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Paul Mitiguy                                                      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**@file
 * Tests for the classes defined in Rotation.h.
 */

//-----------------------------------------------------------------------------
#include "SimTKcommon.h"
#include <iostream>

//-----------------------------------------------------------------------------
using namespace SimTK;

bool  doRequiredTasks();
void  WriteStringToScreen( const char outputString[] )  { std::cout << outputString; }  

//-------------------------------------------------------------------
int main() {

   // Default value is program failed
   bool programSucceeded = false;

   // It is a good programming practice to do little in the main function of a program.
   // The try-catch code in this main routine catches exceptions thrown by functions in the
   // try block, e.g., catching an exception that occurs when a NULL pointer is de-referenced.
   try
   {
      // Do the required programming tasks
      programSucceeded = doRequiredTasks();
   }
   // This catch statement handles certain types of exceptions
   catch( const std::exception& e )
   {
      WriteStringToScreen( "\n\n Error: Programming error encountered.\n The exception thrown is: " );
      WriteStringToScreen( e.what() );
      WriteStringToScreen( "\n\n" );
   }
   // The exception-declaration statement (...) handles any type of exception,
   // including C exceptions and system/application generated exceptions.
   // This includes exceptions such as memory protection and floating-point violations.
   // An ellipsis catch handler must be the last handler for its try block.
   catch( ... )
   {
      WriteStringToScreen( "\n\n Error: Programming error encountered.\n An unhandled exception was thrown.\n\n" );
   }

   // The value returned by the main function is the exit status of the program.
   // A normal program exit returns 0 (other return values usually signal an error).
   return programSucceeded == true ? 0 : 1;
}


//-------------------------------------------------------------------
bool  doRequiredTasks( ) {
  
    return true;
}



