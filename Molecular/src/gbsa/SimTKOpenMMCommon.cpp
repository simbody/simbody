
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

// static settings

// initialization of static data members

const std::string SimTKOpenMMCommon::NotSet                  = std::string( "NotSet" );
const std::string SimTKOpenMMCommon::Comment                 = std::string( "#" );
const std::string SimTKOpenMMCommon::Tab                     = std::string( "\t" ); 

const int SimTKOpenMMCommon::DefaultReturn                   = 0;
const int SimTKOpenMMCommon::ErrorReturn                     = -1;
const RealOpenMM SimTKOpenMMCommon::BigCutoffValue           = 1.0e+05;

      // units

const int SimTKOpenMMCommon::MdUnits                         = 1;
const int SimTKOpenMMCommon::KcalAngUnits                    = 2;

// specify RealOpenMM number format
  
const int SimTKOpenMMCommon::HighStringStreamNumberWidth     = 20;   
const int SimTKOpenMMCommon::HighStringStreamNumberPrecision = 12;

const RealOpenMM SimTKOpenMMCommon::DegreeToRadians          = (RealOpenMM) 0.017453292;

