/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * This file contains the code used to build the various constraint types.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "ConstraintNode.h"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;

///////////////////////////////////////////////
// Implementation of ConstraintNode methods. //
///////////////////////////////////////////////

/*
 * How to specify a constraint equation:
 *
 * Required info:
 *
 *    Dependencies:
 *      - constraint level: position, velocity, acceleration
 *      - has time dependence?
 *    A list of bodies
 *      - this is just those bodies to which constraint
 *        forces & torques are *directly* applied to 
 *        enforce the constraint (can be spatial force
 *        or mobility force)
 *
 * All constraints:
 *    Acceleration error (given udot)
 *    Constraint forces (given lambda)
 * Position or velocity level:
 *    Velocity error (given u)
 * Position level only:
 *    Position error (given q)
 */
