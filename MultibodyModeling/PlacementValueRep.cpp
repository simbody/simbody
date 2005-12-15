/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * Implementation of non-inline PlacementValueRep methods.
 */

#include "SimbodyCommon.h"
#include "PlacementRep.h"
#include "SubsystemRep.h"

#include <iostream> 
using std::cout;
using std::endl;

namespace simtk {


    // PLACEMENT VALUE REP //

void PlacementValueRep::checkPlacementValueConsistency(const Subsystem* expOwner, 
                                                       int              expIndexInOwner,
                                                       const Subsystem& expRoot) const
{
    cout << "CHECK PLACEMENT VALUE CONSISTENCY FOR PlacementValueRep@" << this << endl;
    if (!myHandle) 
        cout << "*** NO HANDLE ***" << endl;
    else if (myHandle->rep != this)
        cout << "*** Handle->rep=" << myHandle->rep << " which is *** WRONG ***" << endl;
    if (owner != expOwner)
        cout << "*** WRONG OWNER@" << owner << "; should have been " << expOwner << endl;
    if (indexInOwner != expIndexInOwner)
        cout << "*** WRONG INDEX " << indexInOwner << "; should have been " << expIndexInOwner << endl;


    if (expOwner == 0) {
        if (client)
          cout << "*** UNOWNED PLACEMENT VALUE HAD CLIENT ***" << endl;
    } else {
        if (!client) 
            cout << "*** NO CLIENT ***" << endl;
        else if (!client->getRep().hasValueSlot())
            cout << "*** CLIENT HAS NO VALUE SLOT??? ***" << endl;
        else if (&(client->getRep().getValueSlot().getRep()) != this)
            cout << "*** CLIENT HAS WRONG PLACEMENT VALUE SLOT@" 
                << &client->getRep().getValueSlot().getRep() << endl;
    }
}


} // namespace simtk
