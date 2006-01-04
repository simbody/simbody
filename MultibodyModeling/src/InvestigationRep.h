#ifndef SIMTK_SIMBODY_INVESTIGATION_REP_H_
#define SIMTK_SIMBODY_INVESTIGATION_REP_H_

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
 * Declaration of the library side (rep) of the client Investigation class. Rebuilding this
 * code does not require recompiling client-side code, and if the library is a DLL
 * (shared library) then the client doesn't require relinking either.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Model.h"
#include "simbody/internal/State.h"
#include "simbody/internal/Investigation.h"

namespace simtk {

class InvestigationRep {
public:
    InvestigationRep() : model(0) { }
    InvestigationRep(const Model& m) : model(new Model(m)) {
        current = model->getDefaultState();
    }
    virtual ~InvestigationRep() { }

    virtual InvestigationRep* clone() const = 0;

private:
    Model*  model;
    State   current;
};


} // namespace simtk


#endif // SIMTK_SIMBODY_INVESTIGATION_REP_H_
