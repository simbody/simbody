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
 * Implementations of State and StateRep.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/State.h"

#include "StateRep.h"

#include <string>
#include <iostream> 
#include <sstream>


namespace simtk {

State::State()
  : rep(new StateRep()) {
}

State::State(int nq, int nu)
  : rep(new StateRep()) {
    rep->resize(nq,nu);
}

State::~State() {
    delete rep;
}

State::State(const State& src) : rep(0) {
    if (src.rep)
        rep = src.rep->clone();
}

State& State::operator=(const State& src) {
    if (&src != this) {
        delete rep; rep=0;
        if (src.rep)
            rep = src.rep->clone();
    }
    return *this;
}


const Vector& State::getQ() const {
    return rep->getConfiguration();
}
const Vector& State::getU() const{
    return rep->getMotion();
}

Vector& State::updQ()  {
    return rep->updConfiguration();
}
Vector& State::updU() {
    return rep->updMotion();
}

} // namespace simtk
