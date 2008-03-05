#ifndef SimTK_ELEMENTREP_H_
#define SimTK_ELEMENTREP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
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

#include "SimTKcommon.h"

#include "simbody/internal/Element.h"

namespace SimTK {

class ElementRep : public PIMPLImplementation<Element,ElementRep> {
public:
    ElementRep* clone() const {return new ElementRep(*this);}

    ElementRep()
        : atomicNumber(Element::InvalidAtomicNumber)
    {}

    ElementRep(int an, Element::Name nm, Element::Symbol sm, Mass ms)
        : atomicNumber(an), name(nm), symbol(sm), averageMass(ms)
    {}

    Element::Symbol getSymbol() const {return symbol;}
    Element::Name   getName() const {return name;}

private:
    Element::Symbol symbol;
    Element::Name name;
    // Average mass is used, as opposed to typical mass
    Mass averageMass;
    int atomicNumber;
};

} // namespace SimTK

#endif // SimTK_ELEMENTREP_H_
