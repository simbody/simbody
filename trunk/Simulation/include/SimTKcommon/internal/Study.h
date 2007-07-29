#ifndef SimTK_SimTKCOMMON_STUDY_H_
#define SimTK_SimTKCOMMON_STUDY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"

namespace SimTK {

/**
 * The handle class which serves as the abstract parent of all Studies.
 *
 * TODO!
 *
 * There are two distinct users of this class:
 *   - Study Users: people who are making use of a concrete Study (which will
 *     inherit methods from this class)
 *   - Study Developers: people who are writing concrete Study classes
 *
 * Only methods intended for Study Users and a few bookkeeping methods
 * are in the main Study class, which is a SimTK Handle class, meaning
 * that it consists only of a single pointer, which points to a 
 * Study::Guts class. The Guts class is abstract, and virtual methods
 * to be implemented by Study Developers in the concrete
 * Study are defined there, along
 * with other utilities of use to the concrete Study Developer but
 * not to the end user. The Guts class is declared in a separate
 * header file, and only people who are writing their own Study
 * classes need look there.
 */

class SimTK_SimTKCOMMON_EXPORT Study {
public:
    class Guts; // local; name is Study::Guts
    friend class Guts;
private:
    // This is the only data member in this class. Also, any class derived from
    // Study must have *NO* data members at all (data goes in the Guts class).
    Guts* guts;
public:
    Study() : guts(0) { }
    Study(const Study&);
    Study& operator=(const Study&);
    ~Study();

    explicit Study(const System& sys);

    const String& getName()    const;
    const String& getVersion() const;

    const System& getSystem() const;
    const State&  getState()  const;
    State&        updState();

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;


    // There can be multiple handles on the same Study.
    bool isSameStudy(const Study& otherStudy) const;

    // Internal use only

    // dynamic_cast the returned reference to a reference to your concrete Guts
    // class.
    const Study::Guts& getStudyGuts() const {assert(guts); return *guts;}
    Study::Guts&       updStudyGuts()       {assert(guts); return *guts;}

    // Put new *unowned* Guts into this *empty* handle and take over ownership.
    // If this handle is already in use, or if Guts is already owned this
    // routine will throw an exception.
    void adoptStudyGuts(Study::Guts* g);

    explicit Study(Study::Guts* g) : guts(g) { }
    bool hasGuts() const {return guts!=0;}
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_STUDY_H_
