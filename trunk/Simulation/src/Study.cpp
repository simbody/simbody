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


/**@file
 *
 * Implementation of Study and Study::Guts, and declaration and implementation
 * of Study::GutsRep.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Study.h"
#include "SimTKcommon/internal/StudyGuts.h"

#include <cassert>

namespace SimTK {


// This is the hidden implementation of the Study::Guts abstract base class.
class Study::Guts::GutsRep {
public:
    GutsRep() 
      : studyName("<NONAME>"), studyVersion("0.0.0"), myHandle(0)
    {
        clearAllFunctionPointers();
    }
    GutsRep(const String& name, const String& version) 
      : studyName(name), studyVersion(version), myHandle(0)
    {
    }

    GutsRep(const GutsRep& src) {
        studyName = src.studyName;
        studyVersion = src.studyVersion;
        myHandle = 0;
        copyAllFunctionPointers(src);
    }

    ~GutsRep() {
        clearMyHandle();
    }

    const String& getName()    const {return studyName;}
    const String& getVersion() const {return studyVersion;}

    void setMyHandle(Study& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

protected:
    String studyName;
    String studyVersion;

private:
    friend class Study;
    friend class Study::Guts;
    Study* myHandle;     // the owner handle of these guts

        // POINTERS TO CLIENT-SIDE FUNCTION LOCATORS

        // This is a virtual function table, but the addresses are
        // determined at run time so that we don't have to depend on a
        // particular ordering in the client side virtual function table.

    Study::Guts::DestructImplLocator destructp;
    Study::Guts::CloneImplLocator    clonep;

    void clearAllFunctionPointers() {
        destructp = 0;
        clonep    = 0;
    }

    void copyAllFunctionPointers(const GutsRep& src) {
        destructp = src.destructp;
        clonep    = src.clonep;
    }
};

    ///////////
    // STUDY //
    ///////////

bool Study::isEmptyHandle() const {return guts==0;}
bool Study::isOwnerHandle() const {return guts==0 || &guts->getStudy()==this;}
bool Study::isSameStudy(const Study& otherStudy) const {
    return guts && (guts==otherStudy.guts);
}

Study::Study(const Study& src) : guts(0) {
    if (src.guts) {
        guts = src.guts->clone();
        guts->setOwnerHandle(*this);
    }
}

// Don't use ordinary delete, assignment, or copy here. Must go
// through the library-side VFT to get access to the correct client-side
// virtual destructor and clone method.
Study& Study::operator=(const Study& src) {
    if (!isSameStudy(src)) {
        if (isOwnerHandle())
            Study::Guts::destruct(guts);
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerHandle(*this);
        }
    }
    return *this;
}

Study::~Study() {
    // Must delete using the library-side VFT, so that we can get access
    // to the client side virtual destructor to destruct this client-side
    // Study::Guts object.
    if (guts && isOwnerHandle())
        Study::Guts::destruct(guts);
    guts=0;
}

void Study::adoptStudyGuts(Study::Guts* g) {
    SimTK_ASSERT_ALWAYS(g, "Study::adoptStudyGuts(): can't adopt null Guts");
    SimTK_ASSERT_ALWAYS(!guts,
        "Study::adoptStudyGuts(): this Study handle is already in use");
    guts = g;
    guts->setOwnerHandle(*this);
}

const String& Study::getName()    const {return getStudyGuts().getName();}
const String& Study::getVersion() const {return getStudyGuts().getVersion();}

    /////////////////
    // STUDY::GUTS //
    /////////////////

// Default constructor is inline, but calls librarySideConstuction() here.
void Study::Guts::librarySideConstruction(const String& name, const String& version) {
    rep = new GutsRep(name,version);
    // note that the GutsRep object currently has no owner handle
}

// Destructor is inline, but calls librarySideDestruction() here.
void Study::Guts::librarySideDestruction() {
    delete rep;
    rep=0;
}

// Copy constructor
Study::Guts::Guts(const Guts& src) : rep(0) {
    if (src.rep) {
        rep = new GutsRep(*src.rep);
        // note that the GutsRep object currently has no owner handle
    }
}

// Copy assignment is suppressed
    

const Study& Study::Guts::getStudy() const {
    assert(rep->myHandle);
    return *rep->myHandle;
}

Study& Study::Guts::updStudy() {
    assert(rep->myHandle);
    return *rep->myHandle;
}

void Study::Guts::setOwnerHandle(Study& sys) {
    assert(!rep->myHandle);
    rep->myHandle = &sys;
}

bool Study::Guts::hasOwnerHandle() const {
    return rep->myHandle != 0;
}

const String& Study::Guts::getName()    const {return getRep().getName();}
const String& Study::Guts::getVersion() const {return getRep().getVersion();}

void Study::Guts::registerDestructImpl(DestructImplLocator f) {
    updRep().destructp = f;
}
void Study::Guts::registerCloneImpl(CloneImplLocator f) {
    updRep().clonep = f;
}

Study::Guts* Study::Guts::clone() const {
    return getRep().clonep(*this);
}


/*static*/void Study::Guts::destruct(Study::Guts* gutsp) {
    if (gutsp)
        gutsp->getRep().destructp(gutsp);
}

    //////////////////////////
    // STUDY::GUTS::GUTSREP //
    //////////////////////////

// All inline currently.


} // namespace SimTK

