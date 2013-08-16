/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
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
    }
    GutsRep(const String& name, const String& version) 
      : studyName(name), studyVersion(version), myHandle(0)
    {
    }

    GutsRep(const GutsRep& src)
    :   studyName(src.studyName),
        studyVersion(src.studyVersion),
        myHandle(0)
    {
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

Study& Study::operator=(const Study& src) {
    if (!isSameStudy(src)) {
        if (isOwnerHandle())
            delete guts;
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerHandle(*this);
        }
    }
    return *this;
}

Study::~Study() {
    if (guts && isOwnerHandle())
        delete guts;
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

// This is also the default constructor.
Study::Guts::Guts(const String& name, const String& version) {
    rep = new GutsRep(name,version);
    // note that the GutsRep object currently has no owner handle
}

Study::Guts::~Guts() {
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

Study::Guts* Study::Guts::clone() const {
    return cloneImpl();
}


    //////////////////////////
    // STUDY::GUTS::GUTSREP //
    //////////////////////////

// All inline currently.


} // namespace SimTK

