#ifndef SimTK_SimTKCOMMON_STUDY_GUTS_H_
#define SimTK_SimTKCOMMON_STUDY_GUTS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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

namespace SimTK {

// TODO: more to come.

/**
 * This is the declaration for the Study::Guts class, the abstract object to
 * which a Study handle points. This is in a separate header file from Study
 * because only people who are extending the Study class to make their own
 * Studies need to be aware of the details. End users access only methods from
 * the Study class and classes derived from Study, never anything from
 * Study::Guts or its derived classes.
 *
 * Below is the physical layout of memory for a Study, and which
 * portions are allocated by the client program and which by the
 * binary library code. For binary compatiblity, only the side
 * which allocated a piece of memory can access it. Exception: both
 * the client and library side must agree on the virtual function
 * table (VFT) ordering of the client's virtual functions.
 *
 *               CLIENT SIDE                    .  LIBRARY SIDE
 *                                              .
 *        Study               Study::Guts       . Study::Guts::GutsRep
 *   ---------------       ------------------   .   -------------
 *  | Study::Guts*  | --> | Study::GutsRep*  | --> |   GutsRep   |
 *   ---------------       ------------------   .  |             |
 *          ^             | Concrete Guts    |  .  | Other opaque|
 *          |             | class data and   |  .  |   stuff     |
 *   ===============      | client-side VFT  |  .  |             |
 *   Concrete Study        ------------------   .  |             |
 *    adds no data                              .   -------------
 *       members   
 *
 * If the concrete Study::Guts class also has an opaque implementation,
 * as it will for concrete Studies provided by the SimTK Core, then
 * the Study author should expose only the data-free handle class 
 * derived from Study.
 */
class SimTK_SimTKCOMMON_EXPORT Study::Guts {
    class GutsRep;
    friend class GutsRep;

    // This is the only data member in this class.
    GutsRep* rep; // opaque implementation of Study::Guts base class.
public:
    // Note that this serves as a default constructor since both arguments have defaults.
    explicit Guts(const String& name="<UNNAMED STUDY>", 
                  const String& version="0.0.0");
    virtual ~Guts();

    const String& getName()    const;
    const String& getVersion() const;

    // Obtain the owner handle for this Study::Guts object.
    const Study& getStudy() const;
    Study& updStudy();

    void setOwnerHandle(Study&);
    bool hasOwnerHandle() const;

    explicit Guts(class GutsRep* r) : rep(r) { }
    bool           hasRep() const {return rep!=0;}
    const GutsRep& getRep() const {assert(rep); return *rep;}
    GutsRep&       updRep() const {assert(rep); return *rep;}

    // Wrap the cloneImpl virtual method.
    Study::Guts* clone() const;

protected:
    Guts(const Guts&);  // copies the base class; for use from derived class copy constructors
    
    // The destructor is already virtual; see above.

    virtual Study::Guts* cloneImpl() const = 0;

private:
    Guts& operator=(const Guts&); // suppress default copy assignment operator
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_STUDY_GUTS_H_
