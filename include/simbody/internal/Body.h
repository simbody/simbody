#ifndef SimTK_SIMBODY_BODY_H_
#define SimTK_SIMBODY_BODY_H_

/* Portions copyright (c) 2007 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * This defines the Body class, which represents a reference frame which
 * can be used to describe mass properties and geometry. These are
 * in turn used to build MobilizedBodies which combine a body and a
 * specified mobilizer which defines how the reference frame can move
 * with respect to other MobilizedBodies.
 *
 * Body is an abstract base class, with concrete classes defined
 * for each kind of body. There are a set of built-in body types.
 * TODO: "Custom" bodies(?)
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

class SimTK_SIMBODY_EXPORT Body {
public:
    Body() : rep(0) { }
    ~Body();
    Body(const Body&);
    Body& operator=(const Body&);

    // Every type of Body should provide an initial set of rigid body mass 
    // properties defined at Topology stage (i.e., in the System rather than
    // the State). 
    const MassProperties& getDefaultRigidBodyMassProperties() const;

    // These are the built-in Body types.
    class Ground;     // infinitely massive
    class Massless;   // just a reference frame
    class Particle;   // point mass only; com=0; no inertia
    class Linear;     // point masses along a line; scalar inertia
    class Rigid;      // general rigid body
    class Deformable; // base class for bodies with internal deformation coords

    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    class BodyRep; // local subclass
    explicit Body(class BodyRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const BodyRep& getRep() const {assert(rep); return *rep;}
    BodyRep&       updRep() const {assert(rep); return *rep;}
	void setRep(BodyRep& r) {assert(!rep); rep = &r;}
protected:
    class BodyRep* rep;
};

class SimTK_SIMBODY_EXPORT Body::Rigid : public Body {
public:
    Rigid(); // default mass properties (1,Vec3(0),Inertia(1,1,1))
    Rigid(const MassProperties&);
    Rigid& setDefaultMassProperties(const MassProperties&);

    class RigidRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Rigid, Body);
private:
    RigidRep&       updRep();
    const RigidRep& getRep() const;
};

class SimTK_SIMBODY_EXPORT Body::Linear : public Body {
public:
    Linear(); // default mass properties (1,Vec3(0),Inertia(1,1,0))
    Linear(const MassProperties&);
    Linear& setDefaultMassProperties(const MassProperties&);

    class LinearRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Linear, Body);
private:
    LinearRep&       updRep();
    const LinearRep& getRep() const;
};

class SimTK_SIMBODY_EXPORT Body::Ground : public Body {
public:
    Ground();

    class GroundRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Ground, Body);
private:
    GroundRep&       updRep();
    const GroundRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_BODY_H_



