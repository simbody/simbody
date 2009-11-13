#ifndef SimTK_SIMBODY_BODY_H_
#define SimTK_SIMBODY_BODY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-9 Stanford University and the Authors.         *
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

/** @file
 * This defines the Body class, which represents a reference frame that
 * can be used to describe mass properties and geometry. These are
 * in turn used to build MobilizedBodies which combine a body and a
 * specified mobilizer that defines how the reference frame can move
 * with respect to other MobilizedBodies.
 *
 * Body is an abstract base class handle, with concrete classes defined
 * for each kind of body. There are a set of built-in body types.
 * TODO: "Custom" bodies(?)
 *
 * SimTK Design Patterns used:
 *    - abstract private implementation
 *    - binary compatible interface
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

class DecorativeGeometry;

class SimTK_SIMBODY_EXPORT Body {
public:
    Body() : rep(0) { }
    ~Body();
    Body(const Body&);
    Body& operator=(const Body&);

    /// This is a default conversion from MassProperties to Body. It will 
    /// result in a rigid body (concrete class Body::Rigid) being created 
    /// using these as its default MassProperties.
    Body(const MassProperties&);

    /// Add a piece of decorative geometry fixed at some location on this
    /// Body. This can be used for visualization of the Body's motion.
    /// Returns a reference to the Body so these can be chained like
    /// assignment operators.
    Body& addDecoration(const Transform& X_BD, const DecorativeGeometry&);

    // Every type of Body should provide an initial set of rigid body mass 
    // properties defined at Topology stage (i.e., in the System rather than
    // the State). 
    const MassProperties& getDefaultRigidBodyMassProperties() const;
    Body& setDefaultRigidBodyMassProperties(const MassProperties&);

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
    bool           hasRep() const {return rep!=0;}
    const BodyRep& getRep() const {assert(rep); return *rep;}
    BodyRep&       updRep() const {assert(rep); return *rep;}
	void           setRep(BodyRep& r) {assert(!rep); rep = &r;}
protected:
    class BodyRep* rep;
};

/**
 * A general rigid body. This can represent a body with mass properties
 * that are full, linear, inertialess (e.g. a point), or massless.
 */
class SimTK_SIMBODY_EXPORT Body::Rigid : public Body {
public:
    Rigid(); // default mass properties (1,Vec3(0),Inertia(1,1,1))
    explicit Rigid(const MassProperties&);

    Rigid& addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)Body::addDecoration(X_BD,g);
        return *this;
    }
    Rigid& setDefaultRigidBodyMassProperties(const MassProperties& m) {
        (void)Body::setDefaultRigidBodyMassProperties(m);
        return *this;
    }

    class RigidRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Rigid, Body);
private:
    RigidRep&       updRep();
    const RigidRep& getRep() const;
};

/**
 * A rigid body in the shape of a line, which is inherently inertialess
 * about its axis.
 */
class SimTK_SIMBODY_EXPORT Body::Linear : public Body {
public:
    Linear(); // default mass properties (1,Vec3(0),Inertia(1,1,0))
    explicit Linear(const MassProperties&);

    Linear& addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)Body::addDecoration(X_BD,g);
        return *this;
    }
    Linear& setDefaultRigidBodyMassProperties(const MassProperties& m) {
        (void)Body::setDefaultRigidBodyMassProperties(m);
        return *this;
    }

    class LinearRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Linear, Body);
private:
    LinearRep&       updRep();
    const LinearRep& getRep() const;
};

class SimTK_SIMBODY_EXPORT Body::Ground : public Body {
public:
    Ground();

    Ground& addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)Body::addDecoration(X_BD,g);
        return *this;
    }

    class GroundRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Ground, Body);
private:
    GroundRep&       updRep();
    const GroundRep& getRep() const;
};

class SimTK_SIMBODY_EXPORT Body::Massless : public Body {
public:
    Massless();

    Massless& addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)Body::addDecoration(X_BD,g);
        return *this;
    }

    class MasslessRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Massless, Body);
private:
    MasslessRep&       updRep();
    const MasslessRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_BODY_H_



