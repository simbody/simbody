#ifndef SIMTK_BODY_H_
#define SIMTK_BODY_H_

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

/** @file
 * User-visible definitions for the objects that go into building a multibody system.
 * This is not the data structure used at run time, so the emphasis is on 
 * nice behavior for the caller. We'll have plenty of time for speed later.
 *
 * Feature: Station, Direction, Frame, MassElement, ...
 * Placement: constant, expression or feature
 * Body: is a Frame, has (Feature,Placement) pairs
 *
 */

#include "Feature.h"
#include "MassElement.h"
#include "Placement.h"

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"
#include <iostream>

namespace simtk {

// Declared below.
class Body;
class RigidBody;
class DeformableBody;
class Multibody;
class Joint;


enum JointType {
    UnknownJointType    = 0,
    ThisIsGround        = 1, // Ground's "inboard joint"
    WeldJoint           = 2,
    TorsionJoint        = 3,
    PinJoint            = TorsionJoint,
    SlidingJoint        = 4,
    UJoint              = 5,
    CylinderJoint       = 6,
    PlanarJoint         = 7,
    GimbalJoint         = 8,
    OrientationJoint    = 9,
    BallJoint           = OrientationJoint,
    CartesianJoint      = 10,
    FreeLineJoint       = 11,
    FreeJoint           = 12
};

/**
 * This is an abstract class representing all the bodies.
 */
class Body : public Frame {
public:
    const RealMeasure&    getMassMeasure() const;
    const StationMeasure& getCentroidMeasure() const;

    static bool        isInstanceOf(const Feature&);
    static const Body& downcast(const Feature&);
    static Body&       downcast(Feature&);
};

class RigidBody : public Body {
public:
    explicit RigidBody(const String&);
    RigidBody(const RigidBody&);
    RigidBody& operator=(const RigidBody&);
    ~RigidBody();

    // add RigidBody features
    MassElement& addMassElementLike(const MassElement&, const String&,
                                    const Placement&);
    MassElement& addMassElementLike(const MassElement&, const String&);

    static bool             isInstanceOf(const Feature&);
    static const RigidBody& downcast(const Feature&);
    static RigidBody&       downcast(Feature&);
};

class DeformableBody : public Body {
public:
    explicit DeformableBody(const String&);
    DeformableBody(const DeformableBody&);
    DeformableBody& operator=(const DeformableBody&);
    ~DeformableBody();

    static bool                  isInstanceOf(const Feature&);
    static const DeformableBody& downcast(const Feature&);
    static DeformableBody&       downcast(Feature&);
};

class Multibody : public Body {
public:
    explicit Multibody(const String&);
    Multibody(const Multibody&);
    Multibody& operator=(const Multibody&);
    ~Multibody();

    const Frame& getGroundFrame() const;

    // add Multibody features
    RigidBody& addGroundBody();
    RigidBody& addRigidBody(const String&);
    RigidBody& addRigidBodyLike(const RigidBody&, const String&);
    Body& addBodyLike(const Body&, const String&);
    Joint& addJoint(JointType, const String&);
    Joint& addJoint(JointType, const String&,
                    const FramePlacement& reference,
                    const FramePlacement& moving);

    static bool             isInstanceOf(const Feature&);
    static const Multibody& downcast(const Feature&);
    static Multibody&       downcast(Feature&);
};

class Joint : public Feature {
public:
    explicit Joint(JointType, const String&);
    Joint(const Joint&);
    Joint& operator=(const Joint&);
    ~Joint();

    static bool         isInstanceOf(const Feature&);
    static const Joint& downcast(const Feature&);
    static Joint&       downcast(Feature&);
};

class MultibodySystem {
    MultibodySystem() : rep(0) { }
    explicit MultibodySystem(const String&);
    MultibodySystem(const MultibodySystem&);
    MultibodySystem& operator=(const MultibodySystem&);
    ~MultibodySystem();
public:
    class MultibodySystemRep* rep;
    friend class MultibodySystemRep;
};
std::ostream& operator<<(std::ostream& o, const MultibodySystem&);

} // namespace simtk

#endif // SIMTK_BODY_H_
