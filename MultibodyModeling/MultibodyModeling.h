#ifndef SIMTK_MULTIBODY_MODELING_H_
#define SIMTK_MULTIBODY_MODELING_H_

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
 */

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

#include <iostream>

namespace simtk {

class Station;
class Direction;
class Frame;
class Body;
class Joint;
class RigidBody;
class DeformableBody;
class Multibody;
class MultibodySystem;
class MassProperties;
class Inertia;

class Placement {
public:
    Placement() : rep(0) { }
    Placement(const Placement&);
    Placement& operator=(const Placement&);
    ~Placement();

    bool hasPlacement() const;
    const Frame& getParentFrame() const;
private:
    class PlacementRep* rep;
    friend class PlacementRep;
};
std::ostream& operator<<(std::ostream& o, const Placement&);

class Frame {
public:
    Frame() : rep(0) { }
    explicit Frame(const String&);
    Frame(const Frame&);    // placements are not copied
    Frame& operator=(const Frame&);
    ~Frame(); 

    String getFullName() const;
private:
    class FrameRep* rep;
    friend class FrameRep;
};
std::ostream& operator<<(std::ostream& o, const Frame&);

class Station {
public:
    Station() : rep(0) { }
    explicit Station(const String&);
    Station(const Station&);    // placements are not copied
    Station& operator=(const Station&);
    ~Station();

    Vec3 getMeasureNumbers() const;

    void setName(const String&);
    void setPlacement(const Frame&, const Vec3&);

private:
    class StationRep* rep;
    friend class StationRep;
};
std::ostream& operator<<(std::ostream& o, const Station&);

class Direction {
public:
    Direction() : rep(0) { }
    explicit Direction(const String&);
    Direction(const Direction&);    // placements are not copied
    Direction& operator=(const Direction&);
    ~Direction();

    Vec3 getMeasureNumbers() const;
    
    void setName(const String&);
    void setPlacement(const Frame&, const Vec3&);
private:
    class DirectionRep* rep;
    friend class DirectionRep;
};
std::ostream& operator<<(std::ostream& o, const Direction&);

// This is an abstract handle class.
class MassElement {
public:
    MassElement() : rep(0) { }
    MassElement(const MassElement&);    // placements are not copied
    MassElement& operator=(const MassElement&);
    ~MassElement();
private:
    class MassElementRep* rep;
    friend class MassElementRep;
};
std::ostream& operator<<(std::ostream& o, const MassElement&);

class PointMassElement : public MassElement {
public:
    PointMassElement();
    explicit PointMassElement(const String&);
    explicit PointMassElement(const Real&);  // default mass
    PointMassElement(const Real&, const String&);   
};

class Body : public Frame {
public:
    Body() : rep(0) { }
    explicit Body(const String&);
    Body(const Body&);
    Body& operator=(const Body&);
    ~Body();
private:
    class BodyRep* rep;
    friend class BodyRep;
};
std::ostream& operator<<(std::ostream& o, const Body&);

class RigidBody : public Body {
    RigidBody() : rep(0) { }
    explicit RigidBody(const String&);
    RigidBody(const RigidBody&);
    RigidBody& operator=(const RigidBody&);
    ~RigidBody();
public:
    class RigidBodyRep* rep;
    friend class RigidBodyRep;
};
std::ostream& operator<<(std::ostream& o, const RigidBody&);

class DeformableBody : public Body {
    DeformableBody() : rep(0) { }
    explicit DeformableBody(const String&);
    DeformableBody(const DeformableBody&);
    DeformableBody& operator=(const DeformableBody&);
    ~DeformableBody();
public:
    class DeformableBodyRep* rep;
    friend class DeformableBodyRep;
};
std::ostream& operator<<(std::ostream& o, const DeformableBody&);

class Multibody : public Body {
    Multibody() : rep(0) { }
    explicit Multibody(const String&);
    Multibody(const Multibody&);
    Multibody& operator=(const Multibody&);
    ~Multibody();
public:
    class MultibodyRep* rep;
    friend class MultibodyRep;
};
std::ostream& operator<<(std::ostream& o, const Multibody&);

class Joint {
public:
    Joint() : rep(0) { }
    explicit Joint(const String&);
    Joint(const Joint&);
    Joint& operator=(const Joint&);
    ~Joint();
private:
    class JointRep* rep;
    friend class JointRep;
};
std::ostream& operator<<(std::ostream& o, const Joint&);

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

#endif // SIMTK_MULTIBODY_MODELING_H_
