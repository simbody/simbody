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

#include <string>

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

class Vec3 {
    double d[3];
public:
    Vec3() { }
    Vec3(double xx, double yy, double zz) { d[0]=xx; d[1]=yy; d[2]=zz; }
};

class Mat33 {
    Vec3 cols[3];
public:
    Mat33() { }
    Mat33(Vec3 c1, Vec3 c2, Vec3 c3) { cols[0]=c1; cols[1]=c2; cols[2]=c3; }
};

class String : public std::string {
public:
    String() { }
    String(const char* c) : std::string(c) { }
};

class MassProperties;
class Inertia;

class Placement {
public:
    Placement();
    ~Placement();

    bool hasPlacement() const;
    const Frame& getParentFrame() const;
private:
    class PlacementRep* rep;
};

class Frame {
public:
    Frame();
    explicit Frame(const String&);
    ~Frame(); 

    const String getFullName() const;
private:
    class FrameRep* rep;
};

class Station {
public:
    Station();
    explicit Station(const String&);
    ~Station();

    void setName(const String&);
    void setPlacement(const Frame&, const Vec3&);

private:
    class StationRep* rep;
};

class Direction {
public:
    Direction();
    explicit Direction(const String&);
    ~Direction();
    
    void setName(const String&);
    void setPlacement(const Frame&, const Vec3&);
private:
    class DirectionRep* rep;
};

class MassElement {
public:
    MassElement();
    ~MassElement();
private:
    class MassElementRep* rep;
};

class Body : public Frame {
public:
    Body();
    ~Body();
private:
    class BodyRep* rep;
};

class RigidBody : public Body {
    RigidBody();
    ~RigidBody();
public:
    class RigidBodyRep* rep;
};

class DeformableBody : public Body {
    DeformableBody();
    ~DeformableBody();
public:
    class DeformableBodyRep* rep;
};

class Multibody : public Body {
    Multibody();
    ~Multibody();
public:
    class MultibodyRep* rep;
};

class Joint {
public:
    Joint();
    ~Joint();
private:
    class JointRep* rep;
};

class MultibodySystem {
    MultibodySystem();
    ~MultibodySystem();
public:
    class MultibodySystemRep* rep;
};

} // namespace simtk

#endif // SIMTK_MULTIBODY_MODELING_H_
