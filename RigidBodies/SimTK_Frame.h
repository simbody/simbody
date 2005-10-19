#ifndef SIMTK_FRAME_H_
#define SIMTK_FRAME_H_

/** @file
 * User-visible definitions for the objects that go into building a multibody system.
 * This is not the data structure used at run time, so the emphasis is on 
 * nice behavior for the caller. We'll have plenty of time for speed later.
 */

namespace simtk {

class Vec3;
class Mat33;
class MassProperties;
class Inertia;

class Frame {
public:
private:
    class FrameRep* rep;
};

class Station {
public:
    Station();
    ~Station();

    void setName(const String&);
    void setPlacement(const Frame&, const Vec3&);

private:
    class StationRep* rep;
};

class Direction {
public:
    Direction();
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

#endif /* SIMTK_FRAME_H_ */
