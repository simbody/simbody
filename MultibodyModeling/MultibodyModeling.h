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
 *
 * Feature: Station, Direction, Frame, MassElement, ...
 * Placement: constant, expression or feature
 * Body: is a Frame, has (Feature,Placement) pairs
 *
 */

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

#include <iostream>

namespace simtk {

class Placement;
class StationPlacement;
class DirectionPlacement;
class FramePlacement;

class Station;
class Direction;

class Frame;
class Body;
class RigidBody;
class DeformableBody;
class Multibody;

class Joint;
class MultibodySystem;
class MassProperties;
class Inertia;

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

// 
// <placement>     ::= constant | parameter | feature.placement
//                     | func(<placement> ...)
//
// Features provide a list of the placement types they need. Each of
// these list elements must be associated with a <placementExpr> eventually.
// Then these are evaluated at the appropriate runtime stage.
//
// <placement>'s have value types scalar, station, direction,
// orientation, or frame.
//                     

class Placement {
public:
    Placement() : rep(0) { }
    Placement(const Placement&);
    Placement& operator=(const Placement&);
    ~Placement();

    const Feature& getOwner() const;

    // Create a Placement which is evaluated by returning the
    // value of the indicated Feature's Placement. This is
    // deferred until runtime; the Feature may not even have
    // a placement yet.
    explicit Placement(const Feature&);

private:
    class PlacementRep* rep;
    friend class PlacementRep;
};
std::ostream& operator<<(std::ostream& o, const Placement&);

class ScalarPlacement : public Placement {
public:
    ScalarPlacement(const Real&);
    static ScalarPlacement plus  (const ScalarPlacement& l,
                                  const ScalarPlacement& r);
    static ScalarPlacement minus (const ScalarPlacement& l,
                                  const ScalarPlacement& r);
    static ScalarPlacement times (const ScalarPlacement& l,
                                  const ScalarPlacement& r);
    static ScalarPlacement divide(const ScalarPlacement& l,
                                  const ScalarPlacement& r);
    static ScalarPlacement length(const DirectionPlacement&);
private:
};

class StationPlacement : public Placement {
public:
    StationPlacement(const Vec3&);    // implicit conversion
    StationPlacement(const Station&); //   "
    StationPlacement(const Frame&);   //   "

    static StationPlacement   plus(const StationPlacement&,
                                   const DirectionPlacement&);
private:
};

class DirectionPlacement : public Placement {
public:
    DirectionPlacement(const Vec3&);      // implicit conversion
    DirectionPlacement(const Direction&); //   "

    static DirectionPlacement plus (const DirectionPlacement& l,
                                    const DirectionPlacement& r);
    static DirectionPlacement minus(const StationPlacement& head,
                                    const StationPlacement& tail);
    static DirectionPlacement minus(const DirectionPlacement& l,
                                    const DirectionPlacement& r);
private:
};

// Three, mutually orthogonal, right handed directions.
class OrientationPlacement : public Placement {
public:
    OrientationPlacement(); // identity
    OrientationPlacement(const Mat33&);      // implicit conversion
private:
};

class FramePlacement : public Placement {
public:
    FramePlacement(const Frame&);   // implicit conversions
    FramePlacement(const Station&); //   orientation inherited from 
                                    //   Station's owner
    FramePlacement(const OrientationPlacement&, 
                   const StationPlacement&);
private:
};

/**
 * Attachments are the "loose wires" dangling off features before they
 * mounted on anything. Attachments need to be "placed", which is done
 * by assigning a Placement to the Attachment. The Placement type must
 * be appropriate for this Attachment type or exceptions will fly.
 *
 * Attachments can be viewed as variables associated with a Feature, 
 * Placements are the expressions which are used to assign a value to
 * the Attachment variables. At run time the Placements are evaluated
 * to yield a numerical value.
 */
class Attachment {
public:
    Attachment() : rep(0) { }

    std::string    getName()  const {return std::string(getNamePtr());}
    int            getId()    const;
    const Feature& getOwner() const;

    bool isOKPlacement(const Placement&) const;
    void setPlacement(const Placement&);
    void setDefaultPlacement(const Placement&);
    bool hasDefaultPlacement() const;
    bool hasPlacement() const;
    const Placement& getPlacement() const;
    
private:
    // std::string class is not opaque so can't cross API reliably
    const char* getNamePtr() const;

    class AttachmentRep* rep;
    friend class AttachmentRep;
};
std::ostream& operator<<(std::ostream& o, const Attachment&);


/**
 * An Attachment that wants a scalar-valued placement.
 */
class ScalarAttachment : public Attachment {
public:
    ScalarAttachment();
    explicit ScalarAttachment(const char* name);
    void setPlacement(const ScalarPlacement&);
};

/**
 * An Attachment that wants a Station-valued placement.
 */
class StationAttachment : public Attachment {
public:
    StationAttachment();
    explicit StationAttachment(const char* name);
    void setPlacement(const StationPlacement&);
};

/**
 * An Attachment that wants a Direction-valued placement.
 */
class DirectionAttachment : public Attachment {
public:
    DirectionAttachment();
    explicit DirectionAttachment(const char* name);
    void setPlacement(const DirectionPlacement&);
};

/**
 * An Attachment that wants an Orientation-valued placement.
 */
class OrientationAttachment : public Attachment {
public:
    OrientationAttachment();
    explicit OrientationAttachment(const char* name);
    void setPlacement(const OrientationPlacement&);
};

/**
 * An Attachment that wants an Frame-valued placement.
 */
class FrameAttachment : public Attachment {
public:
    FrameAttachment();
    explicit FrameAttachment(const char* name);
    void setPlacement(const FramePlacement&);
};

/**
 * Features form a tree because many Features have child Features.
 * Parent features own their children (destructing the parent 
 * destructs all the children).
 * Most features require placement in order to be useful (e.g.,
 * a Station has to have a location). Placement expressions can
 * be constants or can involve parent or ancestor features, but
 * not children or siblings.
 * Every feature has a name and an index by which it is known
 * to its parent. Parentless features can exist but they can't
 * be placed. Parentless features still have a name but they do
 * not have an index.
 */
class Feature {
public:
    Feature() : rep(0) { }
    Feature(const Feature&);    // placements are not copied
    Feature& operator=(const Feature&);
    ~Feature();

    const String&  getName()  const; // name and index as known to parent
    int            getIndex() const; // -1 if no parent
    const Feature& getParentFeature() const;

    // Subfeatures of this feature
    int            getNFeatures() const;
    const Feature& getFeature(const String&) const;
    const Feature& getFeature(int) const;

    // true if Feature & all children have been placed
    bool hasPlacement() const;
    void place(const Placement&);
    const Placement& getPlacement() const;

private:
    class FeatureRep* rep;
    friend class FeatureRep;
};
std::ostream& operator<<(std::ostream& o, const Feature&);

class Parameter : public Feature {
public:
    explicit Parameter(const String& name);
    Parameter(const String& name, const Real& defaultValue);
    Parameter(const Parameter&);
    Parameter& operator=(const Parameter&);
    ~Parameter();

    void place(const ScalarPlacementExpr&);
};

class Station : public Feature {
public:
    explicit Station(const String& name);
    Station(const String& name, const Vec3& defaultValue);
    Station(const Station&);
    Station& operator=(const Station&);
    ~Station();

    void place(const StationPlacementExpr&);
};
std::ostream& operator<<(std::ostream& o, const Station&);

class Direction : public Feature {
public:
    explicit Direction(const String& name);
    Direction(const String& name, const Vec3& defaultValue);
    Direction(const Direction&);
    Direction& operator=(const Direction&);
    ~Direction();

    void place(const DirectionPlacementExpr&);
};
std::ostream& operator<<(std::ostream& o, const Direction&);

class Frame : public Feature {
public:
    explicit Frame(const String& name);
    Frame(const String& name);
    Frame(const Frame&);
    Frame& operator=(const Frame&);
    ~Frame();

    void place(const FramePlacementExpr&);
};
std::ostream& operator<<(std::ostream& o, const Frame&);


class Body : public Frame {
public:
    Body() : rep(0) { }
    explicit Body(const String&);
    Body(const Body&);
    Body& operator=(const Body&);
    ~Body();

    // Add features to this Frame
    Station&   addStation(const String&);
    Station&   addStation(const String&, const StationPlacement&);

    Direction& addDirection(const String&);
    Direction& addDirection(const String&, const DirectionPlacement&);

    Orientation& addOrientation(const String&);
    Orientation& addOrientation(const String&, const OrientationPlacement&);

    Frame& addFrame(const String&);
    Frame& addFrame(const String&, const FramePlacement&);
    Frame& addFrame(const String&, const OrientationPlacement&, 
                                   const StationPlacement&);

    const Station& getStation(int) const;
    const Station& getStation(const String&) const;

    const Direction& getDirection(int) const;
    const Direction& getDirection(const String&) const;

    const Orientation& getOrientation(int) const;
    const Orientation& getOrientation(const String&) const;

    const Frame& getFrame(int) const;
    const Frame& getFrame(const String&) const;

private:
    class BodyRep* rep;
    friend class BodyRep;
};
std::ostream& operator<<(std::ostream& o, const Body&);

class RigidBody : public Body {
public:
    RigidBody() : rep(0) { }
    explicit RigidBody(const String&);
    RigidBody(const RigidBody&);
    RigidBody& operator=(const RigidBody&);
    ~RigidBody();

    // add RigidBody features
    MassElement& addMassElementLike(const MassElement&, const Placement&,
                                    const char* name=0);
private:
    class RigidBodyRep* rep;
    friend class RigidBodyRep;
};
std::ostream& operator<<(std::ostream& o, const RigidBody&);

class DeformableBody : public Body {
public:
    DeformableBody() : rep(0) { }
    explicit DeformableBody(const String&);
    DeformableBody(const DeformableBody&);
    DeformableBody& operator=(const DeformableBody&);
    ~DeformableBody();
private:
    class DeformableBodyRep* rep;
    friend class DeformableBodyRep;
};
std::ostream& operator<<(std::ostream& o, const DeformableBody&);

class Multibody : public Body {
public:
    Multibody() : rep(0) { }
    explicit Multibody(const String&);
    Multibody(const Multibody&);
    Multibody& operator=(const Multibody&);
    ~Multibody();

    const Frame& getGroundFrame() const;

    // add Multibody features
    Body& addGroundBody(const char* name=0);
    Body& addBody(const char* name=0);
    Body& addBodyLike(const Body&, const char* name=0);
    Joint& addJoint(const Placement& reference,
                    const Placement& moving,
                    int jointType,
                    const char* name=0);
private:
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
