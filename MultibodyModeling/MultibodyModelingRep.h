#ifndef SIMTK_MULTIBODY_MODELING_REP_H_
#define SIMTK_MULTIBODY_MODELING_REP_H_

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

/**@file
 * Declarations for the *real* Multibody Modeling objects. These are opaque to
 * users.
 */

#include "MultibodyModeling.h"

#include <string>
#include <vector>
#include <cassert>

namespace simtk {

// Declare the 'rep' classes which actually represent the objects pointed
// to by the wrapper classes:
//    Frame,Station,Direction,MassElement,Body,Joint,RigidBody,DeformableBody,
//    Multibody,MultibodySystem

class PlacementRep {
public:
    PlacementRep() : parentFrame(0) { }
    virtual ~PlacementRep() { }

    bool hasPlacement() const { return parentFrame != 0; }
    const Frame& getParentFrame() const {assert(parentFrame); return *parentFrame;}
private:
    Frame*  parentFrame;    // reference to an existing Frame, don't delete!
};

// TODO should be parameterizable.
class StationPlacementRep : public PlacementRep {
public:
    StationPlacementRep() { }
    virtual ~StationPlacementRep() { }

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers()     const = 0;
};

class DirectionPlacementRep : public PlacementRep {
public:
    DirectionPlacementRep() { }
    virtual ~DirectionPlacementRep() { }

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers()     const = 0;
};

class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep() { }
    virtual ~FramePlacementRep() { }

    // These should allow for state to be passed in.
    virtual Vec3  getStation()     const = 0;   // origin measure numbers
    virtual Mat33 getOrientation() const = 0;   // columns are x,y,z measure numbers
};

// A concrete Placement in which there are no variables.
class FixedFramePlacementRep : public FramePlacementRep {
public:
    FixedFramePlacementRep(const Vec3& sta, const Mat33& ori)
      : station(sta), orientation(ori)
    { }

    // Implementations of pure virtuals.
    Vec3  getStation()     const { return station; }
    Mat33 getOrientation() const { return orientation; }

private:
    Vec3    station;
    Mat33   orientation;
};

class FrameRep {
public:
    FrameRep(const Frame& f, const char* nm) {
        name = nm;
        childStations.push_back  (Station("O"));
        childDirections.push_back(Direction("x"));
        childDirections.push_back(Direction("y"));
        childDirections.push_back(Direction("z"));

        childStations[0].setPlacement  (f, Vec3(0,0,0));
        childDirections[0].setPlacement(f, Vec3(1,0,0));
        childDirections[1].setPlacement(f, Vec3(0,1,0));
        childDirections[2].setPlacement(f, Vec3(0,0,1));
    }

    const std::string getFullName() const { 
        std::string s;
        if (place.hasPlacement())
            s = place.getParentFrame().getFullName() + "/";
        return s + getLocalName(); 
    }

    const std::string getLocalName() const { return name; }

    const Station& getOrigin() const { return childStations[0]; }
    const Direction& getAxis(int i) const 
      { assert(0<=i && i<=2); return childDirections[i]; }
    const Direction& x() const { return childDirections[0]; }
    const Direction& y() const { return childDirections[1]; }
    const Direction& z() const { return childDirections[2]; }

    const Station&   addStation  (const Placement&, const char* nm);
    const Direction& addDirection(const Placement&, const char* nm);
    const Frame&     addFrame    (const Placement&, const char* nm);

private:
    std::string       name;

    // Objects wholly owned by this Frame and placed on it. The first
    // station is this Frame's origin, and the first three Directions
    // are this Frame's x,y,z axes, respectively. The rest are in no
    // particular order.
    std::vector<Station>   childStations;
    std::vector<Direction> childDirections;
    std::vector<Frame>     childFrames;

    // If this Frame has been placed with respect to some other
    // Frame, this is the placement information.
    Placement place;
};

class StationRep {
public:
    StationRep(const Station&, const char* nm) : name(nm) { }
private:
    std::string name;
    Placement   station;
};

class DirectionRep {
public:
    DirectionRep(const Direction&, const char* nm) : name(nm) { }
private:
    std::string name;
    Placement   direction;
};

class MassElementRep {
public:
    MassElementRep(const MassElement&, const char* nm) : name(nm) { }
private:
    std::string name;
};

class JointRep {
public:
    JointRep(const Joint&, const char* nm) : name(nm) { }
private:
    std::string name;
};

class MultibodySystemRep {
public:
    MultibodySystemRep(const MultibodySystem&, const char* nm) : name(nm) { }
private:
    std::string name;
};

// Body extends Frame.
class BodyRep {
public:
    BodyRep(const Body&) { }
    void addMassElementLike(const MassElement&, const Placement&);
private:
    std::vector<MassElement> childMassElements;
};

// RigidBody extends Body.
class RigidBodyRep {
public:
    RigidBodyRep(const RigidBody&) { }
private:
};

// DeformableBody extends Body.
class DeformableBodyRep {
public:
    DeformableBodyRep(const DeformableBody&) { }
private:
    std::string name;
};

// Multibody extends Body.
class MultibodyRep {
public:
    MultibodyRep(const Multibody&) { }
private:
    std::string name;
};

} // namespace simtk

#endif // SIMTK_MULTIBODY_MODELING_REP_H_
