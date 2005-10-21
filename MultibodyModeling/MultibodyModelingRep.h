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
#include <sstream>

namespace simtk {

// Declare the 'rep' classes which actually represent the objects pointed
// to by the wrapper classes:
//    Frame,Station,Direction,MassElement,Body,Joint,RigidBody,DeformableBody,
//    Multibody,MultibodySystem

class PlacementRep {
public:
    PlacementRep(const Frame& f) : parentFrame(&f) { }
    virtual ~PlacementRep() { }

    bool hasPlacement() const { return parentFrame != 0; }
    const Frame& getParentFrame() const {assert(parentFrame); return *parentFrame;}

    virtual PlacementRep* clone() const = 0;
    virtual std::string   getPlacementTypeName() const = 0;
    virtual std::string   toString() const = 0;

    // Placement helpers
    static const PlacementRep* getRep(const Placement& p)   { return p.rep; }
    static PlacementRep*       updRep(Placement& p)         { return p.rep; }
    static void setRep(Placement& p, PlacementRep* rep)     
    { assert(p.rep==0); p.rep=rep; }
    static void replaceRep(Placement& p, PlacementRep* rep) { delete p.rep; p.rep=rep; }
    static void clearRep(Placement& p) { replaceRep(p,0); }
private:
    const Frame*  parentFrame;    // reference to an existing Frame, don't delete!
};

// TODO should be parameterizable.
class StationPlacementRep : public PlacementRep {
public:
    StationPlacementRep(const Frame& f) : PlacementRep(f) { }
    virtual ~StationPlacementRep() { }

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers()     const = 0;
};
// A concrete StationPlacement in which there are no variables.
class FixedStationPlacementRep : public StationPlacementRep {
public:
    FixedStationPlacementRep(const Frame& f, const Vec3& v)
      : StationPlacementRep(f), station(v) { }

    // Implementations of pure virtuals.
    Vec3 getMeasureNumbers() const { return station; }
    PlacementRep* clone() const { return new FixedStationPlacementRep(*this); }
    std::string getPlacementTypeName() const { return "Constant Station"; }
    std::string toString() const {
        std::stringstream s;
        s << getPlacementTypeName() << " Placement: " << station;
        return s.str();
    }
private:
    Vec3 station;
};

class DirectionPlacementRep : public PlacementRep {
public:
    DirectionPlacementRep(const Frame& f) : PlacementRep(f) { }
    virtual ~DirectionPlacementRep() { }

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers()     const = 0;
};

// A concrete StationPlacement in which there are no variables.
class FixedDirectionPlacementRep : public DirectionPlacementRep {
public:
    FixedDirectionPlacementRep(const Frame& f, const Vec3& v)
      : DirectionPlacementRep(f), direction(v) { }

    // Implementations of pure virtuals.
    Vec3 getMeasureNumbers() const { return direction; }
    PlacementRep* clone() const { return new FixedDirectionPlacementRep(*this); }
    std::string getPlacementTypeName() const { return "Constant Direction"; }
    std::string toString() const {
        std::stringstream s;
        s << getPlacementTypeName() << " Placement: " << direction;
        return s.str();
    }
private:
    Vec3 direction;
};

class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep(const Frame& f) : PlacementRep(f) { }
    virtual ~FramePlacementRep() { }

    // These should allow for state to be passed in.
    virtual Vec3  getStation()     const = 0;   // origin measure numbers
    virtual Mat33 getOrientation() const = 0;   // columns are x,y,z measure numbers
};

// A concrete Placement in which there are no variables.
class FixedFramePlacementRep : public FramePlacementRep {
public:
    FixedFramePlacementRep(const Frame& f, const Vec3& sta, const Mat33& ori)
      : FramePlacementRep(f), station(sta), orientation(ori) { }

    // Implementations of pure virtuals.
    Vec3  getStation()     const { return station; }
    Mat33 getOrientation() const { return orientation; }
    PlacementRep* clone() const { return new FixedFramePlacementRep(*this); }
    std::string getPlacementTypeName() const { return "Constant Frame"; }
    std::string toString() const {
        std::stringstream s;
        s << getPlacementTypeName() << " Placement: station=" 
          << station << " orientation=" << std::endl;
        s << orientation;
        return s.str();
    }
private:
    Vec3    station;
    Mat33   orientation;
};

class FrameRep {
public:
    FrameRep(const Frame& f) : name("frame") {
        initializeFeatures(f);
    }
    FrameRep(const Frame& f, const char* nm) : name(nm) {
        initializeFeatures(f);
    }
    const Placement& getPlacement() const { return place; }
    void setPlacement(const Frame& f, const Vec3& org, const Mat33& ori) {
        PlacementRep::setRep(place, new FixedFramePlacementRep(f,org,ori));
    }
    void removePlacement() {
        PlacementRep::clearRep(place);
    }

    const std::vector<Station>& getStations() const { return childStations; }
    const std::vector<Direction>& getDirections() const { return childDirections; }
    const std::vector<Frame>& getFrames() const { return childFrames; }

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

    // Frame helpers
    static const FrameRep* getRep(const Frame& f)   { return f.rep; }
    static FrameRep*       updRep(Frame& f)         { return f.rep; }
    static void setRep(Frame& f, FrameRep* rep)     { assert(f.rep==0); f.rep=rep; }
    static void replaceRep(Frame& f, FrameRep* rep) { delete f.rep; f.rep=rep; }
    static void clearRep(Frame& f) { replaceRep(f,0); }

private:
    void initializeFeatures(const Frame& f) {
        childStations.clear();
        childDirections.clear();
        childFrames.clear();

        childStations.push_back  (Station("O"));
        childDirections.push_back(Direction("x"));
        childDirections.push_back(Direction("y"));
        childDirections.push_back(Direction("z"));

        childStations[0].setPlacement  (f, Vec3(0,0,0));
        childDirections[0].setPlacement(f, Vec3(1,0,0));
        childDirections[1].setPlacement(f, Vec3(0,1,0));
        childDirections[2].setPlacement(f, Vec3(0,0,1));
    }

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
    StationRep(const Station&) : name("station") { }
    StationRep(const Station&, const char* nm) : name(nm) { }
    const Placement& getPlacement() const { return station; }
    void setPlacement(const Frame& f, const Vec3& v) {
        PlacementRep::setRep(station, new FixedStationPlacementRep(f,v));
    }
    void removePlacement() {
        PlacementRep::clearRep(station);
    }
    const std::string& getName() const { return name; }

    
    // Station helpers
    static const StationRep* getRep(const Station& s)   { return s.rep; }
    static StationRep*       updRep(Station& s)         { return s.rep; }
    static void setRep(Station& s, StationRep* rep)     { assert(s.rep==0); s.rep=rep; }
    static void replaceRep(Station& s, StationRep* rep) { delete s.rep; s.rep=rep; }
    static void clearRep(Station& s) { replaceRep(s,0); }
private:
    std::string name;
    Placement   station;
};

class DirectionRep {
public:
    DirectionRep(const Direction&) : name("direction") { }
    DirectionRep(const Direction&, const char* nm) : name(nm) { }
    const Placement& getPlacement() const { return direction; }
    void setPlacement(const Frame& f, const Vec3& v) {
        PlacementRep::setRep(direction, new FixedDirectionPlacementRep(f,v));
    }
    void removePlacement() {
        PlacementRep::clearRep(direction);
    }
    const std::string& getName() const { return name; }

    // Direction helpers
    static const DirectionRep* getRep(const Direction& d)   { return d.rep; }
    static DirectionRep*       updRep(Direction& d)         { return d.rep; }
    static void setRep(Direction& d, DirectionRep* rep)     { assert(d.rep==0); d.rep=rep; }
    static void replaceRep(Direction& d, DirectionRep* rep) { delete d.rep; d.rep=rep; }
    static void clearRep(Direction& d) { replaceRep(d,0); }
private:
    std::string name;
    Placement   direction;
};

// Abstract base class for MassElements.
class MassElementRep {
public:
    virtual Real  getMass()         const = 0;
    virtual Vec3  getCenterOfMass() const = 0;
    virtual Mat33 getInertia()      const = 0;
    virtual MassElementRep* clone() const = 0;

    const std::string& getName() const { return name; }
    void setName(const char* nm) { name=nm; }
    const MassElement& getHandle() const { return handle; }
    MassElement&       updHandle()       { return handle; }

    const Placement& getPlacement() const { return where; }
    void setPlacement(const Placement& pl) {
        PlacementRep::setRep(where, PlacementRep::getRep(pl)->clone());
    }
    void removePlacement() {
        PlacementRep::clearRep(where);
    }

    // MassElement helpers
    static const MassElementRep* getRep(const MassElement& me)   { return me.rep; }
    static MassElementRep*       updRep(MassElement& me)         { return me.rep; }
    static void setRep(MassElement& me, MassElementRep* rep)     { assert(me.rep==0); me.rep=rep; }
    static void replaceRep(MassElement& me, MassElementRep* rep) { delete me.rep; me.rep=rep; }
    static void clearRep(MassElement& me) { replaceRep(me,0); }
protected:
    MassElementRep(MassElement& me) : handle(me) { }

    MassElement&  handle;
    std::string   name;
    Placement     where;
};

class PointMassElementRep : public MassElementRep {
public:
    PointMassElementRep(PointMassElement& pme)
      : MassElementRep(pme), mass(-1) { }

    void setMass(const Real& m) { mass=m; }

    // Override for convenience.
    const PointMassElement& getHandle() const {
        return reinterpret_cast<const PointMassElement&>(handle);
    }
    PointMassElement& updHandle() {
        return reinterpret_cast<PointMassElement&>(handle);
    }

    // Implementations of virtuals.
    Real  getMass()         const { assert(mass >= 0.); return mass; }
    Vec3  getCenterOfMass() const { return Vec3(0,0,0); }
    Mat33 getInertia()      const {
        const Vec3 z(0,0,0);
        return Mat33(z,z,z);
    }
    MassElementRep* clone() const { return new PointMassElementRep(*this); }
private:
    Real mass;
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
