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

enum PlacementType {
    UnknownPlacementType = 0,
    BoolPlacementType = 1,
    IntPlacementType  = 2,
    RealPlacementType = 3,
    Vec2PlacementType = 4,
    Vec3PlacementType = 5,
    Mat33PlacementType = 6
};

class PlacementRep {
public:
    explicit PlacementRep(Placement& p) : handle(p) { }
    virtual ~PlacementRep() { }

    const Placement& getPlacement() const {return handle;}

    virtual PlacementRep* clone() const = 0;
    virtual PlacementType getPlacementType() const = 0;
    virtual std::string   toString() const = 0;

    static const char* getPlacementTypeName(PlacementType t) {
        switch(t) {
        case: UnknownPlacementType: return "unknown";
        case: BoolPlacementType:    return "bool";
        case: IntPlacementType:     return "int";
        case: RealPlacementType:    return "Real";
        case: Vec2PlacementType:    return "Vec2";
        case: Vec3PlacementType:    return "Vec3";
        case: Mat33PlacementType:   return "Mat33";
        default: return "ILLEGAL PLACEMENT TYPE";
        };
    }

    SIMTK_REP_HELPERS(Placement,PlacementRep)
private:
    Placement&      handle;    // the Placement whose rep this is
};

/**
 * A PlacementRep whose value is a constant.
 */
class ConstantPlacementRep : public PlacementRep {
public:
    ConstantPlacementRep(ConstantPlacement& cp, bool b) 
      : PlacementRep(cp), bv(b),  t(BoolPlacementType) { }
    ConstantPlacementRep(ConstantPlacement& cp, int  i) 
      : PlacementRep(cp), iv(i),  t(IntPlacementType) { }
    ConstantPlacementRep(ConstantPlacement& cp, Real r) 
      : PlacementRep(cp), rv(r),  t(RealPlacementType) { }
    ConstantPlacementRep(ConstantPlacement& cp, const Vec2&  v) 
      : PlacementRep(cp), v2(v),  t(Vec2PlacementType) { }
    ConstantPlacementRep(ConstantPlacement& cp, const Vec3&  v) 
      : PlacementRep(cp), v3(v),  t(Vec3PlacementType) { }
    ConstantPlacementRep(ConstantPlacement& cp, const Mat33& m) 
      : PlacementRep(cp), m33(m), t(Mat33PlacementType) { }

    ConstantPlacementRep(const ConstantPlacementRep& src) : t(src.t) {
        copyInValidValue(src);
    }
    ~ConstantPlacementRep() { }

    PlacementType getPlacementType() const { return t; }
    bool getValueBool() const {assert(t==BoolPlacementType); return bv;}
    bool getValueInt () const {assert(t==IntPlacementType);  return iv;}
    bool getValueReal() const {assert(t==RealPlacementType); return rv;}
    bool getValueVec2() const {assert(t==Vec2PlacementType); return v2;}
    bool getValueVec3() const {assert(t==Vec3PlacementType); return v3;}
    bool getValueMat33()const {assert(t==Mat33PlacementType);return m33;}

private:
    // avoid references to uninitialized memory
    void copyInValidValue(const ConstantPlacementRep& src) {
        switch(t) {
        case: BoolPlacementType:    bv  = src.bv; break;
        case: IntPlacementType:     iv  = src.iv; break;
        case: RealPlacementType:    rv  = src.rv; break;
        case: Vec2PlacementType:    v2  = src.v2; break;
        case: Vec3PlacementType:    v3  = src.v3; break;
        case: Mat33PlacementType:   m33 = src.m33;break;
        default: assert(false);
        };
    }

private:
    const PlacementType t;
    bool bv;    // only one of these is valid
    int  iv;
    Real rv;
    Vec2 v2;
    Vec3 v3;
    Mat33 m33;
};

/**
 * A PlacementRep whose value is the placement of some feature.
 */
class FeaturePlacementRep : public PlacementRep {
public:
    explicit FeaturePlacementRep(FeaturePlacement& fp, const Feature& f)
        : PlacementRep(fp), feature(f) { }
    ~FeaturePlacementRep() { }
    PlacementRep* clone() const {return FeaturePlacementRep(*this);}
private:
    Feature&    feature;
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
