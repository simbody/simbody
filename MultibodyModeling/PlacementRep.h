#ifndef SIMTK_PLACEMENT_REP_H_
#define SIMTK_PLACEMENT_REP_H_

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
 * The opaque implementation of Placements.
 */

#include "SimbodyCommon.h"
#include "Placement.h"
#include "Feature.h"

#include <string>
#include <vector>
#include <cassert>
#include <sstream>

namespace simtk {

/**
 * Abstract class which represents an operator which returns a Real result
 * when applied to an argument list of Placements.
 */
class RealPlacementOp {
public:
    virtual ~RealPlacementOp() { }
    virtual RealPlacementOp* clone() const = 0;
    virtual bool checkArgs(const std::vector<const Placement*>&) const = 0;

    virtual std::string getOpName() const = 0;

    // Run time
    virtual Real apply(/*State,*/const std::vector<const Placement*>&) const = 0;
};


enum PlacementType {
    UnknownPlacementType = 0,
    BoolPlacementType,
    IntPlacementType,
    RealPlacementType,
    StationPlacementType,
    DirectionPlacementType,
    OrientationPlacementType,
    FramePlacementType,
    Vec2PlacementType,
    Vec3PlacementType,
    Mat33PlacementType
};

class PlacementRep {
public:
    explicit PlacementRep(Placement& p) : handle(&p), owner(0) { }
    virtual ~PlacementRep() { }

    const Placement& getHandle() const {return *handle;}
    const Feature&   getOwner()     const {assert(owner); return *owner;}

    bool hasOwner() const { return owner != 0; }
    void setOwner(const Feature& f, int index) {owner = &f; indexInOwner=index;}

    virtual PlacementType getPlacementType() const = 0;
    virtual PlacementRep* clone(Placement&) const = 0;  // clone but with new handle
    virtual std::string   toString(const std::string& linePrefix) const = 0;

    static const char* getPlacementTypeName(PlacementType t) {
        switch(t) {
        case UnknownPlacementType:      return "unknown";
        case BoolPlacementType:         return "bool";
        case IntPlacementType:          return "int";
        case RealPlacementType:         return "Real";
        case StationPlacementType:      return "Station";
        case DirectionPlacementType:    return "Direction";
        case OrientationPlacementType:  return "Orientation";
        case Vec2PlacementType:         return "Vec2";
        case Vec3PlacementType:         return "Vec3";
        case Mat33PlacementType:        return "Mat33";
        default: return "ILLEGAL PLACEMENT TYPE";
        };
    }

    SIMTK_REP_HELPERS(Placement,PlacementRep)
protected:
    void cleanUpAfterClone(Placement& p) {
        handle = &p;
        // TODO more to come
    }
private:
    Placement*      handle;    // the Placement whose rep this is
    const Feature*  owner;
    size_t          indexInOwner;
};

/**
 * A PlacementRep with a Real value. This is still abstract.
 */
class RealPlacementRep : public PlacementRep {
public:
    RealPlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~RealPlacementRep() { }

    PlacementType getPlacementType() const { return RealPlacementType; }
    // clone and toString are still missing

    // This should allow for state to be passed in.
    virtual Real getValue(/*State*/) const = 0;
};

/**
 * A concrete PlacementRep whose value is a Real constant.
 */
class RealConstantPlacementRep : public RealPlacementRep {
public:
    RealConstantPlacementRep(Placement& p, const Real& r) 
      : RealPlacementRep(p), value(r) { }
    ~RealConstantPlacementRep() { }

    PlacementRep* clone(Placement& p) const {
        RealConstantPlacementRep* copy = new RealConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
        return copy;
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Real[" << value << "]";   
        return s.str();
    }

    Real getValue(/*State*/) const { return value; }
private:
    Real value;
};

class RealBinaryOpRR : public RealPlacementOp {
public:
    enum OpKind { Plus, Minus, Times, Divide };
    explicit RealBinaryOpRR(OpKind k) : op(k) { }

    RealPlacementOp* clone() const { return new RealBinaryOpRR(*this);}
    bool checkArgs(const std::vector<const Placement*>& args) const {
        return args.size() == 2 
               && PlacementRep::getRep(*args[0])->getPlacementType()==RealPlacementType
               && PlacementRep::getRep(*args[1])->getPlacementType()==RealPlacementType;
    }
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Plus:   p="add<Real>"; break;
            case Minus:  p="sub<Real>"; break;
            case Times:  p="mul<Real>"; break;
            case Divide: p="dvd<Real>"; break;
            default: p="UNKNOWN RealBinaryOpRR";
        };
        return std::string(p);
    }

    // XXX not yet
    Real apply(/*State,*/ const std::vector<const Placement*>&) const {assert(false);}

private:
    OpKind op;
};

/**
 * A concrete PlacementRep whose value is a Real expression. This
 * is always Func(List<Placement>). 
 */
class RealExprPlacementRep : public RealPlacementRep {
public:
    RealExprPlacementRep(Placement& p, const RealPlacementOp&  f, 
                                       const std::vector<const Placement*>& a) 
      : RealPlacementRep(p), func(f), args(a) {
        assert(f.checkArgs(args));
    }
    ~RealExprPlacementRep() { }

    PlacementRep* clone(Placement& p) const {
        RealExprPlacementRep* copy = new RealExprPlacementRep(*this);
        copy->cleanUpAfterClone(p);
        return copy;
    }
    std::string toString(const std::string& linePrefix) const {
        std::stringstream s;
        s << func.getOpName() << "(";
        for (size_t i=0; i<args.size(); ++i)
            s << (i>0?", ":"") 
              << PlacementRep::getRep(*args[i])->toString(linePrefix);
        s << ")";
        return s.str();
    }

    Real getValue(/*State*/) const { return func.apply(/*State,*/args); }
private:
    const RealPlacementOp&        func;
    std::vector<const Placement*> args;
};

class StationPlacementRep : public PlacementRep {
public:
    StationPlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~StationPlacementRep() { }

    PlacementType getPlacementType() const { return StationPlacementType; }
    // clone and toString are still missing

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers(/*State*/)     const = 0;
};

// A concrete StationPlacement in which there are no variables.
class StationConstantPlacementRep : public StationPlacementRep {
public:
    StationConstantPlacementRep(Placement& p, const Vec3& v)
      : StationPlacementRep(p), loc(v) { }

    // Implementations of pure virtuals.

    PlacementRep* clone(Placement& p) const {
        StationConstantPlacementRep* copy = new StationConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
        return copy;
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Station[" << loc << "]";
        return s.str();
    }

    Vec3 getMeasureNumbers(/*State*/) const { return loc; }
private:
    Vec3 loc;
};

class DirectionPlacementRep : public PlacementRep {
public:
    DirectionPlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~DirectionPlacementRep() { }

    PlacementType getPlacementType() const { return DirectionPlacementType; }
    // clone and toString are still missing

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers(/*State*/)     const = 0;
};


// A concrete DirectionPlacement in which there are no variables.
class DirectionConstantPlacementRep : public DirectionPlacementRep {
public:
    DirectionConstantPlacementRep(Placement& p, const Vec3& v)
      : DirectionPlacementRep(p), dir(v) {
        const Real len = dir.norm();
        dir /= len;
    }

    // Implementations of pure virtuals.

    PlacementRep* clone(Placement& p) const {
        DirectionConstantPlacementRep* copy = new DirectionConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
        return copy;
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Direction[" << dir << "]";
        return s.str();
    }

    Vec3 getMeasureNumbers(/*State*/) const { return dir; }
private:
    Vec3 dir;
};

class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~FramePlacementRep() { }

    PlacementType getPlacementType() const { return FramePlacementType; }
    // clone and toString are still missing
};


} // namespace simtk

#endif // SIMTK_PLACEMENT_REP_H_
