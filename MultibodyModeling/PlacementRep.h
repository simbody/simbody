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
    InvalidPlacementType = 0,
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
    explicit PlacementRep(Placement& p) : myHandle(&p), owner(0), indexInOwner(-1) { }
    virtual ~PlacementRep() { }

    const Placement& getMyHandle()     const {assert(myHandle); return *myHandle;}
    Placement&       updMyHandle()           {assert(myHandle); return *myHandle;}          
    const Feature&   getOwner()        const {assert(owner);    return *owner;}
    int              getIndexInOwner() const {assert(owner);    return indexInOwner;}

    bool hasOwner() const { return owner != 0; }
    void setOwner(const Feature& f, int index) {owner = &f; indexInOwner=index;}

    virtual PlacementType getPlacementType() const = 0;
    virtual void          clone(Placement&)  const = 0;  // clone but with new handle
    virtual std::string   toString(const std::string& linePrefix) const = 0;

    virtual bool isConstant() const { return false; }
    virtual bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { offender=0; return true; }
    virtual void repairFeatureReferences(const Feature& oldRoot, 
                                         const Feature& newRoot) { }
    virtual bool dependsOn(const Feature&) const {return false;}

    // A non-constant Placement may reference many Features, however we expect
    // all of them to be on a common Feature tree. Here we are given the root of
    // the expected tree and return the youngest feature in that tree which is
    // an ancestor of *all* the features in the Placement. Don't call this
    // method on a constant Placement (you can check first with isConstant());
    virtual const Feature* findAncestorFeature(const Feature& root) const = 0;

    static const char* getPlacementTypeName(PlacementType t) {
        switch(t) {
        case InvalidPlacementType:      return "INVALID";
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

protected:
    void cleanUpAfterClone(Placement& p) {
        myHandle = &p;
        p.setRep(this);
    }
private:
    Placement*      myHandle;    // the Placement whose rep this is
    const Feature*  owner;
    int             indexInOwner;
};

/**
 * A concrete PlacementRep whose value is the Placement of some Feature.
 */
class FeaturePlacementRep : public PlacementRep {
public:
    FeaturePlacementRep(FeaturePlacement& p, const Feature& f) 
      : PlacementRep(p), feature(&f) { }
    ~FeaturePlacementRep() { }

    // Check that this feature is on the feature subtree rooted by "ancestor". If
    // not return a pointer to this feature for use in a friendly error message.
    // If this is the right tree, we return true with offender==NULL.
    bool isLimitedToSubtree(const Feature& ancestor, const Feature*& offender) const;
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);

    bool dependsOn(const Feature& f) const {
        assert(feature);
        return feature->dependsOn(f);
    }

    const Feature* findAncestorFeature(const Feature& root) const;
    PlacementType  getPlacementType() const;

    void clone(Placement& handle) const {
        // Note that pointer gets copied as-is. This will require repair when 
        // we copy a Feature tree.
        FeaturePlacementRep* copy = new FeaturePlacementRep(*this);
        copy->cleanUpAfterClone(handle);
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Feature[";
        s << (feature ? feature->getFullName()
                      : std::string("NULL FEATURE"));
        s << "]";   
        return s.str();
    }
    SIMTK_DOWNCAST(FeaturePlacementRep,PlacementRep);
private:
    const Feature* feature;
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
    SIMTK_DOWNCAST(RealPlacementRep,PlacementRep);
};

/**
 * A concrete PlacementRep whose value is a Real constant.
 */
class RealConstantPlacementRep : public RealPlacementRep {
public:
    RealConstantPlacementRep(Placement& p, const Real& r) 
      : RealPlacementRep(p), value(r) { }
    ~RealConstantPlacementRep() { }

    bool isConstant() const { return true; }

    void clone(Placement& p) const {
        RealConstantPlacementRep* copy = new RealConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Real[" << value << "]";   
        return s.str();
    }

    const Feature* findAncestorFeature(const Feature&) const {
        assert(false); // not allowed for constants
        return 0;
    }

    Real getValue(/*State*/) const { return value; }

    SIMTK_DOWNCAST2(RealConstantPlacementRep,RealPlacementRep,PlacementRep);
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
               && args[0]->getRep().getPlacementType()==RealPlacementType
               && args[1]->getRep().getPlacementType()==RealPlacementType;
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

    bool isConstant() const {
        for (size_t i=0; i < args.size(); ++i) {
            assert(args[i]);
            if (!args[i]->isConstant()) return false;
        }
        return true;
    }

    bool dependsOn(const Feature& f) const {
        for (size_t i=0; i < args.size(); ++i) {
            assert(args[i]);
            if (args[i]->dependsOn(f))
                return true;
        }
        return false;
    }

    const Feature* findAncestorFeature(const Feature& root) const;

    void clone(Placement& p) const {
        RealExprPlacementRep* copy = new RealExprPlacementRep(*this);
        copy->cleanUpAfterClone(p);
    }
    std::string toString(const std::string& linePrefix) const {
        std::stringstream s;
        s << func.getOpName() << "(";
        for (size_t i=0; i<args.size(); ++i)
            s << (i>0?", ":"") 
              << args[i]->getRep().toString(linePrefix);
        s << ")";
        return s.str();
    }

    Real getValue(/*State*/) const { return func.apply(/*State,*/args); }

    SIMTK_DOWNCAST2(RealExprPlacementRep,RealPlacementRep,PlacementRep);
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
    SIMTK_DOWNCAST(StationPlacementRep,PlacementRep);
};

// A concrete StationPlacement in which there are no variables.
class StationConstantPlacementRep : public StationPlacementRep {
public:
    StationConstantPlacementRep(Placement& p, const Vec3& v)
      : StationPlacementRep(p), loc(v) { }

    // Implementations of pure virtuals.
    bool isConstant() const { return true; }

    void clone(Placement& p) const {
        StationConstantPlacementRep* copy = new StationConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Station[";
        if (loc == Vec3(0)) s << "0";
        else s << loc;
        s << "]";
        return s.str();
    }

    const Feature* findAncestorFeature(const Feature&) const {
        assert(false); // not allowed for constants
        return 0;
    }

    Vec3 getMeasureNumbers(/*State*/) const { return loc; }

    SIMTK_DOWNCAST2(StationConstantPlacementRep,StationPlacementRep,PlacementRep);
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
    SIMTK_DOWNCAST(DirectionPlacementRep,PlacementRep);
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
    bool isConstant() const { return true; }

    void clone(Placement& p) const {
        DirectionConstantPlacementRep* copy = new DirectionConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Direction[";
        if      (dir == Vec3(1,0,0)) s << "X";
        else if (dir == Vec3(0,1,0)) s << "Y";
        else if (dir == Vec3(0,0,1)) s << "Z";
        else s << dir;
        s << "]";
        return s.str();
    }

    const Feature* findAncestorFeature(const Feature&) const {
        assert(false); // not allowed for constants
        return 0;
    }

    Vec3 getMeasureNumbers(/*State*/) const { return dir; }

    SIMTK_DOWNCAST2(DirectionConstantPlacementRep,DirectionPlacementRep,PlacementRep);
private:
    Vec3 dir;
};


class OrientationPlacementRep : public PlacementRep {
public:
    OrientationPlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~OrientationPlacementRep() { }

    PlacementType getPlacementType() const { return OrientationPlacementType; }
    // clone and toString are still missing

    // These should allow for state to be passed in.
    virtual Mat33  getMeasureNumbers(/*State*/)     const = 0;
    SIMTK_DOWNCAST(OrientationPlacementRep,PlacementRep);
};

// A concrete OrientationPlacement in which there are no variables.
class OrientationConstantPlacementRep : public OrientationPlacementRep {
public:
    OrientationConstantPlacementRep(Placement& p, const Mat33& m)
      : OrientationPlacementRep(p), ori(m) {
        // TODO: check orientation matrix validity
    }

    // Implementations of pure virtuals.
    bool isConstant() const { return true; }

    void clone(Placement& p) const {
        OrientationConstantPlacementRep* copy = new OrientationConstantPlacementRep(*this);
        copy->cleanUpAfterClone(p);
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Orientation[";
        if (ori == Mat33(1)) s << "I";
        else s << ori(0) << ori(1) << ori(2);
        s << "]";
        return s.str();
    }

    const Feature* findAncestorFeature(const Feature&) const {
        assert(false); // not allowed for constants
        return 0;
    }

    Mat33 getMeasureNumbers(/*State*/) const { return ori; }

    SIMTK_DOWNCAST2(OrientationConstantPlacementRep,OrientationPlacementRep,PlacementRep);
private:
    Mat33 ori;
};

/**
 * FramePlacementRep is concrete because it always consists of references
 * to an Orientation Feature and a Station Feature.
 */
class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep(Placement& p, const Orientation& o, const Station& s)
      : PlacementRep(p), orientation(&o), station(&s) { }
    ~FramePlacementRep() { }

    bool dependsOn(const Feature& f) const {
        return orientation->dependsOn(f) || station->dependsOn(f);
    }
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const;
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);

    const Feature* findAncestorFeature(const Feature& root) const;

    PlacementType getPlacementType() const { return FramePlacementType; }
    void clone(Placement& p) const {
        // Note that pointers are copied as-is. These will need repair if
        // we're copying a whole Feature tree.
        FramePlacementRep* copy = new FramePlacementRep(*this);
        copy->cleanUpAfterClone(p);
    }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Frame[";
        s << (orientation ? orientation->getFullName()
                          : std::string("NULL ORIENTATION FEATURE"));
        s << ", ";
        s << (station ? station->getFullName()
                      : std::string("NULL ORIGIN FEATURE"));
        s << "]";
        return s.str();
    }

    SIMTK_DOWNCAST(FramePlacementRep,PlacementRep);
private:
    const Orientation* orientation;
    const Station*     station;
};



} // namespace simtk

#endif // SIMTK_PLACEMENT_REP_H_
