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
 * Abstract class representing an operator which acts on a list
 * of Placement arguments to produce a placement expression. The
 * result type is unspecified here but concrete classes will
 * return specific types.
 */
class PlacementOp {
public:
    virtual ~PlacementOp() { }
    virtual PlacementOp* clone() const = 0;
    virtual bool checkArgs(const std::vector<Placement>&) const = 0;
    virtual std::string getOpName() const = 0;
};

/**
 * Abstract class which represents an operator which returns a Real result
 * when applied to an argument list of Placements.
 */
class RealPlacementOp : public PlacementOp {
public:
    virtual ~RealPlacementOp() { }
    // Run time
    virtual Real apply(/*State,*/const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(RealPlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Real result when applied to two Real placements.
 */
class RealOps : public RealPlacementOp {
public:
    enum OpKind { Negate, Plus, Minus, Times, Divide };
    explicit RealOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new RealOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Negate: p="negate<Real>"; break;
            case Plus:   p="add<Real>"; break;
            case Minus:  p="sub<Real>"; break;
            case Times:  p="mul<Real>"; break;
            case Divide: p="dvd<Real>"; break;
            default: p="UNKNOWN RealOp";
        };
        return std::string(p);
    }

    // XXX not yet
    Real apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return 0.;}

private:
    OpKind op;
};

/**
 * Abstract class which represents an operator which returns a Vec3 result
 * when applied to an argument list of Placements.
 */
class Vec3PlacementOp : public PlacementOp {
public:
    virtual ~Vec3PlacementOp() { }
    // Run time
    virtual Vec3 apply(/*State,*/const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(Vec3PlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Vec3 result when applied to Placement
 * arguments of whatever number and type is appropriate for the operator.
 */
class Vec3Ops : public Vec3PlacementOp {
public:
    enum OpKind { Scale, Plus, Minus, Cast };
    explicit Vec3Ops(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new Vec3Ops(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Scale: p="scale<Vec3>"; break;
            case Plus:  p="add<Vec3>";   break;
            case Minus: p="sub<Vec3>";   break;
            case Cast:  p="cast<Vec3>";  break;
            default:    p="UNKNOWN Vec3Op";
        };
        return std::string(p);
    }

    // XXX not yet
    Vec3 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Vec3(0);}

private:
    OpKind op;
};

/**
 * Abstract class which represents an operator which returns a Station result
 * when applied to an argument list of Placements.
 */
class StationPlacementOp : public PlacementOp {
public:
    virtual ~StationPlacementOp() { }
    // Run time
    virtual Vec3 apply(/*State,*/const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(StationPlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Station result when applied to two Placements
 * of whatever type is appropriate for the operator.
 */
class StationOps : public StationPlacementOp {
public:
    enum OpKind { Plus, Minus, Cast };
    explicit StationOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new StationOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Plus:  p="add<Station>";  break;
            case Minus: p="sub<Station>";  break;
            case Cast:  p="cast<Station>"; break;
            default: p="UNKNOWN StationOp";
        };
        return std::string(p);
    }

    // XXX not yet
    Vec3 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Vec3(0);}

private:
    OpKind op;
};

/**
 * Abstract class which represents an operator which returns a Direction result
 * when applied to an argument list of Placements.
 */
class DirectionPlacementOp : public PlacementOp {
public:
    virtual ~DirectionPlacementOp() { }
    // Run time
    virtual Vec3 apply(/*State,*/const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(DirectionPlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Direction result when applied to Placement
 * arguments of whatever number and type is appropriate for the operator.
 */
class DirectionOps : public DirectionPlacementOp {
public:
    enum OpKind { Normalize };
    explicit DirectionOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new DirectionOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Normalize: p="normalize<Direction>"; break;
            default:    p="UNKNOWN DirectionOp";
        };
        return std::string(p);
    }

    // XXX not yet
    Vec3 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Vec3(0);}

private:
    OpKind op;
};

/**
 * This class captures the methods common to all Placement expressions, regardless
 * of the specific return type.
 */
class PlacementExpr {
public:
    PlacementExpr(const PlacementOp&  f, const std::vector<const Placement*>& a) 
      : func(f.clone()), args(a.size()) {
        for (size_t i=0; i<a.size(); ++i)
            args[i] = *a[i];
        assert(f.checkArgs(args));
    }
    ~PlacementExpr() { }

    bool exprIsConstant() const {
        for (size_t i=0; i < args.size(); ++i) {
            if (!args[i].isConstant()) return false;
        }
        return true;
    }

    bool exprDependsOn(const Feature& f) const {
        for (size_t i=0; i < args.size(); ++i) {
            if (args[i].dependsOn(f))
                return true;
        }
        return false;
    }

    const Feature* exprFindAncestorFeature(const Feature& root) const;
    std::string    exprToString(const std::string& linePrefix) const;
    bool           exprIsLimitedToSubtree(const Feature& root, const Feature*& offender) const; 
    void           exprRepairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);
protected:
    Concretize<PlacementOp>   func;
    std::vector<Placement>    args;
};

// TODO: shouldn't be necessary to have this enum; use virtual methods instead
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

    void             setMyHandle(Placement& p) {myHandle = &p;}
    const Placement& getMyHandle()     const {assert(myHandle); return *myHandle;}
    Placement&       updMyHandle()           {assert(myHandle); return *myHandle;} 

    void             setOwner(const Feature& f, int index) {owner = &f; indexInOwner=index;}
    bool             hasOwner()        const {return owner != 0;}
    const Feature&   getOwner()        const {assert(owner);    return *owner;}
    int              getIndexInOwner() const {assert(owner);    return indexInOwner;}

    // Note that this copies all feature & placement reference pointers verbatim.
    // The copy will require repair if we are copying a whole Feature tree
    // to get the references to refer to objects in the new tree. 
    void cloneUnownedWithNewHandle(Placement& p) const {
        PlacementRep* pr = clone();
        pr->myHandle = &p; pr->owner = 0; pr->indexInOwner = -1;
        p.setRep(pr);
    }

    virtual Placement negate() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "-", getPlacementTypeName(getPlacementType()));
    }
    virtual RealPlacement length() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "length()", getPlacementTypeName(getPlacementType()));
    }
    virtual DirectionPlacement normalize() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "normalize()", getPlacementTypeName(getPlacementType()));
    }
    virtual Placement add(const Placement& r) const {
        SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                    getPlacementTypeName(getPlacementType()),
                    "+", getPlacementTypeName(r.getRep().getPlacementType()));
    }
    virtual Placement sub(const Placement& r) const {
        SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                    getPlacementTypeName(getPlacementType()),
                    "-", getPlacementTypeName(r.getRep().getPlacementType()));
    }
    virtual Placement mul(const Placement& r) const {
        SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                    getPlacementTypeName(getPlacementType()),
                    "*", getPlacementTypeName(r.getRep().getPlacementType()));
    }
    virtual Placement dvd(const Placement& r) const {
        SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                    getPlacementTypeName(getPlacementType()),
                    "/", getPlacementTypeName(r.getRep().getPlacementType()));
    }

    virtual RealPlacement castToRealPlacement() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "castToRealPlacement()", getPlacementTypeName(getPlacementType()));
    }
    virtual Vec3Placement castToVec3Placement() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "castToVec3Placement()", getPlacementTypeName(getPlacementType()));
    }
    virtual StationPlacement castToStationPlacement() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "castToStationPlacement()", getPlacementTypeName(getPlacementType()));
    }
    virtual DirectionPlacement castToDirectionPlacement() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "castToDirectionPlacement()", getPlacementTypeName(getPlacementType()));
    }
    virtual OrientationPlacement castToOrientationPlacement() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "castToOrientationPlacement()", getPlacementTypeName(getPlacementType()));
    }
    virtual FramePlacement castToFramePlacement() const {
        SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                     "castToFramePlacement()", getPlacementTypeName(getPlacementType()));
    }

    virtual PlacementType getPlacementType() const = 0;
    virtual PlacementRep* clone()            const = 0;
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
        case StationPlacementType:      return "Station";     // 3 Reals
        case DirectionPlacementType:    return "Direction";   // 3 Reals
        case OrientationPlacementType:  return "Orientation"; // 3 Directions
        case FramePlacementType:        return "Frame"; // Orientation, Direction
        case Vec2PlacementType:         return "Vec2";  // 2 Reals
        case Vec3PlacementType:         return "Vec3";  // 3 Reals
        case Mat33PlacementType:        return "Mat33"; // 3 Vec3's (columns)
        default: return "ILLEGAL PLACEMENT TYPE";
        };
    }

    static int getNIndicesAllowed(PlacementType t) {
        switch(t) {
        case StationPlacementType:      return 3; // 3 Reals
        case DirectionPlacementType:    return 3; // 3 Reals
        case OrientationPlacementType:  return 3; // 3 Directions
        case FramePlacementType:        return 2; // Orientation, Station
        case Vec2PlacementType:         return 2; // 2 Reals
        case Vec3PlacementType:         return 3; // 3 Reals
        case Mat33PlacementType:        return 3; // 3 Vec3's (columns)
        default: return 0;
        };
    }

    // If a PlacementType is indexed, what is the resulting PlacementType?
    static PlacementType getIndexedPlacementType(PlacementType t, int i) {
        if (i == -1) 
            return t;   // -1 means not indexed, i.e., the whole thing

        assert(0 <= i && i <= getNIndicesAllowed(t));
        switch(t) {
        case StationPlacementType:
        case DirectionPlacementType: 
        case Vec2PlacementType:
        case Vec3PlacementType:         
            return RealPlacementType;  

        case OrientationPlacementType: 
            return DirectionPlacementType;

        case FramePlacementType:
            return i==0 ? OrientationPlacementType : StationPlacementType;

        case Mat33PlacementType:
            return Vec3PlacementType;

        default: assert(false);
            //NOTREACHED
        };
        //NOTREACHED
        return InvalidPlacementType;
    }

private:
    Placement*      myHandle;    // the Placement whose rep this is
    const Feature*  owner;
    int             indexInOwner;
};

/**
 * A concrete PlacementRep whose value is the Placement of some Feature,
 * or an indexed subcomponent of such a Placement.
 */
class FeaturePlacementRep : public PlacementRep {
public:
//    FeaturePlacementRep(const Feature& f) : PlacementRep(), feature(&f), index(-1) { }
    FeaturePlacementRep(FeaturePlacement& p, const Feature& f) 
      : PlacementRep(p), feature(&f), index(-1) { }
    FeaturePlacementRep(FeaturePlacement& p, const Feature& f, int ix) 
      : PlacementRep(p), feature(&f), index(ix) { }
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

    // Note that pointer gets copied as-is. This will require repair when 
    // we copy a Feature tree.
    PlacementRep* clone() const {return new FeaturePlacementRep(*this);}

    std::string toString(const std::string&) const;
    SIMTK_DOWNCAST(FeaturePlacementRep,PlacementRep);
private:
    const Feature* feature;
    const int      index;
};


/**
 * A PlacementRep with a Real value. This is still abstract.
 */
class RealPlacementRep : public PlacementRep {
public:
    RealPlacementRep(RealPlacement& p) : PlacementRep(p) { }
    virtual ~RealPlacementRep() { }

    PlacementType getPlacementType() const { return RealPlacementType; }
    // clone, toString, findAncestorFeature are still missing

    Placement negate() const {
        return RealPlacement::negate(getMyRealHandle());
    }

    Placement add(const Placement& r) const {
        if (!RealPlacement::isInstanceOf(r))
            SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                        getPlacementTypeName(getPlacementType()),
                        "+", getPlacementTypeName(r.getRep().getPlacementType()));
        return RealPlacement::plus(getMyRealHandle(), RealPlacement::downcast(r));
    }
    Placement sub(const Placement& r) const {
        if (!RealPlacement::isInstanceOf(r))
            SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                         getPlacementTypeName(getPlacementType()),
                         "-", getPlacementTypeName(r.getRep().getPlacementType()));
        return RealPlacement::minus(getMyRealHandle(), RealPlacement::downcast(r));
    }

    Placement mul(const Placement& r) const;

    Placement dvd(const Placement& r) const {
        if (!RealPlacement::isInstanceOf(r))
            SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                         getPlacementTypeName(getPlacementType()),
                         "/", getPlacementTypeName(r.getRep().getPlacementType()));
        return RealPlacement::divide(getMyRealHandle(), RealPlacement::downcast(r));
    }

    RealPlacement castToRealPlacement() const {
        return RealPlacement(getMyRealHandle()); // we can use it as is!
    }

    // This should allow for state to be passed in.
    virtual Real getValue(/*State*/) const = 0;
    SIMTK_DOWNCAST(RealPlacementRep,PlacementRep);
private:
    const RealPlacement& getMyRealHandle() const
      { return RealPlacement::downcast(getMyHandle()); }
};

/**
 * A concrete PlacementRep whose value is a Real constant.
 */
class RealConstantPlacementRep : public RealPlacementRep {
public:
    RealConstantPlacementRep(RealPlacement& p, const Real& r) 
      : RealPlacementRep(p), value(r) { }
    ~RealConstantPlacementRep() { }

    bool isConstant() const { return true; }

    PlacementRep* clone() const {return new RealConstantPlacementRep(*this);}

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


/**
 * A concrete PlacementRep whose value is a Real expression. This
 * is always Func(List<Placement>). 
 */
class RealExprPlacementRep : public RealPlacementRep, public PlacementExpr {
public:
    RealExprPlacementRep(RealPlacement& p, const RealPlacementOp&  f, 
                                           const std::vector<const Placement*>& a) 
      : RealPlacementRep(p), PlacementExpr(f,a)
    { }
    ~RealExprPlacementRep() { }

    static RealExprPlacementRep* unop
        (RealPlacement& handle, RealOps::OpKind, const Placement& a);
    static RealExprPlacementRep* binop
        (RealPlacement& handle, RealOps::OpKind, const Placement& l, const Placement& r);
    
    PlacementRep*  clone() const {return new RealExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return exprFindAncestorFeature(f);}
    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }


    Real getValue(/*State*/) const 
      { return RealPlacementOp::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(RealExprPlacementRep,RealPlacementRep,PlacementRep);
};



/**
 * A PlacementRep with a Vec3 value. This is still abstract.
 */
class Vec3PlacementRep : public PlacementRep {
public:
    Vec3PlacementRep(Vec3Placement& p) : PlacementRep(p) { }
    virtual ~Vec3PlacementRep() { }

    PlacementType getPlacementType() const { return Vec3PlacementType; }
    // clone, toString, findAncestorFeature are still missing

    Vec3Placement castToVec3Placement() const {
        return Vec3Placement(getMyVec3Handle()); // we can use it as is!
    }
    DirectionPlacement castToDirectionPlacement() const {
        return DirectionPlacement::normalize(getMyVec3Handle());
    }

    // Default inserts a cast operator, but we can do better than that for constants
    // so we'll leave this virtual.
    virtual StationPlacement castToStationPlacement() const {
        return StationPlacement::cast(getMyVec3Handle());
    }

    // This should allow for state to be passed in.
    virtual Vec3 getValue(/*State*/) const = 0;
    SIMTK_DOWNCAST(Vec3PlacementRep,PlacementRep);
private:
    const Vec3Placement& getMyVec3Handle() const {return Vec3Placement::downcast(getMyHandle());}
};

/**
 * A concrete PlacementRep whose value is a Vec3 constant.
 */
class Vec3ConstantPlacementRep : public Vec3PlacementRep {
public:
    Vec3ConstantPlacementRep(Vec3Placement& p, const Vec3& r) 
      : Vec3PlacementRep(p), value(r) { }
    ~Vec3ConstantPlacementRep() { }

    StationPlacement castToStationPlacement() const {
        return StationPlacement(value);
    }

    bool isConstant() const { return true; }

    PlacementRep* clone() const {return new Vec3ConstantPlacementRep(*this);}

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Vec3[" << value << "]";   
        return s.str();
    }

    const Feature* findAncestorFeature(const Feature&) const {
        assert(false); // not allowed for constants
        return 0;
    }

    Vec3 getValue(/*State*/) const { return value; }

    SIMTK_DOWNCAST2(Vec3ConstantPlacementRep,Vec3PlacementRep,PlacementRep);
private:
    Vec3 value;
};


/**
 * A concrete PlacementRep whose value is a Vec3 expression. This
 * is always Func(List<Placement>). 
 */
class Vec3ExprPlacementRep : public Vec3PlacementRep, public PlacementExpr {
public:
    Vec3ExprPlacementRep(Vec3Placement& p, const Vec3PlacementOp&  f, 
                                           const std::vector<const Placement*>& a) 
      : Vec3PlacementRep(p), PlacementExpr(f,a)
    { }
    ~Vec3ExprPlacementRep() { }

    static Vec3ExprPlacementRep* scale
        (Vec3Placement& handle, const Placement& l, const Placement& r);
    static Vec3ExprPlacementRep* plus
        (Vec3Placement& handle, const Placement& l, const Placement& r);
    static Vec3ExprPlacementRep* minus
        (Vec3Placement& handle, const Placement& l, const Placement& r);

    static Vec3ExprPlacementRep* cast
        (Vec3Placement& handle, const Placement&);

    PlacementRep*  clone() const {return new Vec3ExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return exprFindAncestorFeature(f);}
    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }


    Vec3 getValue(/*State*/) const 
      { return Vec3PlacementOp::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(Vec3ExprPlacementRep,Vec3PlacementRep,PlacementRep);
};

class StationPlacementRep : public PlacementRep {
public:
    StationPlacementRep(StationPlacement& p) : PlacementRep(p) { }
    virtual ~StationPlacementRep() { }

    PlacementType getPlacementType() const { return StationPlacementType; }
    // clone, toString, findAncestorFeature are still missing

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers(/*State*/)     const = 0;
    SIMTK_DOWNCAST(StationPlacementRep,PlacementRep);
};

// A concrete StationPlacement in which there are no variables.
class StationConstantPlacementRep : public StationPlacementRep {
public:
    StationConstantPlacementRep(StationPlacement& p, const Vec3& v)
      : StationPlacementRep(p), loc(v) { }

    // Implementations of pure virtuals.
    bool isConstant() const { return true; }

    PlacementRep* clone() const {return new StationConstantPlacementRep(*this);}

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

/**
 * A concrete PlacementRep whose value is a Station expression. This
 * is always Func(List<Placement>). 
 */
class StationExprPlacementRep : public StationPlacementRep, public PlacementExpr {
public:
    StationExprPlacementRep(StationPlacement& p, const StationPlacementOp&  f, 
                                                 const std::vector<const Placement*>& a) 
      : StationPlacementRep(p), PlacementExpr(f,a)
    { }
    ~StationExprPlacementRep() { }

    static StationExprPlacementRep* plus
        (StationPlacement& handle, const StationPlacement&, const Vec3Placement&);
    static StationExprPlacementRep* minus
        (StationPlacement& handle, const StationPlacement&, const Vec3Placement&);
    static StationExprPlacementRep* cast
        (StationPlacement& handle, const Vec3Placement&);

    PlacementRep*  clone() const {return new StationExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return exprFindAncestorFeature(f);}
    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    Vec3 getMeasureNumbers(/*State*/) const 
      { return StationPlacementOp::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(StationExprPlacementRep,StationPlacementRep,PlacementRep);
};

class DirectionPlacementRep : public PlacementRep {
public:
    DirectionPlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~DirectionPlacementRep() { }

    PlacementType getPlacementType() const { return DirectionPlacementType; }
    // clone, toString, findAncestorFeature are still missing

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

    PlacementRep* clone() const {return new DirectionConstantPlacementRep(*this);}

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

/**
 * A concrete PlacementRep whose value is a Direction expression. This
 * is always Func(List<Placement>). 
 */
class DirectionExprPlacementRep : public DirectionPlacementRep, public PlacementExpr {
public:
    DirectionExprPlacementRep(DirectionPlacement& p, const DirectionPlacementOp&  f, 
                                                     const std::vector<const Placement*>& a) 
      : DirectionPlacementRep(p), PlacementExpr(f,a)
    { }
    ~DirectionExprPlacementRep() { }

    static DirectionExprPlacementRep* normalize
        (DirectionPlacement& handle, const Vec3Placement& v);
    static DirectionExprPlacementRep* normalize
        (DirectionPlacement& handle, const StationPlacement& v);

    PlacementRep*  clone() const {return new DirectionExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return exprFindAncestorFeature(f);}
    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }


    Vec3 getMeasureNumbers(/*State*/) const 
      { return DirectionPlacementOp::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(DirectionExprPlacementRep,DirectionPlacementRep,PlacementRep);
};

class OrientationPlacementRep : public PlacementRep {
public:
    OrientationPlacementRep(Placement& p) : PlacementRep(p) { }
    virtual ~OrientationPlacementRep() { }

    PlacementType getPlacementType() const { return OrientationPlacementType; }
    // clone, toString, findAncestorFeature are still missing

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

    PlacementRep* clone() const {return new OrientationConstantPlacementRep(*this);}

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
    FramePlacementRep(FramePlacement& p, const Orientation& o, const Station& s)
      : PlacementRep(p), orientation(&o), station(&s) { }
    ~FramePlacementRep() { }

    bool dependsOn(const Feature& f) const {
        return orientation->dependsOn(f) || station->dependsOn(f);
    }
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const;
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);

    const Feature* findAncestorFeature(const Feature& root) const;

    PlacementType getPlacementType() const { return FramePlacementType; }

    // Note that pointers are copied as-is. These will need repair if
    // we're copying a whole Feature tree.
    PlacementRep* clone() const {return new FramePlacementRep(*this);}

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
