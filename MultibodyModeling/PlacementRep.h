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

// defined below
class PlacementRep;
class   BoolPlacementRep;
class     BoolConstantPlacementRep;
class     BoolFeaturePlacementRep;
class     BoolExprPlacementRep;
class   IntPlacementRep;
class     IntConstantPlacementRep;
class     IntFeaturePlacementRep;
class     IntExprPlacementRep;
class   RealPlacementRep;
class     RealConstantPlacementRep;
class     RealFeaturePlacementRep;
class     RealExprPlacementRep;
class   Vec2PlacementRep;
class     Vec2ConstantPlacementRep;
class     Vec2FeaturePlacementRep;
class     Vec2ExprPlacementRep;
class   Vec3PlacementRep;
class     Vec3ConstantPlacementRep;
class     Vec3FeaturePlacementRep;
class     Vec3ExprPlacementRep;
class   StationPlacementRep;
class     StationConstantPlacementRep;
class     StationFeaturePlacementRep;
class     StationExprPlacementRep;
class   DirectionPlacementRep;
class     DirectionConstantPlacementRep;
class     DirectionFeaturePlacementRep;
class     DirectionExprPlacementRep;
class   OrientationPlacementRep;
class     OrientationConstantPlacementRep;
class     OrientationFeaturePlacementRep;
class     OrientationExprPlacementRep;
class   PlacementListRep;
class     PlacementListFeatureRep;
class   FramePlacementRep;
class     FrameConstantPlacementRep;
class     FrameFeaturePlacementRep;

// TODO: shouldn't be necessary to have this enum; use virtual methods instead
enum PlacementType {
    InvalidPlacementType = 0,
    VoidPlacementType,
    BoolPlacementType,
    IntPlacementType,
    RealPlacementType,
    Vec2PlacementType,
    Vec3PlacementType,
    Mat33PlacementType,
    StationPlacementType,
    DirectionPlacementType,
    OrientationPlacementType,
    FramePlacementType
};

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
 * Concrete class for operators producing Real placements.
 */
class RealOps : public RealPlacementOp {
public:
    enum OpKind { Negate, Length,                           // unary 
                  Plus, Minus, Times, Divide, Distance };   // binary
    explicit RealOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new RealOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Negate:        p="negate<Real>"; break;
            case Length:        p="length<Real>"; break;
            case Plus:          p="add<Real>"; break;
            case Minus:         p="sub<Real>"; break;
            case Times:         p="mul<Real>"; break;
            case Divide:        p="dvd<Real>"; break;
            case Distance:      p="distance<Real>"; break;
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
    enum OpKind { Scale, Cast,          // unary
                  Plus, Minus };        // binary
    explicit Vec3Ops(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new Vec3Ops(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Scale: p="scale<Vec3>"; break;
            case Cast:  p="cast<Vec3>";  break;
            case Plus:  p="add<Vec3>";   break;
            case Minus: p="sub<Vec3>";   break;
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
    enum OpKind { Cast,                 // unary
                  Plus, Minus };        // binary
    explicit StationOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new StationOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Cast:        p="cast<Station>"; break;
            case Plus:        p="add<Station>";  break;
            case Minus:       p="sub<Station>";  break;
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
            case Normalize:     p="normalize<Direction>"; break;
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
 * Abstract class which represents an operator which returns an Orientation
 * matrix result when applied to an argument list of Placements.
 */
class OrientationPlacementOp : public PlacementOp {
public:
    virtual ~OrientationPlacementOp() { }
    // Run time
    virtual Mat33 apply(/*State,*/const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(OrientationPlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Direction result when applied to Placement
 * arguments of whatever number and type is appropriate for the operator.
 */
class OrientationOps : public OrientationPlacementOp {
public:
    enum OpKind { NoneYet };
    explicit OrientationOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new OrientationOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case NoneYet: p="NoneYet<Orientation>"; break;
            default:      p="UNKNOWN DirectionOp";
        };
        return std::string(p);
    }

    // XXX not yet
    Mat33 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Mat33(0);}

private:
    OpKind op;
};

/**
 * Abstract class which represents an operator which returns a Frame
 * result when applied to an argument list of Placements.
 */
class FramePlacementOp : public PlacementOp {
public:
    virtual ~FramePlacementOp() { }
    // Run time XXX TODO wrong return type; should be numerical Frame
    virtual Mat33 apply(/*State,*/const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(FramePlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Direction result when applied to Placement
 * arguments of whatever number and type is appropriate for the operator.
 */
class FrameOps : public FramePlacementOp {
public:
    enum OpKind { NoneYet };
    explicit FrameOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new FrameOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case NoneYet: p="???<Frame>"; break;
            default:      p="UNKNOWN DirectionOp";
        };
        return std::string(p);
    }

    // XXX not yet TODO this is the wrong return type
    Mat33 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Mat33(0);}

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
        for (size_t i=0; i < args.size(); ++i)
            if (!args[i].isConstant()) return false;
        return true;
    }

    bool exprDependsOn(const Feature& f) const {
        for (size_t i=0; i < args.size(); ++i)
            if (args[i].dependsOn(f))
                return true;
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

/**
 * This class captures the methods common to all Placements which are simply
 * references to the Placement of some Feature, or an indexed element
 * of such a Placement.
 */
class FeatureReference {
protected:
    explicit FeatureReference(const Feature& f, int i = -1);
    // default copy, assignment, destructor

    const Feature& getReferencedFeature() const {
        assert(feature);
        return *feature;
    }

    const Placement& getReferencedPlacement() const {
        assert(feature);
        return feature->getPlacement();
    }
    bool isIndexed() const { assert(feature); return index != -1; }
    int getPlacementIndex() const { assert(feature); return index; }

    bool refIsConstant() const { return false; } // might be, but we can't count on it

    bool refDependsOn(const Feature& f) const {
        assert(feature);
        return feature->dependsOn(f);
    }

    const Feature* refFindAncestorFeature(const Feature& root) const;
    bool           refIsLimitedToSubtree(const Feature& root, const Feature*& offender) const; 
    void           refRepairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);
    std::string    refToString(const std::string& linePrefix) const;

    // Return the required placement type for the referenced feature (after indexing).
    // That isn't necessary the same as the type of the enclosing Placement, which
    // may be performing some kind of transformation (Station to Vec3, say.)
    PlacementType refGetPlacementType() const;

private:
    const Feature* feature;
    const int      index;
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
        case VoidPlacementType:         return "void";
        case BoolPlacementType:         return "bool";
        case IntPlacementType:          return "int";
        case RealPlacementType:         return "Real";
        case StationPlacementType:      return "Station";
        case DirectionPlacementType:    return "Direction";
        case OrientationPlacementType:  return "Orientation";
        case FramePlacementType:        return "Frame";
        case Vec2PlacementType:         return "Vec2";
        case Vec3PlacementType:         return "Vec3";
        case Mat33PlacementType:        return "Mat33";
        default: return "ILLEGAL PLACEMENT TYPE";
        };
    }

    static int getNIndicesAllowed(PlacementType t) {
        switch(t) {
        case VoidPlacementType:         return 0; // can't use at all
        case BoolPlacementType:
        case IntPlacementType:
        case RealPlacementType:         return 1; // no index or index==0 OK

        case Vec2PlacementType:         return 2; // 2 Reals

        case Vec3PlacementType:
        case StationPlacementType:
        case DirectionPlacementType:    return 3; // 3 Reals

        case Mat33PlacementType:        return 3; // 3 Vec3's (columns)
        case OrientationPlacementType:  return 3; // 3 Directions

        case FramePlacementType:        return 2; // Orientation, Station

        default: 
            assert(false);
        };
        //NOTREACHED
        return -1;
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

        case Mat33PlacementType:
            return Vec3PlacementType;

        case OrientationPlacementType: 
            return DirectionPlacementType;

        case FramePlacementType:
            return i==0 ? OrientationPlacementType : StationPlacementType;

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
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Real placement.
 */
class RealFeaturePlacementRep : public RealPlacementRep, public FeatureReference {
public:
    RealFeaturePlacementRep(RealPlacement& p, const Feature& f, int index = -1) 
      : RealPlacementRep(p), FeatureReference(f,index)
    { }
    ~RealFeaturePlacementRep() { }
    
    PlacementRep*  clone() const {return new RealFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return refFindAncestorFeature(f);}
    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    Real getValue(/*State*/) const;

    SIMTK_DOWNCAST2(RealFeaturePlacementRep,RealPlacementRep,PlacementRep);
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
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Vec3 placement.
 */
class Vec3FeaturePlacementRep : public Vec3PlacementRep, public FeatureReference {
public:
    Vec3FeaturePlacementRep(Vec3Placement& p, const Feature& f, int index = -1) 
      : Vec3PlacementRep(p), FeatureReference(f,index)
    { }
    ~Vec3FeaturePlacementRep() { }
    
    PlacementRep*  clone() const {return new Vec3FeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return refFindAncestorFeature(f);}
    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    Vec3 getValue(/*State*/) const;

    SIMTK_DOWNCAST2(Vec3FeaturePlacementRep,Vec3PlacementRep,PlacementRep);
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

    Placement sub(const Placement& r) const {
        if (StationPlacement::isInstanceOf(r))
            return Vec3Placement::minus(getMyStationHandle(), StationPlacement::downcast(r));
        if (Vec3Placement::isInstanceOf(r))
            return Vec3Placement::minus(getMyStationHandle(), Vec3Placement::downcast(r));

        SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "-", getPlacementTypeName(r.getRep().getPlacementType()));
    }

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers(/*State*/)     const = 0;
    SIMTK_DOWNCAST(StationPlacementRep,PlacementRep);
private:
    const StationPlacement& getMyStationHandle() const 
      { return StationPlacement::downcast(getMyHandle()); }
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
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Station placement.
 */
class StationFeaturePlacementRep : public StationPlacementRep, public FeatureReference {
public:
    StationFeaturePlacementRep(StationPlacement& p, const Feature& f, int index = -1) 
      : StationPlacementRep(p), FeatureReference(f,index)
    { }
    ~StationFeaturePlacementRep() { }
    
    PlacementRep*  clone() const {return new StationFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return refFindAncestorFeature(f);}
    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    virtual FramePlacement castToFramePlacement() const {
        if (!isIndexed() 
            && Station::isInstanceOf(getReferencedFeature()) 
            && getReferencedFeature().hasParentFeature()
            && Frame::isInstanceOf(getReferencedFeature().getParentFeature()))
        {
            return FramePlacement(Frame::downcast(getReferencedFeature()
                                                    .getParentFeature()).getOrientation(),
                                  Station::downcast(getReferencedFeature()));
        }

        SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getReferencedFeature().getFullName(),
                     getReferencedFeature().getFeatureTypeName(),
                     "Orientation");
    }

    Vec3 getMeasureNumbers(/*State*/) const;

    SIMTK_DOWNCAST2(StationFeaturePlacementRep,StationPlacementRep,PlacementRep);
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
    DirectionPlacementRep(DirectionPlacement& p) : PlacementRep(p) { }
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
    DirectionConstantPlacementRep(DirectionPlacement& p, const Vec3& v)
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
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Direction placement.
 */
class DirectionFeaturePlacementRep : public DirectionPlacementRep, public FeatureReference {
public:
    DirectionFeaturePlacementRep(DirectionPlacement& p, const Feature& f, int index = -1) 
      : DirectionPlacementRep(p), FeatureReference(f,index)
    { }
    ~DirectionFeaturePlacementRep() { }
    
    PlacementRep*  clone() const {return new DirectionFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return refFindAncestorFeature(f);}
    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    Vec3 getMeasureNumbers(/*State*/) const;

    SIMTK_DOWNCAST2(DirectionFeaturePlacementRep,DirectionPlacementRep,PlacementRep);
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
    OrientationPlacementRep(OrientationPlacement& p) : PlacementRep(p) { }
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
    OrientationConstantPlacementRep(OrientationPlacement& p, const Mat33& m)
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
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Orientation placement.
 */
class OrientationFeaturePlacementRep : public OrientationPlacementRep, public FeatureReference {
public:
    OrientationFeaturePlacementRep(OrientationPlacement& p, const Feature& f, int index = -1) 
      : OrientationPlacementRep(p), FeatureReference(f,index)
    { }
    ~OrientationFeaturePlacementRep() { }
    
    PlacementRep*  clone() const {return new OrientationFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return refFindAncestorFeature(f);}
    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    Mat33 getMeasureNumbers(/*State*/) const;

    SIMTK_DOWNCAST2(OrientationFeaturePlacementRep,OrientationPlacementRep,PlacementRep);
};

/**
 * A concrete PlacementRep whose value is a Orientation expression. This
 * is always Func(List<Placement>). 
 */
class OrientationExprPlacementRep : public OrientationPlacementRep, public PlacementExpr {
public:
    OrientationExprPlacementRep(OrientationPlacement& p, const OrientationPlacementOp&  f, 
                                                         const std::vector<const Placement*>& a) 
      : OrientationPlacementRep(p), PlacementExpr(f,a)
    { }
    ~OrientationExprPlacementRep() { }

    PlacementRep*  clone() const {return new OrientationExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return exprFindAncestorFeature(f);}
    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    Mat33 getMeasureNumbers(/*State*/) const 
      { return OrientationPlacementOp::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(OrientationExprPlacementRep,OrientationPlacementRep,PlacementRep);
};

class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep(FramePlacement& p) : PlacementRep(p) { }
    virtual ~FramePlacementRep() { }


    PlacementType getPlacementType() const { return FramePlacementType; }
    // clone, toString, findAncestorFeature are still missing

    virtual Mat33 getOrientationMeasureNumbers(/*State*/) const = 0;
    virtual Vec3  getOriginMeasureNumbers(/*State*/)      const = 0;

    SIMTK_DOWNCAST(FramePlacementRep,PlacementRep);
};

/**
 * FrameExprPlacementRep an expression with two subexpressions.
 */
class FrameExprPlacementRep : public FramePlacementRep {
public:
    FrameExprPlacementRep(FramePlacement& p, 
        const OrientationPlacement& o, const StationPlacement& s) 
      : FramePlacementRep(p), orientation(o), origin(s) { } 
    ~FrameExprPlacementRep() { }

    bool isConstant() const 
      { return orientation.isConstant() && origin.isConstant(); }

    bool dependsOn(const Feature& f) const {
        return orientation.dependsOn(f) || origin.dependsOn(f);
    }
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const;
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);
    const Feature* findAncestorFeature(const Feature& root) const;

    PlacementRep* clone() const {return new FrameExprPlacementRep(*this);}

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Frame[" << orientation.toString() << ", " << origin.toString() << "]";
        return s.str();
    }

    Mat33 getOrientationMeasureNumbers(/*State*/) const {
        return OrientationPlacementRep::downcast(orientation.getRep())
            .getMeasureNumbers(/*State*/);
    }
    Vec3  getOriginMeasureNumbers(/*State*/)      const {
        return StationPlacementRep::downcast(origin.getRep())
            .getMeasureNumbers(/*State*/);
    }

    SIMTK_DOWNCAST2(FrameExprPlacementRep,FramePlacementRep,PlacementRep);
private:
    OrientationPlacement orientation;
    StationPlacement     origin;
};

/**
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Frame placement.
 */
class FrameFeaturePlacementRep : public FramePlacementRep, public FeatureReference {
public:
    FrameFeaturePlacementRep(FramePlacement& p, const Feature& f, int index = -1) 
      : FramePlacementRep(p), FeatureReference(f,index)
    { }
    ~FrameFeaturePlacementRep() { }
    
    PlacementRep*  clone() const {return new FrameFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Feature* findAncestorFeature(const Feature& f) const {return refFindAncestorFeature(f);}
    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Feature& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Feature& oldRoot, const Feature& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    Mat33 getOrientationMeasureNumbers(/*State*/) const;
    Vec3  getOriginMeasureNumbers(/*State*/)      const;

    SIMTK_DOWNCAST2(FrameFeaturePlacementRep,FramePlacementRep,PlacementRep);
};

} // namespace simtk

#endif // SIMTK_PLACEMENT_REP_H_
