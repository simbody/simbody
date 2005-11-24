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
 * Concrete class which represents an operator which returns a Real result
 * when applied to an argument list of Placements.
 */
class RealOps : public PlacementOp {
public:
    enum OpKind { Negate, Abs, Sqrt, Exp, Sin, Cos, Asin, Acos, VectorLength, // unary 
                  Add, Subtract, Multiply, Divide, DotProduct2, DotProduct3,  // binary
                  PointDistance };

    explicit RealOps(OpKind k) : op(k) { }
    // default copy, assignment, destructor

    // implementation of pure virtuals
    PlacementOp* clone() const { return new RealOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Negate:        p="negate";   break;
            case Abs:           p="abs";      break;
            case Sqrt:          p="sqrt";     break;
            case Exp:           p="exp";      break;
            case Sin:           p="sin";      break;
            case Cos:           p="cos";      break;
            case Asin:          p="asin";     break;
            case Acos:          p="acos";     break;
            case VectorLength:  p="length";   break;

            case Add:           p="add";      break;
            case Subtract:      p="sub";      break;
            case Multiply:      p="mul";      break;
            case Divide:        p="dvd";      break;
            case DotProduct3:   p="dot3";     break;
            case DotProduct2:   p="dot2";     break;
            case PointDistance: p="distance"; break;
            default:            p="UNKNOWN OP";
        };
        return std::string(p) + "<Real>";
    }

    // Run time calculation of expression value.
    // TODO
    Real apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return 0.;}

    SIMTK_DOWNCAST(RealOps, PlacementOp);
private:
    OpKind op;
};

/**
 * Concrete class producing a Vec3 result when applied to Placement
 * arguments of whatever number and type is appropriate for the operator.
 */
class Vec3Ops : public PlacementOp {
public:
    enum OpKind { RecastStation, RecastDirection, Negate,       // unary
                  Add, Subtract, StationDifference,             // binary
                  ScalarMultiply, ScalarDivide, CrossProduct };
    explicit Vec3Ops(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new Vec3Ops(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case RecastStation:     p="recastStation"; break;
            case RecastDirection:   p="recastDirection"; break;
            case Negate:            p="negate";     break;
            case Add:               p="add";        break;
            case Subtract:          p="sub";        break;
            case StationDifference: p="stationSub"; break;
            case ScalarMultiply:    p="scalarMul";  break;
            case ScalarDivide:      p="scalarDvd";  break;
            case CrossProduct:      p="cross";      break;
            default:                p="UNKNOWN OP";
        };
        return std::string(p) + "<Vec3>";
    }

    // XXX not yet
    Vec3 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Vec3(0);}
    SIMTK_DOWNCAST(Vec3Ops, PlacementOp);
private:
    OpKind op;
};

/**
 * Concrete class producing a Station result when applied to two Placements
 * of whatever type is appropriate for the operator.
 */
class StationOps : public PlacementOp {
public:
    enum OpKind { RecastVec3,      // unary
                  Add, Subtract }; // binary
    explicit StationOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new StationOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case RecastVec3: p="recastVec3"; break;
            case Add:        p="add";    break;   // station = station + vec3
            case Subtract:   p="sub";    break;   // station = station - vec3
            default:         p="UNKNOWN OP";
        };
        return std::string(p) + "<Station>";
    }

    // XXX not yet
    Vec3 apply(/*State,*/ const std::vector<Placement>&) const {assert(false); return Vec3(0);}
    SIMTK_DOWNCAST(StationOps, PlacementOp);
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
 * Not many operators return a Direction (unit vector)!
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
            case Normalize: p="normalize"; break;
            default:        p="UNKNOWN OP";
        };
        return std::string(p) + "<Direction>";
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
            case NoneYet: p="NoneYet"; break;
            default:      p="UNKNOWN OP";
        };
        return std::string(p) + "<Orientation>";
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
            case NoneYet: p="NoneYet"; break;
            default:      p="UNKNOWN OP";
        };
        return std::string(p) + "<Frame>";
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
    // default copy, assignment, destructor

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
        assert(feature); return *feature;
    }
    const Placement& getReferencedPlacement() const {
        assert(feature); return feature->getPlacement();
    }
    bool isIndexed() const { assert(feature); return index != -1; }
    int getPlacementIndex() const { assert(feature); return index; }

    bool refIsConstant() const { return false; } // might be, but we can't count on it

    bool refDependsOn(const Feature& f) const {
        assert(feature); return feature->dependsOn(f);
    }

    const Feature* refFindAncestorFeature(const Feature& root) const;
    bool           refIsLimitedToSubtree(const Feature& root, const Feature*& offender) const; 
    void           refRepairFeatureReferences(const Feature& oldRoot, const Feature& newRoot);
    std::string    refToString(const std::string& linePrefix) const;

    // Return the required placement type for the referenced feature (after indexing).
    // That isn't necessarily the same as the type of the enclosing Placement, which
    // may be performing some kind of transformation (Station to Vec3, say.)
    PlacementType refGetPlacementType() const;

private:
    const Feature* feature;
    const int      index;
};


class PlacementRep {
public:
    explicit PlacementRep() : myHandle(0), owner(0), indexInOwner(-1) { }
    virtual ~PlacementRep() { }

    void             setMyHandle(Placement& p) {myHandle = &p;}
    bool             hasHandle()       const {return myHandle != 0;}
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

    // These are the generic Placement operators, overloaded by argument type.
    // The default implementations provided here throw an exception saying 
    // that the operator is not supported on the given arguments. Any concrete
    // placement that thinks it knows how to implement one of these should 
    // override the default implementation.

    virtual Placement negate()    const;
    virtual Placement abs()       const;
    virtual Placement sqrt()      const;
    virtual Placement sin()       const;
    virtual Placement cos()       const;
    virtual Placement asin()      const;
    virtual Placement acos()      const;
    virtual Placement length()    const;
    virtual Placement normalize() const;
    virtual Placement add(const Placement& r) const;
    virtual Placement sub(const Placement& r) const ;
    virtual Placement mul(const Placement& r) const;
    virtual Placement dvd(const Placement& r) const;
    virtual Placement distance(const Placement& r) const;
    virtual Placement dotProduct  (const Placement& r) const;
    virtual Placement crossProduct(const Placement& r) const;

    virtual RealPlacement        castToRealPlacement()        const;
    virtual Vec3Placement        castToVec3Placement()        const;
    virtual StationPlacement     castToStationPlacement()     const;
    virtual DirectionPlacement   castToDirectionPlacement()   const;
    virtual OrientationPlacement castToOrientationPlacement() const;
    virtual FramePlacement       castToFramePlacement()       const;

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

    static const char* getPlacementTypeName(PlacementType t);
    static int getNIndicesAllowed(PlacementType t);

    // If a PlacementType is indexed, what is the resulting PlacementType?
    static PlacementType getIndexedPlacementType(PlacementType t, int i);

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
    RealPlacementRep() : PlacementRep() { }
    virtual ~RealPlacementRep() { }
    const RealPlacement& getMyHandle() const 
      { return RealPlacement::downcast(PlacementRep::getMyHandle()); }

    PlacementType getPlacementType() const { return RealPlacementType; }
    // clone, toString, findAncestorFeature are still missing

    // These are unary operators on a RealPlacement or binary operators
    // with a RealPlacement on the left, and a Placement of unknown type
    // on the right.

    Placement negate()  const;
    Placement abs()     const;
    Placement sqrt()    const;
    Placement exp()     const;
    Placement log()     const;
    Placement sin()     const;
    Placement cos()     const;
    Placement asin()    const;
    Placement acos()    const;

    RealPlacement castToRealPlacement() const {return RealPlacement(getMyHandle());}

    Placement add(const Placement& r) const;
    Placement sub(const Placement& r) const;
    Placement mul(const Placement& r) const;
    Placement dvd(const Placement& r) const;


    // This should allow for state to be passed in.
    virtual Real getValue(/*State*/) const = 0;
    SIMTK_DOWNCAST(RealPlacementRep,PlacementRep);
};

/**
 * A concrete PlacementRep whose value is a Real constant.
 */
class RealConstantPlacementRep : public RealPlacementRep {
public:
    explicit RealConstantPlacementRep(const Real& r) 
      : RealPlacementRep(), value(r) { }
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
 * Feature which uses a Real placement, or whose placement is Real after
 * indexing.
 */
class RealFeaturePlacementRep : public RealPlacementRep, public FeatureReference {
public:
    explicit RealFeaturePlacementRep(const Feature& f, int index = -1) 
      : RealPlacementRep(), FeatureReference(f,index)
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
    RealExprPlacementRep(const RealOps& f, const std::vector<const Placement*>& a) 
      : RealPlacementRep(), PlacementExpr(f,a)
    { }
    ~RealExprPlacementRep() { }

    // Supported RealExpr-building operators

    static RealExprPlacementRep* negateOp(const RealPlacement&);
    static RealExprPlacementRep* absOp   (const RealPlacement&);
    static RealExprPlacementRep* sqrtOp  (const RealPlacement&);
    static RealExprPlacementRep* expOp   (const RealPlacement&);
    static RealExprPlacementRep* logOp   (const RealPlacement&);
    static RealExprPlacementRep* sinOp   (const RealPlacement&);
    static RealExprPlacementRep* cosOp   (const RealPlacement&);
    static RealExprPlacementRep* asinOp  (const RealPlacement&);
    static RealExprPlacementRep* acosOp  (const RealPlacement&);
    static RealExprPlacementRep* lengthOp(const Vec3Placement&);

    static RealExprPlacementRep* addOp(const RealPlacement& l, const RealPlacement& r);
    static RealExprPlacementRep* subOp(const RealPlacement& l, const RealPlacement& r);
    static RealExprPlacementRep* mulOp(const RealPlacement& l, const RealPlacement& r);
    static RealExprPlacementRep* dvdOp(const RealPlacement& l, const RealPlacement& r);

    static RealExprPlacementRep* distanceOp(const StationPlacement& l, const StationPlacement& r);
    static RealExprPlacementRep* dot2Op    (const Vec2Placement& l,    const Vec2Placement& r);
    static RealExprPlacementRep* dot3Op    (const Vec3Placement& l,    const Vec3Placement& r);
    
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
      { return RealOps::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(RealExprPlacementRep,RealPlacementRep,PlacementRep);
private:
    static RealExprPlacementRep* unaryOp (RealOps::OpKind, const Placement&);
    static RealExprPlacementRep* binaryOp(RealOps::OpKind, const Placement& l, const Placement& r);
};



/**
 * A PlacementRep with a Vec3 value. This is still abstract.
 */
class Vec3PlacementRep : public PlacementRep {
public:
    Vec3PlacementRep() : PlacementRep() { }
    virtual ~Vec3PlacementRep() { }
    const Vec3Placement& getMyHandle() const 
      { return Vec3Placement::downcast(PlacementRep::getMyHandle()); }

    Placement negate() const;
    Placement length() const;
    Placement normalize() const;

    Placement add  (const Placement& r) const;
    Placement sub  (const Placement& r) const;
    Placement mul  (const Placement& r) const;
    Placement dvd  (const Placement& r) const;
    Placement dotProduct  (const Placement& r) const;
    Placement crossProduct(const Placement& r) const;

    PlacementType getPlacementType() const { return Vec3PlacementType; }
    // clone, toString, findAncestorFeature are still missing

    Vec3Placement castToVec3Placement() const {
        return Vec3Placement(getMyHandle()); // we can use it as is!
    }
    DirectionPlacement castToDirectionPlacement() const {
        return DirectionPlacement(normalize());
    }

    // Default inserts a cast operator, but we can do better than that for constants
    // so we'll leave this virtual.
    virtual StationPlacement castToStationPlacement() const 
      { return StationPlacement(getMyHandle()); }

    // This should allow for state to be passed in.
    virtual Vec3 getValue(/*State*/) const = 0;
    SIMTK_DOWNCAST(Vec3PlacementRep,PlacementRep);
private:
};

/**
 * A concrete PlacementRep whose value is a Vec3 constant.
 */
class Vec3ConstantPlacementRep : public Vec3PlacementRep {
public:
    explicit Vec3ConstantPlacementRep(const Vec3& r) : Vec3PlacementRep(), value(r) { }
    ~Vec3ConstantPlacementRep() { }

    StationPlacement castToStationPlacement() const {return StationPlacement(value);}

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
    explicit Vec3FeaturePlacementRep(const Feature& f, int index = -1) 
      : Vec3PlacementRep(), FeatureReference(f,index)
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
    Vec3ExprPlacementRep(const Vec3Ops& f, const std::vector<const Placement*>& a) 
      : Vec3PlacementRep(), PlacementExpr(f,a)
    { }
    ~Vec3ExprPlacementRep() { }

    
    // Supported Vec3Expr-building operators

    static Vec3ExprPlacementRep* negateOp(const Vec3Placement&);

    static Vec3ExprPlacementRep* recastOp(const StationPlacement&);
    static Vec3ExprPlacementRep* recastOp(const DirectionPlacement&);

    static Vec3ExprPlacementRep* addOp (const Vec3Placement& l, const Vec3Placement& r);
    static Vec3ExprPlacementRep* subOp (const Vec3Placement& l, const Vec3Placement& r);
    static Vec3ExprPlacementRep* stationSubOp (const StationPlacement& head, 
                                               const StationPlacement& tail);

    static Vec3ExprPlacementRep* smulOp(const Vec3Placement& l,      const RealPlacement& r);
    static Vec3ExprPlacementRep* smulOp(const StationPlacement& l,   const RealPlacement& r);
    static Vec3ExprPlacementRep* smulOp(const DirectionPlacement& l, const RealPlacement& r);

    static Vec3ExprPlacementRep* sdvdOp(const Vec3Placement& l,      const RealPlacement& r);
    static Vec3ExprPlacementRep* sdvdOp(const StationPlacement& l,   const RealPlacement& r);
    static Vec3ExprPlacementRep* sdvdOp(const DirectionPlacement& l, const RealPlacement& r);

    static Vec3ExprPlacementRep* crossOp(const Vec3Placement& l, const Vec3Placement& r);


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
      { return Vec3Ops::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(Vec3ExprPlacementRep,Vec3PlacementRep,PlacementRep);
private:
    static Vec3ExprPlacementRep* unaryOp (Vec3Ops::OpKind, const Placement&);
    static Vec3ExprPlacementRep* binaryOp(Vec3Ops::OpKind, const Placement& l, const Placement& r);
};

class StationPlacementRep : public PlacementRep {
public:
    StationPlacementRep() : PlacementRep() { }
    virtual ~StationPlacementRep() { }
    const StationPlacement& getMyHandle() const 
      { return StationPlacement::downcast(PlacementRep::getMyHandle()); }

    PlacementType getPlacementType() const { return StationPlacementType; }
    // clone, toString, findAncestorFeature are still missing

    Placement negate()    const;
    Placement length()    const;
    Placement normalize() const;

    Placement add     (const Placement& r) const;
    Placement sub     (const Placement& r) const;
    Placement mul     (const Placement& r) const;
    Placement dvd     (const Placement& r) const;
    Placement dot     (const Placement& r) const;
    Placement cross   (const Placement& r) const;
    Placement distance(const Placement& r) const;

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers(/*State*/)     const = 0;
    SIMTK_DOWNCAST(StationPlacementRep,PlacementRep);

};

// A concrete StationPlacement in which there are no variables.
class StationConstantPlacementRep : public StationPlacementRep {
public:
    StationConstantPlacementRep(const Vec3& v) : StationPlacementRep(), loc(v) { }

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
    explicit StationFeaturePlacementRep(const Feature& f, int index = -1) 
      : StationPlacementRep(), FeatureReference(f,index)
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
    StationExprPlacementRep(const StationOps& f, const std::vector<const Placement*>& a) 
      : StationPlacementRep(), PlacementExpr(f,a)
    { }
    ~StationExprPlacementRep() { }
    
    // Supported StationExpr-building operators

    static StationExprPlacementRep* recastOp(const Vec3Placement&);

    static StationExprPlacementRep* addOp (const StationPlacement& l, const Vec3Placement& r);
    static StationExprPlacementRep* subOp (const StationPlacement& l, const Vec3Placement& r);

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
      { return StationOps::downcast(func).apply(/*State,*/args); }

    SIMTK_DOWNCAST2(StationExprPlacementRep,StationPlacementRep,PlacementRep);
private:
    static StationExprPlacementRep* unaryOp (StationOps::OpKind, const Placement&);
    static StationExprPlacementRep* binaryOp(StationOps::OpKind, const Placement& l, const Placement& r);
};

class DirectionPlacementRep : public PlacementRep {
public:
    DirectionPlacementRep() : PlacementRep() { }
    virtual ~DirectionPlacementRep() { }

    Placement negate()    const;

    Placement add     (const Placement& r) const;
    Placement sub     (const Placement& r) const;
    Placement mul     (const Placement& r) const;
    Placement dvd     (const Placement& r) const;
    Placement dot     (const Placement& r) const;
    Placement cross   (const Placement& r) const;

    PlacementType getPlacementType() const { return DirectionPlacementType; }
    // clone, toString, findAncestorFeature are still missing

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers(/*State*/)     const = 0;
    SIMTK_DOWNCAST(DirectionPlacementRep,PlacementRep);
};


// A concrete DirectionPlacement in which there are no variables.
class DirectionConstantPlacementRep : public DirectionPlacementRep {
public:
    explicit DirectionConstantPlacementRep(const Vec3& v)
      : DirectionPlacementRep(), dir(v) {
        const Real len = dir.norm();
        dir /= len; // let there be NaN's!
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
    explicit DirectionFeaturePlacementRep(const Feature& f, int index = -1) 
      : DirectionPlacementRep(), FeatureReference(f,index)
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
    DirectionExprPlacementRep(const DirectionPlacementOp&  f, 
                              const std::vector<const Placement*>& a) 
      : DirectionPlacementRep(), PlacementExpr(f,a)
    { }
    ~DirectionExprPlacementRep() { }

    // Supported DirectionExpr-building operators

    static DirectionExprPlacementRep* normalizeOp (const StationPlacement&);
    static DirectionExprPlacementRep* normalizeOp (const Vec3Placement&);

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
private:
    static DirectionExprPlacementRep* unaryOp (DirectionOps::OpKind, const Placement&);
    static DirectionExprPlacementRep* binaryOp(DirectionOps::OpKind, const Placement& l, const Placement& r);
};

class OrientationPlacementRep : public PlacementRep {
public:
    OrientationPlacementRep() : PlacementRep() { }
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
    explicit OrientationConstantPlacementRep(const Mat33& m)
      : OrientationPlacementRep(), ori(m) {
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
    explicit OrientationFeaturePlacementRep(const Feature& f, int index = -1) 
      : OrientationPlacementRep(), FeatureReference(f,index)
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
    OrientationExprPlacementRep(const OrientationPlacementOp&  f, 
                                const std::vector<const Placement*>& a) 
      : OrientationPlacementRep(), PlacementExpr(f,a)
    { }
    ~OrientationExprPlacementRep() { }

    // Supported OrientationExpr-building operators
    // NONE YET

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
private:
    static OrientationExprPlacementRep* unaryOp (OrientationOps::OpKind, const Placement&);
    static OrientationExprPlacementRep* binaryOp(OrientationOps::OpKind, 
                                                 const Placement& l, const Placement& r);
};

class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep() : PlacementRep() { }
    virtual ~FramePlacementRep() { }


    PlacementType getPlacementType() const { return FramePlacementType; }
    // clone, toString, findAncestorFeature are still missing

    virtual Mat33 getOrientationMeasureNumbers(/*State*/) const = 0;
    virtual Vec3  getOriginMeasureNumbers(/*State*/)      const = 0;

    SIMTK_DOWNCAST(FramePlacementRep,PlacementRep);
};

/**
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a Frame placement.
 */
class FrameFeaturePlacementRep : public FramePlacementRep, public FeatureReference {
public:
    explicit FrameFeaturePlacementRep(const Feature& f, int index = -1) 
      : FramePlacementRep(), FeatureReference(f,index)
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

/**
 * FrameExprPlacementRep an expression with two subexpressions.
 */
class FrameExprPlacementRep : public FramePlacementRep {
public:
    FrameExprPlacementRep(const OrientationPlacement& o, const StationPlacement& s) 
      : FramePlacementRep(), orientation(o), origin(s) { } 
    ~FrameExprPlacementRep() { }

    // Supported FrameExpr-building operators
    // NONE YET

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

} // namespace simtk

#endif // SIMTK_PLACEMENT_REP_H_
