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
#include "BasicPlacements.h"
#include "PlacementValue.h"
#include "PlacementValueRep.h"
#include "Subsystem.h"
#include "Feature.h"
#include "BasicFeatures.h"

#include <string>
#include <vector>
#include <cassert>
#include <sstream>

namespace simtk {

class PlacementValueSlot;

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
class   PlacementListRep; // TODO
class     PlacementListFeatureRep;
class   FramePlacementRep;
class     FrameConstantPlacementRep;
class     FrameFeaturePlacementRep;

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
    enum OpKind { Negate, Abs, Sqrt, Exp, Log, Sin, Cos, Asin, Acos, VectorLength, // unary 
                  Add, Subtract, Multiply, Divide, DotProduct2, DotProduct3,  // binary
                  PointDistance, AngleBetweenDirections };

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
            case Log:           p="log";      break;
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
            case AngleBetweenDirections: 
                                p="angle"; break;
            default:            p="UNKNOWN OP";
        };
        return std::string(p) + "<Real>";
    }

    // Run time calculation of expression value.
    Real apply(const std::vector<Placement>&) const;

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

    Vec3 apply(const std::vector<Placement>&) const;

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

    Vec3 apply(const std::vector<Placement>&) const;

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
    virtual Vec3 apply(const std::vector<Placement>&) const = 0;

    SIMTK_DOWNCAST(DirectionPlacementOp, PlacementOp);
};

/**
 * Concrete class producing a Direction result when applied to Placement
 * arguments of whatever number and type is appropriate for the operator.
 * Not many operators return a Direction (unit vector)!
 */
class DirectionOps : public DirectionPlacementOp {
public:
    enum OpKind { Negate, NormalizeVec3, NormalizeStation };
    explicit DirectionOps(OpKind k) : op(k) { }

    PlacementOp* clone() const { return new DirectionOps(*this);}
    bool checkArgs(const std::vector<Placement>& args) const;
    std::string getOpName() const {
        char *p = 0;
        switch(op) {
            case Negate:           p="negate";    break;
            case NormalizeVec3:    p="normalizeVec3"; break;
            case NormalizeStation: p="normalizeStation"; break;
            default:               p="UNKNOWN OP";
        };
        return std::string(p) + "<Direction>";
    }

    Vec3 apply(const std::vector<Placement>&) const;

    SIMTK_DOWNCAST(DirectionOps, PlacementOp);
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
    virtual Mat33 apply(const std::vector<Placement>&) const = 0;

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

    Mat33 apply(const std::vector<Placement>&) const;

    SIMTK_DOWNCAST(OrientationOps, PlacementOp);
private:
    OpKind op;
};

/**
 * Abstract class which represents an operator which returns a FrameFeature
 * result when applied to an argument list of Placements.
 */
class FramePlacementOp : public PlacementOp {
public:
    virtual ~FramePlacementOp() { }
    // Run time
    virtual Mat34 apply(const std::vector<Placement>&) const = 0;

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
        return std::string(p) + "<FrameFeature>";
    }

    Mat34 apply(const std::vector<Placement>&) const;

    SIMTK_DOWNCAST(FrameOps, PlacementOp);
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
      : func(f.clone()), args(a.size())
    {
        for (size_t i=0; i<a.size(); ++i)
            args[i] = *a[i];
        assert(f.checkArgs(args));
    }
    // default copy, assignment, destructor

    const PlacementOp&            exprGetFunc() const {return func.getRef();}
    const std::vector<Placement>& exprGetArgs() const {return args;}

    // Make sure all the arguments are realized.
    void exprRealize(/*State*/Stage) const;

    // Return true if all the arguments are constant.
    bool exprIsConstant() const;

    // Return true if any of the arguments depend on f.
    bool exprDependsOn(const Feature& f) const;

    const Subsystem* exprFindAncestorSubsystem(const Subsystem& youngestAllowed) const;
    const Subsystem* exprFindPlacementValueOwnerSubsystem(const Subsystem& youngestAllowed) const;
    std::string      exprToString(const std::string& linePrefix) const;
    bool             exprIsLimitedToSubtree(const Subsystem& root, const Feature*& offender) const; 
    void             exprRepairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot);
protected:
    const Concretize<PlacementOp>   func;
    std::vector<Placement>          args; // logically const also
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

    const Feature& refGetReferencedFeature() const {
        assert(feature); return *feature;
    }

    bool isIndexed() const { assert(feature); return index != -1; }
    int getPlacementIndex() const { assert(feature); return index; }

    // Make sure this Feature's placement has been realized.
    void refRealize(Stage) const;

    bool refIsConstant() const { return false; } // might be, but we can't count on it
    bool refDependsOn(const Feature& f) const;

    const Subsystem* refFindAncestorSubsystem(const Subsystem& youngestAllowed) const;
    const Subsystem* refFindPlacementValueOwnerSubsystem(const Subsystem& youngestAllowed) const;
    bool             refIsLimitedToSubtree(const Subsystem& root, const Feature*& offender) const; 
    void             refRepairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot);
    std::string      refToString(const std::string& linePrefix) const;

private:
    const Feature* feature;
    const int      index;
};


class PlacementRep {
public:
    explicit PlacementRep() : myHandle(0) { }
    virtual ~PlacementRep() { }

    // This gets the Placement expression ready to evaluate.
    virtual void realize(Stage) const = 0;

    // This performs the calculation and stuffs the result into the
    // supplied PlacementValue.
    virtual void evaluate(PlacementValue&) const = 0;

    // Can we create an instance of this type of Placement by conversion
    // from the Placement p? I.e., Vec3 might say yes if p were a
    // Station and Station might say yes if p were a Frame. Note that these
    // must be resolved at the parent class of a family of PlacementReps
    // which are of the same type.
    virtual PlacementRep* 
        createPlacementFrom(const Placement&, bool dontThrow=false) const = 0;
    virtual bool isSamePlacementType(const Placement&) const = 0;

    virtual std::string    getPlacementTypeName() const = 0;
    virtual int            getNIndicesAllowed()   const = 0;
    virtual PlacementValue createEmptyPlacementValue() const = 0;


    virtual bool isConstant() const { return false; }
    virtual bool isFeatureReference() const { return false; }
    virtual const Feature& getReferencedFeature() const {
        SIMTK_THROW(Exception::NotAFeatureReferencePlacement);
    }

    virtual PlacementRep* clone()            const = 0;
    virtual std::string   toString(const std::string& linePrefix) const = 0;

    virtual Placement genericNegate()    const;
    virtual Placement genericAbs()       const;
    virtual Placement genericSqrt()      const;
    virtual Placement genericExp()       const;
    virtual Placement genericLog()       const;
    virtual Placement genericSin()       const;
    virtual Placement genericCos()       const;
    virtual Placement genericAsin()      const;
    virtual Placement genericAcos()      const;
    virtual Placement genericLength()    const;
    virtual Placement genericNormalize() const;

    virtual Placement genericAdd         (const Placement& rhs) const;
    virtual Placement genericSub         (const Placement& rhs) const;
    virtual Placement genericMul         (const Placement& rhs) const;
    virtual Placement genericDvd         (const Placement& rhs) const;
    virtual Placement genericDistance    (const Placement& rhs) const;
    virtual Placement genericAngle       (const Placement& rhs) const;
    virtual Placement genericDotProduct  (const Placement& rhs) const;
    virtual Placement genericCrossProduct(const Placement& rhs) const;


    virtual bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { offender=0; return true; }
    virtual void repairFeatureReferences(const Subsystem& oldRoot, 
                                         const Subsystem& newRoot) { }
    virtual bool dependsOn(const Feature&) const {return false;}

    // A non-constant Placement may reference many Features, however we expect
    // all of them to be on a common Subsystem tree. Here we are given a subsystem in
    // the expected tree and return the youngest subsystem in that tree which is
    // an ancestor of *all* the features in the Placement, but no younger
    // than the 'youngestAllowed' subsystem supplied in the call, which provides
    // a subsystem that can be used when the Placement doesn't reference any.
    // NOTE: this is to find a suitable owner for the *Placement*, not its
    // *value*. So we are only interested in directly referenced Features; we're
    // not going to recurse through their placements at all.
    virtual const Subsystem* findAncestorSubsystem(const Subsystem& youngestAllowed) const = 0;

    // OK, but this one is a different story. Now we want to find the
    // youngest Subsystem which can hold our PlacementValue cache entry. That
    // Subsystem will be the oldest owner Subsystem for any placement involved
    // in *evaluating* this Placement. That means we DO recurse through the
    // referenced Features' placements. If we encounter any unplaced feature
    // during this search, we'll stop and return NULL because the placement
    // can't yet be evaluated.
    virtual const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestAllowed) const = 0;

    void             setMyHandle(Placement& p) {myHandle = &p;}
    bool             hasHandle()       const {return myHandle != 0;}
    const Placement& getMyHandle()     const {assert(myHandle); return *myHandle;}
    Placement&       updMyHandle()           {assert(myHandle); return *myHandle;} 

    // If a PlacementType is indexed, what is the resulting PlacementType?
    //static PlacementType getIndexedPlacementType(PlacementType t, int i);

private:
    Placement*       myHandle;     // the Placement whose rep this is
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

    PlacementValue createEmptyPlacementValue() const {return PlacementValue_<Real>();}
    std::string    getPlacementTypeName()      const {return "Real";}
    int            getNIndicesAllowed()        const {return 1;} // no index or index==0

    static PlacementRep* createRealPlacementFrom(const Placement&, bool dontThrow=false);
    PlacementRep* createPlacementFrom(const Placement& p, bool dontThrow) const 
      { return createRealPlacementFrom(p,dontThrow); }
    bool isSamePlacementType(const Placement& p) const {return RealPlacement::isInstanceOf(p);}

    // realize, clone, toString, findAncestorSubsystem are still missing

    // These are unary operators on a RealPlacement or binary operators
    // with a RealPlacement on the left and a Placement of unknown type
    // on the right.

    Placement genericNegate()  const;
    Placement genericAbs()     const;
    Placement genericSqrt()    const;
    Placement genericExp()     const;
    Placement genericLog()     const;
    Placement genericSin()     const;
    Placement genericCos()     const;
    Placement genericAsin()    const;
    Placement genericAcos()    const;

    Placement genericAdd(const Placement& r) const;
    Placement genericSub(const Placement& r) const;
    Placement genericMul(const Placement& r) const;
    Placement genericDvd(const Placement& r) const;

    void evaluate(PlacementValue& pv) const {
        PlacementValue_<Real>::initializeToValueType(pv);
        evaluateReal(PlacementValue_<Real>::downcast(pv).upd());
    }
    virtual void evaluateReal(Real&) const = 0;
    Real calcRealValue() const {Real r;evaluateReal(r);return r;}
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

    void realize(Stage) const { } // always ready to evaluate
    void evaluateReal(Real& r) const {r=value;}

    bool isConstant() const { return true; }

    PlacementRep* clone() const {return new RealConstantPlacementRep(*this);}

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Real[" << value << "]";   
        return s.str();
    }

    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }


    SIMTK_DOWNCAST(RealConstantPlacementRep, PlacementRep);
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
      : RealPlacementRep(), FeatureReference(f,index) { }
    ~RealFeaturePlacementRep() { }
    
    void realize(Stage g) const {refRealize(g);}
    void evaluateReal(Real& r) const {r=getReferencedValue();}

    bool isFeatureReference() const {return true;}
    const Feature& getReferencedFeature() const {return refGetReferencedFeature();}

    PlacementRep*  clone() const {return new RealFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)     const {return refToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return refFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return refFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return refIsConstant();}
    bool dependsOn(const Feature& f) const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(RealFeaturePlacementRep, PlacementRep);
private:
    // Get the numerical value of the referenced placement, after indexing.
    const Real& getReferencedValue() const;
};

/**
 * A concrete PlacementRep whose value is a Real expression. This
 * is always Func(List<Placement>). 
 */
class RealExprPlacementRep : public RealPlacementRep, public PlacementExpr {
public:
    RealExprPlacementRep(const RealOps& f, const std::vector<const Placement*>& a) 
      : RealPlacementRep(), PlacementExpr(f,a) { }
    ~RealExprPlacementRep() { }

    void realize(Stage g) const {exprRealize(g);}

    void evaluateReal(Real& r) const {
        r = RealOps::downcast(exprGetFunc()).apply(exprGetArgs());
    }

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

    static RealExprPlacementRep* distanceOp(const StationPlacement&   l, const StationPlacement&   r);
    static RealExprPlacementRep* dot2Op    (const Vec2Placement&      l, const Vec2Placement&      r);
    static RealExprPlacementRep* dot3Op    (const Vec3Placement&      l, const Vec3Placement&      r);
    static RealExprPlacementRep* angleOp   (const DirectionPlacement& l, const DirectionPlacement& r);
    
    PlacementRep*  clone() const {return new RealExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return exprFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return exprFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return exprIsConstant();}
    bool dependsOn(const Feature& f) const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(RealExprPlacementRep, PlacementRep);
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

    PlacementValue createEmptyPlacementValue() const {return PlacementValue_<Vec3>();}
    std::string    getPlacementTypeName()      const {return "Vec3";}
    int            getNIndicesAllowed()        const {return 3;} // 3 reals

    static PlacementRep* createVec3PlacementFrom(const Placement&, bool dontThrow=false);
    PlacementRep* createPlacementFrom(const Placement& p, bool dontThrow) const 
      { return createVec3PlacementFrom(p,dontThrow); }
    bool isSamePlacementType(const Placement& p) const {return Vec3Placement::isInstanceOf(p);}

    Placement genericNegate() const;
    Placement genericLength() const;
    Placement genericNormalize() const;

    Placement genericAdd  (const Placement& rhs) const;
    Placement genericSub  (const Placement& rhs) const;
    Placement genericMul  (const Placement& rhs) const;
    Placement genericDvd  (const Placement& rhs) const;
    Placement genericDotProduct  (const Placement& rhs) const;
    Placement genericCrossProduct(const Placement& rhs) const;
    Placement genericAngle(const Placement& rhs) const;

    // clone, toString, findAncestorSubsystem are still missing

    void evaluate(PlacementValue& pv) const {
        PlacementValue_<Vec3>::initializeToValueType(pv);
        evaluateVec3(PlacementValue_<Vec3>::downcast(pv).upd());
    }
    virtual void evaluateVec3(Vec3&) const = 0;
    Vec3 calcVec3Value() const {Vec3 v;evaluateVec3(v);return v;}
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

    void realize(Stage) const { } // always ready to evaluate
    void evaluateVec3(Vec3& v) const {v=value;}

    bool isConstant() const { return true; }

    PlacementRep* clone() const {return new Vec3ConstantPlacementRep(*this);}

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Vec3[" << value << "]";   
        return s.str();
    }

    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }

    SIMTK_DOWNCAST(Vec3ConstantPlacementRep, PlacementRep);
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
    
    void realize(Stage g) const {refRealize(g);}
    void evaluateVec3(Vec3& v) const {v=getReferencedValue();}

    bool isFeatureReference() const {return true;}
    const Feature& getReferencedFeature() const {return refGetReferencedFeature();}

    PlacementRep*  clone() const {return new Vec3FeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return refFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return refFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return refIsConstant();}
    bool dependsOn(const Feature& f) const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(Vec3FeaturePlacementRep, PlacementRep);
private:
    // Get the numerical value of the referenced placement, after indexing.
    const Vec3& getReferencedValue() const;};

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

    static Vec3ExprPlacementRep* recastStationOp(const StationPlacement&);
    static Vec3ExprPlacementRep* recastDirectionOp(const DirectionPlacement&);

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

    void realize(Stage g) const {exprRealize(g);}
    void evaluateVec3(Vec3& v) const
      { v = Vec3Ops::downcast(exprGetFunc()).apply(exprGetArgs()); }

    PlacementRep*  clone() const {return new Vec3ExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return exprFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return exprFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return exprIsConstant();}
    bool dependsOn(const Feature& f) const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(Vec3ExprPlacementRep, PlacementRep);
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

    PlacementValue createEmptyPlacementValue() const {return PlacementValue_<Vec3>();}
    std::string    getPlacementTypeName()      const {return "Station";}
    int            getNIndicesAllowed()        const {return 3;} // 3 reals

    static PlacementRep* createStationPlacementFrom(const Placement&, bool dontThrow=false);
    PlacementRep* createPlacementFrom(const Placement& p, bool dontThrow) const 
      { return createStationPlacementFrom(p,dontThrow); }
    bool isSamePlacementType(const Placement& p) const {return StationPlacement::isInstanceOf(p);}

    // clone, toString, findAncestorSubsystem are still missing

    Placement genericNegate()    const;
    Placement genericLength()    const;
    Placement genericNormalize() const;

    Placement genericAdd         (const Placement& rhs) const;
    Placement genericSub         (const Placement& rhs) const;
    Placement genericMul         (const Placement& rhs) const;
    Placement genericDvd         (const Placement& rhs) const;
    Placement genericDotProduct  (const Placement& rhs) const;
    Placement genericCrossProduct(const Placement& rhs) const;
    Placement genericDistance    (const Placement& rhs) const;
    Placement genericAngle       (const Placement& rhs) const;

    void evaluate(PlacementValue& pv) const {
        PlacementValue_<Vec3>::initializeToValueType(pv);
        evaluateVec3(PlacementValue_<Vec3>::downcast(pv).upd());
    }
    virtual void evaluateVec3(Vec3&) const = 0;
    Vec3 calcVec3Value() const {Vec3 v;evaluateVec3(v);return v;}

    SIMTK_DOWNCAST(StationPlacementRep,PlacementRep);
};

// A concrete StationPlacement in which there are no variables.
class StationConstantPlacementRep : public StationPlacementRep {
public:
    StationConstantPlacementRep(const Vec3& v) : StationPlacementRep(), loc(v) { }

    // Implementations of pure virtuals.

    void realize(Stage) const { } // always ready to evaluate
    void evaluateVec3(Vec3& v) const {v=loc;}

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

    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }

    SIMTK_DOWNCAST(StationConstantPlacementRep, PlacementRep);
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
    
    void realize(Stage g) const {refRealize(g);}
    void evaluateVec3(Vec3& v) const {v = getReferencedValue();}

    bool isFeatureReference() const {return true;}
    const Feature& getReferencedFeature() const {return refGetReferencedFeature();}

    PlacementRep*  clone() const {return new StationFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return refFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return refFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return refIsConstant();}
    bool dependsOn(const Feature& f) const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(StationFeaturePlacementRep, PlacementRep);
private:
    // Get the numerical value of the referenced placement, after indexing.
    const Vec3& getReferencedValue() const;
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

    static StationExprPlacementRep* recastVec3Op(const Vec3Placement&);

    static StationExprPlacementRep* addOp (const StationPlacement& l, const Vec3Placement& r);
    static StationExprPlacementRep* subOp (const StationPlacement& l, const Vec3Placement& r);

    void realize(Stage g) const {exprRealize(g);}
    void evaluateVec3(Vec3& v) const 
      {  v = StationOps::downcast(exprGetFunc()).apply(exprGetArgs()); }

    PlacementRep*  clone() const {return new StationExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return exprFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return exprFindPlacementValueOwnerSubsystem(s); }

    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(StationExprPlacementRep, PlacementRep);
private:
    static StationExprPlacementRep* unaryOp (StationOps::OpKind, const Placement&);
    static StationExprPlacementRep* binaryOp(StationOps::OpKind, const Placement& l, const Placement& r);
};

class DirectionPlacementRep : public PlacementRep {
public:
    DirectionPlacementRep() : PlacementRep() { }
    virtual ~DirectionPlacementRep() { }
    const DirectionPlacement& getMyHandle() const 
      { return DirectionPlacement::downcast(PlacementRep::getMyHandle()); }

    PlacementValue createEmptyPlacementValue() const {return PlacementValue_<Vec3>();}
    std::string    getPlacementTypeName() const { return "Direction"; }
    int            getNIndicesAllowed()        const {return 3;} // 3 reals

    static PlacementRep* createDirectionPlacementFrom(const Placement&, bool dontThrow=false);
    PlacementRep* createPlacementFrom(const Placement& p, bool dontThrow) const 
      { return createDirectionPlacementFrom(p,dontThrow); }
    bool isSamePlacementType(const Placement& p) const {return DirectionPlacement::isInstanceOf(p);}

    // Negating a direction yields another direction
    Placement genericNegate() const;

    // Direction {+-} Placement doesn't make sense: cast to
    // Vec3 first if that's what you meant.

    // Scaling a Direction by a scalar to produce a Vec3 is OK.
    Placement genericMul         (const Placement& rhs) const;
    Placement genericDvd         (const Placement& rhs) const;

    // Dot(direction,placement) yields Real, Cross yields Vec3
    Placement genericDotProduct  (const Placement& rhs) const;
    Placement genericCrossProduct(const Placement& rhs) const;

    Placement genericAngle       (const Placement& rhs) const;

    // clone, toString, findAncestorSubsystem are still missing

    void evaluate(PlacementValue& pv) const {
        PlacementValue_<Vec3>::initializeToValueType(pv);
        evaluateVec3(PlacementValue_<Vec3>::downcast(pv).upd());
    }
    virtual void evaluateVec3(Vec3&) const = 0;
    Vec3 calcVec3Value() const {Vec3 v;evaluateVec3(v);return v;}

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
    void realize(Stage) const { }   // always ready to evaluate
    void evaluateVec3(Vec3& v) const {v=dir;}

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

    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }

    SIMTK_DOWNCAST(DirectionConstantPlacementRep, PlacementRep);
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

    void realize(Stage g) const {refRealize(g);}
    void evaluateVec3(Vec3& v) const {v = getReferencedValue();}

    bool isFeatureReference() const {return true;}
    const Feature& getReferencedFeature() const {return refGetReferencedFeature();}

    PlacementRep*  clone() const {return new DirectionFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return refFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return refFindPlacementValueOwnerSubsystem(s); }

    bool           isConstant()                          const {return refIsConstant();}
    bool           dependsOn(const Feature& f)           const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(DirectionFeaturePlacementRep, PlacementRep);
private:
    // Get the numerical value of the referenced placement, after indexing.
    const Vec3& getReferencedValue() const;
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

    static DirectionExprPlacementRep* negateOp    (const DirectionPlacement&);
    static DirectionExprPlacementRep* normalizeOp (const StationPlacement&);
    static DirectionExprPlacementRep* normalizeOp (const Vec3Placement&);

    void realize(Stage g) const {exprRealize(g);}
    void evaluateVec3(Vec3& v) const
      { v = DirectionOps::downcast(exprGetFunc()).apply(exprGetArgs()); }

    PlacementRep*  clone() const {return new DirectionExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return exprFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return exprFindPlacementValueOwnerSubsystem(s); }

    bool           isConstant()                          const {return exprIsConstant();}
    bool           dependsOn(const Feature& f)           const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(DirectionExprPlacementRep, PlacementRep);
private:
    static DirectionExprPlacementRep* unaryOp (DirectionOps::OpKind, const Placement&);
    static DirectionExprPlacementRep* binaryOp(DirectionOps::OpKind, const Placement& l, const Placement& r);
};

class OrientationPlacementRep : public PlacementRep {
public:
    OrientationPlacementRep() : PlacementRep() { }
    virtual ~OrientationPlacementRep() { }
    const OrientationPlacement& getMyHandle() const 
      { return OrientationPlacement::downcast(PlacementRep::getMyHandle()); }

    PlacementValue createEmptyPlacementValue() const {return PlacementValue_<Mat33>();}
    std::string    getPlacementTypeName()      const {return "Orientation";}
    int            getNIndicesAllowed()        const {return 3;} // 3 Directions

    static PlacementRep* createOrientationPlacementFrom(const Placement&, bool dontThrow=false);
    PlacementRep* createPlacementFrom(const Placement& p, bool dontThrow) const 
      { return createOrientationPlacementFrom(p,dontThrow); }
    bool isSamePlacementType(const Placement& p) const {return OrientationPlacement::isInstanceOf(p);}

    // clone, toString, findAncestorSubsystem are still missing

    void evaluate(PlacementValue& pv) const {
        PlacementValue_<Mat33>::initializeToValueType(pv);
        evaluateMat33(PlacementValue_<Mat33>::downcast(pv).upd());
    }
    virtual void evaluateMat33(Mat33&) const = 0;
    Mat33 calcMat33Value() const {Mat33 m;evaluateMat33(m);return m;}

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

    void realize(Stage) const { }   // always ready to evaluate
    void evaluateMat33(Mat33& m) const {m=ori;}

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

    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }

    SIMTK_DOWNCAST(OrientationConstantPlacementRep, PlacementRep);
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
      
    void realize(Stage g) const {refRealize(g);}
    void evaluateMat33(Mat33& m) const {m = getReferencedValue();}

    bool isFeatureReference() const {return true;}
    const Feature& getReferencedFeature() const {return refGetReferencedFeature();}

    PlacementRep*  clone() const {return new OrientationFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return refFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return refFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return refIsConstant();}
    bool dependsOn(const Feature& f) const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(OrientationFeaturePlacementRep, PlacementRep);
private:
    // Get the numerical value of the referenced placement, after indexing.
    const Mat33& getReferencedValue() const;
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

    void realize(Stage g) const {exprRealize(g);}
    void evaluateMat33(Mat33& m) const
      { m = OrientationOps::downcast(exprGetFunc()).apply(exprGetArgs()); }

    PlacementRep*  clone() const {return new OrientationExprPlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return exprToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return exprFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return exprFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return exprIsConstant();}
    bool dependsOn(const Feature& f) const {return exprDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return exprIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return exprRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(OrientationExprPlacementRep, PlacementRep);
private:
    static OrientationExprPlacementRep* unaryOp (OrientationOps::OpKind, const Placement&);
    static OrientationExprPlacementRep* binaryOp(OrientationOps::OpKind, 
                                                 const Placement& l, const Placement& r);
};

class FramePlacementRep : public PlacementRep {
public:
    FramePlacementRep() : PlacementRep() { }
    virtual ~FramePlacementRep() { }
    const FramePlacement& getMyHandle() const 
      { return FramePlacement::downcast(PlacementRep::getMyHandle()); }

    PlacementValue createEmptyPlacementValue() const {return PlacementValue_<Mat34>();}
    std::string    getPlacementTypeName()      const {return "Frame";}
    int            getNIndicesAllowed()        const {return 2;} // Orientation, Station

    static PlacementRep* createFramePlacementFrom(const Placement&, bool dontThrow=false);
    PlacementRep* createPlacementFrom(const Placement& p, bool dontThrow) const 
      { return createFramePlacementFrom(p,dontThrow); }
    bool isSamePlacementType(const Placement& p) const {return FramePlacement::isInstanceOf(p);}

    // clone, toString, findAncestorSubsystem are still missing

    void evaluate(PlacementValue& pv) const {
        PlacementValue_<Mat34>::initializeToValueType(pv);
        evaluateMat34(PlacementValue_<Mat34>::downcast(pv).upd());
    }
    virtual void evaluateMat34(Mat34&) const = 0;
    Mat34 calcMat34Value() const {Mat34 m;evaluateMat34(m);return m;}

    SIMTK_DOWNCAST(FramePlacementRep,PlacementRep);
};

// A concrete FramePlacement in which there are no variables.
class FrameConstantPlacementRep : public FramePlacementRep {
public:
    explicit FrameConstantPlacementRep(const Mat34& m)
      : FramePlacementRep(), frame(m) {
        // TODO: check orientation matrix validity
    }

    void realize(Stage) const { }   // always ready to evaluate
    void evaluateMat34(Mat34& m) const {m=frame;}

    bool isConstant() const { return true; }

    PlacementRep* clone() const {return new FrameConstantPlacementRep(*this);}

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "Frame[axes={";
        if (extractOrientation() == Mat33(1)) s << "I";
        else s << extractOrientation()(0) << extractOrientation()(1) << extractOrientation()(2);
        s << "},origin=";
        if (extractOrigin() == Vec3(0)) s << "0";
        else s << extractOrigin();
        s << "]";
        return s.str();
    }

    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const
      { return &youngestSubsystem; }

    SIMTK_DOWNCAST(FrameConstantPlacementRep, PlacementRep);
private:
    Mat34 frame;

    const Mat33& extractOrientation() const {return *reinterpret_cast<const Mat33*>(&frame);}
    const Vec3&  extractOrigin()      const {return frame(3);} // i.e., 4th column
};

/**
 * A concrete PlacementRep whose value is the same as that of a specified
 * Feature which uses a FrameFeature placement.
 */
class FrameFeaturePlacementRep : public FramePlacementRep, public FeatureReference {
public:
    explicit FrameFeaturePlacementRep(const Feature& f, int index = -1) 
      : FramePlacementRep(), FeatureReference(f,index)
    { }
    ~FrameFeaturePlacementRep() { }

    void realize(Stage g) const {refRealize(g);}
    void evaluateMat34(Mat34& m) const {m = getReferencedValue();}

    bool isFeatureReference() const {return true;}
    const Feature& getReferencedFeature() const {return refGetReferencedFeature();}

    PlacementRep*  clone() const {return new FrameFeaturePlacementRep(*this);}
    std::string    toString(const std::string& indent)   const {return refToString(indent);}
    const Subsystem* findAncestorSubsystem(const Subsystem& s) const {return refFindAncestorSubsystem(s);}
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& s) const 
      { return refFindPlacementValueOwnerSubsystem(s); }

    bool isConstant()                const {return refIsConstant();}
    bool dependsOn(const Feature& f) const {return refDependsOn(f);}
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const 
      { return refIsLimitedToSubtree(root,offender); }
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot)
      { return refRepairFeatureReferences(oldRoot, newRoot); }

    SIMTK_DOWNCAST(FrameFeaturePlacementRep, PlacementRep);
private:
    // Get the numerical value of the referenced placement, after indexing.
    const Mat34& getReferencedValue() const;
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
    bool isLimitedToSubtree(const Subsystem& root, const Feature*& offender) const;
    void repairFeatureReferences(const Subsystem& oldRoot, const Subsystem& newRoot);
    const Subsystem* findAncestorSubsystem(const Subsystem& youngestSubsystem) const;
    const Subsystem* findPlacementValueOwnerSubsystem(const Subsystem& youngestSubsystem) const;

    void realize(Stage g) const {
        orientation.getRep().realize(g);
        origin.getRep().realize(g);
    }
    void evaluateMat34(Mat34& m) const {
        Mat33& axesv = reinterpret_cast<Mat33&>(m);
        orientation.getRep().evaluateMat33(axesv);
        origin.getRep().evaluateVec3(m(3));
    }

    PlacementRep* clone() const {return new FrameExprPlacementRep(*this);}

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << "FrameFeature[" << orientation.toString() << ", " << origin.toString() << "]";
        return s.str();
    }

    SIMTK_DOWNCAST(FrameExprPlacementRep, PlacementRep);
private:
    OrientationPlacement orientation;
    StationPlacement     origin;
};


} // namespace simtk

#endif // SIMTK_PLACEMENT_REP_H_
