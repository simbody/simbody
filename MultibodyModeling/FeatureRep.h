#ifndef SIMTK_FEATURE_REP_H_
#define SIMTK_FEATURE_REP_H_

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

#include "SimbodyCommon.h"
#include "Feature.h"
#include "Placement.h"
#include "PlacementRep.h"

#include <string>
#include <cassert>
#include <sstream>
#include <cctype>

namespace simtk {

/**
 * This is a "helper" class which is just a Feature with different constructor
 * and assignment behavior. It exists because there are times we must treat
 * the "top level" Feature differently than its children, particularly when
 * we are copying subtrees around. Copying a Feature results in a copy with
 * all its external placements removed and all its internal pointers tidied up.
 * Copying a SubFeature essentially copies all the bits leaving pointers pointing
 * at whatever old junk they were pointing at before.
 */
class SubFeature : public Feature {
public:
    SubFeature() : Feature() { }
    // Copy & assign do *not* invoke the Feature copy constructor.
    inline SubFeature(const SubFeature& sf);
    inline SubFeature& operator=(const SubFeature& sf);
    ~SubFeature() { }
};

/**
 * See SubFeature for an explanation of this class, which is just a Placement
 * with modified copy & assignment behavior.
 */
class SubPlacement : public Placement {
public:
    SubPlacement() : Placement() { }
    // Copy & assign do *not* invoke the Placement copy constructor.
    inline SubPlacement(const SubPlacement& sf);
    inline SubPlacement& operator=(const SubPlacement& sf);
    ~SubPlacement() { }
};

/**
 * A Feature and its FeatureRep are logically part of the same object. There
 * is always a Feature handle (just a pointer) for every FeatureRep, and it
 * must be the case that this->handle.rep == this!
 *
 * FeatureRep is an abstract base class from which the reps for all the
 * concrete Features (e.g. StationRep) derive. The concrete features themselves
 * (e.g. Station) are dataless classes which derive from the concrete Feature
 * class (which contains only a single pointer).
 */
class FeatureRep {
public:
    FeatureRep(Feature& f, const std::string& nm) 
      : myHandle(&f), name(nm), parent(0), indexInParent(-1), placement(0) 
    { }

    // Some Feature types must have a set of Subfeatures installed to complete
    // their construction.
    virtual void initializeStandardSubfeatures() { }

    // Copying a Feature is tricky. The result should have all the child features
    // and the *internal* placements and *internal* placement values.
    // External placements and values should evaporate. Note that the index
    // numbers for features, placements, and values we own must stay the
    // same so that internal references in the copy are the same as in the original.
    // However, this is all handled in the Feature copy & assignment methods --
    // FetureRep copying is elementwise and dumb and thus dangerous. The idea is
    // to get a straight copy and then go clean up the mess afterwards.

    // default (bitwise) copy constructor and assignment -- look out!

    // This is the guts of the smart Feature copy constructor that knows
    // how to clean up all the bad pointers.
    void cloneWithoutParentOrExternalPlacements(Feature& newHandle) const;
    
    virtual ~FeatureRep() { }

    void realize(/*State,*/ Stage g) const;

    void           setMyHandle(Feature& f) {myHandle = &f;}
    const Feature& getMyHandle() const     {assert(myHandle); return *myHandle;}
    Feature&       updMyHandle()           {assert(myHandle); return *myHandle;}

    virtual PlacementType getRequiredPlacementType()      const = 0;
    virtual std::string   getFeatureTypeName()            const = 0;
    virtual FeatureRep*   clone()                         const = 0;

    // Create the appropriate concrete PlacementRep for a reference to the 
    // Placement of this kind of Feature, or to one of its Placement elements
    // if we're given an index (-1 means the whole Placement).
    virtual PlacementRep* createFeatureReference(Placement&, int i = -1) const = 0;

    // Given a proposed placement for this feature, alter it if necessary
    // and return either (1) a Placement that is acceptable, or (2) a
    // Placement with a null rep indicating that the proposed one was no good.
    virtual Placement recastPlacement(const Placement&) const = 0;

    // If this Feature can be used as the indicated placement type, return
    // a new, unowned Placement of the right type. Most commonly, the returned
    // Placement will just be a feature-reference Placement of the same
    // type as the whole Feature, however, for composite Features this may
    // be a reference to one of its subfeatures instead.
    // For example, if a Frame is used as a StationPlacement, we return a
    // reference to the Frame's origin feature.
    // The newly created PlacementRep will refer to the provided Placement handle, but
    // the handle's rep will not be set (otherwise disaster would ensue if
    // we throw an exception somewhere along the way). Be sure to put the
    // returned pointer into the same handle you pass in.

    virtual PlacementRep* useFeatureAsRealPlacement(RealPlacement&) const;
    virtual PlacementRep* useFeatureAsVec3Placement(Vec3Placement&) const;
    virtual PlacementRep* useFeatureAsStationPlacement(StationPlacement&) const;
    virtual PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement&) const;
    virtual PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement&) const;
    virtual PlacementRep* useFeatureAsFramePlacement(FramePlacement&) const;


    void               setName(const std::string& nm) {name = nm;}
    const std::string& getName() const                {return name;}

    void           setParentFeature(Feature& p, int ix) {parent = &p; indexInParent=ix;}
    bool           hasParentFeature() const {return parent != 0;}
    const Feature& getParentFeature() const {assert(hasParentFeature()); return *parent;}
    Feature&       updParentFeature() const {assert(hasParentFeature()); return *parent;}
    int            getIndexInParent() const {assert(hasParentFeature()); return indexInParent;}

    bool             hasPlacement() const {return placement != 0;}

    const Placement& getPlacement() const {
        if (!placement) 
            SIMTK_THROW1(Exception::RepLevelException, "Feature has no placement");
        return *placement;
    }

    void place(const Placement& p);

    int getNSubfeatures()          const {return subfeatures.size();}
    int getNPlacementExpressions() const {return placementExpressions.size();}
    int getNPlacementValues()      const {return placementValues.size();}

    const Feature&        getSubfeature(size_t i)          const {return subfeatures[i];}
    Feature&              updSubfeature(size_t i)                {return subfeatures[i];}
    const Placement&      getPlacementExpression(size_t i) const {return placementExpressions[i];}
    const PlacementValue& getPlacementValue(size_t i)      const {return placementValues[i];}

    const Feature* findSubfeature(const std::string& nm) const {
        size_t index; return findSubfeatureIndex(nm,index) ? &subfeatures[index] : 0;
    }
    Feature* findUpdSubfeature(const std::string& nm) {
        size_t index; return findSubfeatureIndex(nm,index) ? &subfeatures[index] : 0;
    }

    const Feature& getSubfeature(const std::string& nm) const {
        const Feature* f = findSubfeature(nm);
        if (!f) SIMTK_THROW2(Exception::SubfeatureNameNotFound,nm,getFullName());
        return *f;
    }

    Feature& updSubfeature(const std::string& nm) {
        Feature* f = findUpdSubfeature(nm);
        if (!f) SIMTK_THROW2(Exception::SubfeatureNameNotFound,nm,getFullName());
        return *f;
    }

    std::string getFullName() const;

    Feature&        addSubfeatureLike(const Feature& f, const std::string& nm);
    Placement&      addPlacementLike(const Placement& p);
    PlacementValue& addPlacementValueLike(const PlacementValue& v);

    // Does the *placement* of this feature depend on the indicated one?
    // Note that we don't care about our child features' placements.
    bool dependsOn(const Feature& f) const 
        { return placement && placement->dependsOn(f); }

    const Feature& findRootFeature() const;
    Feature&       findUpdRootFeature();       

    // Is Feature f in the tree rooted at oldRoot? If so, optionally return the 
    // series of indices required to get to this Feature from the root.
    // Complexity is O(log n) where n is tree depth.
    static bool isFeatureInFeatureTree(const Feature& oldRoot, const Feature& f,
                                       std::vector<int>* trace=0);

    // Is Placement p owned by a Feature in the tree rooted at oldRoot?
    static bool isPlacementInFeatureTree(const Feature& oldRoot, const Placement& p);

    // If Feature f is a member of the Feature tree rooted at oldRoot, find
    // the corresponding Feature in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Feature* findCorrespondingFeature
        (const Feature& oldRoot, const Feature& f, const Feature& newRoot);

    // If Placement p's owner Feature is a member of the Feature tree rooted at oldRoot,
    // find the corresponding Placement in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Placement* findCorrespondingPlacement
        (const Feature& oldRoot, const Placement& p, const Feature& newRoot);

    // If PlacementValue v's owner Feature is a member of the Feature tree rooted at oldRoot,
    // find the corresponding PlacementValue in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const PlacementValue* findCorrespondingPlacementValue
        (const Feature& oldRoot, const PlacementValue& v, const Feature& newRoot);

    // Given two features, run up the tree towards the root to find
    // their "least common denominator", i.e. the first shared node
    // on the path back to the root. Return a pointer to that node
    // if found, otherwise NULL meaning that the features aren't on
    // the same tree. If the features are the same, then
    // that feature is the answer.
    // Complexity is O(log n) where n is depth of Feature tree.
    static const Feature* findYoungestCommonAncestor(const Feature& f1, const Feature& f2);
    static Feature*       findUpdYoungestCommonAncestor(Feature& f1, const Feature& f2);

    // For debugging
    void checkFeatureConsistency(const Feature* expParent,
                                 int expIndexInParent,
                                 const Feature& expRoot) const;
private:
    // Return true and ix==feature index if a feature of the given name is found.
    // Otherwise return false and ix==childFeatures.size().
    bool findSubfeatureIndex(const std::string& nm, size_t& ix) const;

    // We have just copied a Feature subtree so all the parent pointers are
    // still pointing to the old tree. Recursively repair them to point into
    // the new tree.
    void reparentMyChildren();

    // We have just created at newRoot a copy of the tree rooted at oldRoot, and the
    // current Feature (for which this is the Rep) is a node in the newRoot tree
    // (with correct myHandle). However, the 'parent' and 'placement' pointers 
    // still retain the values they had in the oldRoot tree; they must be 
    // changed to point to the corresponding entities in the newRoot tree.
    // If these pointers point outside the oldRoot tree, however, we'll just
    // set them to 0 in the newRoot copy.
    void fixPlacements(const Feature& oldRoot, const Feature& newRoot);

    
    static const Feature* getParentPtr(const Feature& f) {
        return f.rep ? f.rep->parent : 0;
    }
private:
    Feature*                  myHandle; // the Feature whose rep this is
    std::string               name;

    // Owner information. If parent is 0 then this Feature is an unplaced
    // prototype. Otherwise this Feature is contained in the parent's
    // childFeatures list, with the given index.
    Feature*                  parent;
    int                       indexInParent;

    // If this Feature has been placed, this is the placement information.
    // If present, this Placement must be owned by this Feature, its parent,
    // or one of its ancestors.
    const Placement*          placement;

    // Subfeatures wholly owned by this Feature.
    StableArray<SubFeature>   subfeatures;

    // Placement expressions wholly owned by this Feature. These expressions
    // can involve only Subfeatures of this Feature (or further descendents). But
    // note that we stop at the Subfeature references -- we don't care where they 
    // are placed or even *whether* they are placed. Only when Features are realized
    // do we chase through Placements to calculate values.
    StableArray<SubPlacement> placementExpressions;

    // This is like a State cache except it holds values for Placement expressions
    // where the highest placement dependency is resolved at this Feature (i.e., is
    // one of the placementExpressions stored above. But note that these values
    // do not in general correspond to the placementExpression; they can be the
    // values of lower-level placement expressions.
    StableArray<PlacementValue> placementValues;
};

class RealParameterRep : public FeatureRep {
public:
    RealParameterRep(RealParameter& p, const std::string& nm) : FeatureRep(p,nm) { }
    // no standard Subfeatures

    ~RealParameterRep() { }

    std::string getFeatureTypeName() const { return "RealParameter"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* clone() const { return new RealParameterRep(*this); }
    PlacementRep* createFeatureReference(Placement& p, int i) const {
        if (!(i==-1 || i==0)) {
            SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
                getFullName(), getFeatureTypeName(), i);
            //NOTREACHED
        }
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(p); p.setRep(prep);
        return prep;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToRealPlacement();
    }

    PlacementRep* useFeatureAsRealPlacement(RealPlacement& handle) const {
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(RealParameterRep,FeatureRep);
};

class Vec3ParameterRep : public FeatureRep {
public:
    Vec3ParameterRep(Vec3Parameter& p, const std::string& nm) : FeatureRep(p,nm) { }
    // no standard Subfeatures

    ~Vec3ParameterRep() { }

    std::string getFeatureTypeName() const { return "Vec3Parameter"; }
    PlacementType getRequiredPlacementType() const { return Vec3PlacementType; }
    FeatureRep* clone() const { return new Vec3ParameterRep(*this); }

    PlacementRep* createFeatureReference(Placement& p, int i) const {
        PlacementRep* prep=0;
        if (i == -1) 
            prep = new Vec3FeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);
        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToVec3Placement();
    }

    PlacementRep* useFeatureAsVec3Placement(Vec3Placement& handle) const {
        PlacementRep* prep = new Vec3FeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(Vec3ParameterRep,FeatureRep);
};

class StationParameterRep : public FeatureRep {
public:
    StationParameterRep(StationParameter& p, const std::string& nm) : FeatureRep(p,nm) { }
    // no standard Subfeatures

    ~StationParameterRep() { }

    std::string getFeatureTypeName() const { return "StationParameter"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationParameterRep(*this); }

    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep=0;
        if (i == -1) 
            prep = new StationFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);
        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationParameterRep,FeatureRep);
};

class RealMeasureRep : public FeatureRep {
public:
    RealMeasureRep(RealMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~RealMeasureRep() { }

    std::string getFeatureTypeName() const { return "RealMeasure"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* clone() const { return new RealMeasureRep(*this); }

    PlacementRep* createFeatureReference(Placement& p, int i) const {
        if (!(i==-1 || i==0)) {
            SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
                getFullName(), getFeatureTypeName(), i);
            //NOTREACHED
        }
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(p); p.setRep(prep);
        return prep;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToRealPlacement();
    }

    PlacementRep* useFeatureAsRealPlacement(RealPlacement& handle) const {
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(RealMeasureRep,FeatureRep);
};

class Vec3MeasureRep : public FeatureRep {
public:
    Vec3MeasureRep(Vec3Measure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~Vec3MeasureRep() { }

    std::string getFeatureTypeName() const { return "Vec3Measure"; }
    PlacementType getRequiredPlacementType() const { return Vec3PlacementType; }
    FeatureRep* clone() const { return new Vec3MeasureRep(*this); }

    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep = 0;
        if (i == -1) 
            prep = new Vec3FeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToVec3Placement();
    }

    PlacementRep* useFeatureAsVec3Placement(Vec3Placement& handle) const {
        PlacementRep* prep = new Vec3FeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(Vec3MeasureRep,FeatureRep);
};

class StationMeasureRep : public FeatureRep {
public:
    StationMeasureRep(StationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~StationMeasureRep() { }

    std::string getFeatureTypeName() const { return "StationMeasure"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationMeasureRep(*this); }
 
    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep = 0;
        if (i == -1) 
            prep = new StationFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationMeasureRep,FeatureRep);
};

class StationRep : public FeatureRep {
public:
    StationRep(Station& s, const std::string& nm) : FeatureRep(s,nm) { }
    // no standard Subfeatures

    ~StationRep() { }

    std::string getFeatureTypeName() const { return "Station"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationRep(*this); }

    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep = 0;
        if (i == -1) 
            prep = new StationFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentFeature() && Frame::isInstanceOf(getParentFeature()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Station", "Orientation");
            //NOTREACHED
        }
        const Frame& parentFrame = Frame::downcast(getParentFeature());
        PlacementRep* prep = new FrameExprPlacementRep(parentFrame.getOrientation(), 
                                                       Station::downcast(getMyHandle()));
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationRep,FeatureRep);
};

class DirectionMeasureRep : public FeatureRep {
public:
    DirectionMeasureRep(DirectionMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~DirectionMeasureRep() { }

    std::string getFeatureTypeName() const { return "DirectionMeasure"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    FeatureRep* clone() const { return new DirectionMeasureRep(*this); }
    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep = 0;
        if (i == -1) 
            prep = new DirectionFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToDirectionPlacement();
    }

    PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement& handle) const {
        PlacementRep* prep = new DirectionFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(DirectionMeasureRep,FeatureRep);
};

class DirectionRep : public FeatureRep {
public:
    DirectionRep(Direction& s, const std::string& nm) : FeatureRep(s,nm) { }
    // no standard Subfeatures

    ~DirectionRep() { }

    std::string getFeatureTypeName() const { return "Direction"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    FeatureRep* clone() const { return new DirectionRep(*this); }
    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep=0;
        if (i == -1) 
            prep = new DirectionFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new RealFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToDirectionPlacement();
    }

    PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement& handle) const {
        PlacementRep* prep = new DirectionFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(DirectionRep,FeatureRep);
};


class OrientationMeasureRep : public FeatureRep {
public:
    OrientationMeasureRep(OrientationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~OrientationMeasureRep() { }

    std::string getFeatureTypeName() const { return "OrientationMeasure"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    FeatureRep* clone() const { return new OrientationMeasureRep(*this); }
    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep=0;
        if (i == -1) 
            prep = new OrientationFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new DirectionFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToOrientationPlacement();
    }

    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(OrientationMeasureRep,FeatureRep);
};

class OrientationRep : public FeatureRep {
public:
    OrientationRep(Orientation& o, const std::string& nm) : FeatureRep(o,nm)
      { axisIndices[0]=axisIndices[1]=axisIndices[2] = -1; }
    // must call initializeStandardSubfeatures() to complete construction.

    ~OrientationRep() { }

    std::string   getFeatureTypeName() const { return "Orientation"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    FeatureRep*   clone() const { return new OrientationRep(*this); }

    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep=0;
        if (i == -1) 
            prep = new OrientationFeaturePlacementRep(getMyHandle());
        else if (0<=i && i<3)
            prep = new DirectionFeaturePlacementRep(getMyHandle(), i);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }
        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToOrientationPlacement();
    }

    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentFeature() && Frame::isInstanceOf(getParentFeature()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Orientation", "Station");
            //NOTREACHED
        }
        const Frame& parentFrame = Frame::downcast(getParentFeature());
        PlacementRep* prep = new FrameExprPlacementRep(Orientation::downcast(getMyHandle()),
                                                       parentFrame.getOrigin());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    const Direction& getAxis(int i) const
      { assert(0<=i&&i<=2); return Direction::downcast(getSubfeature(axisIndices[i])); }
    const Direction& x() const {return Direction::downcast(getSubfeature(axisIndices[0]));}
    const Direction& y() const {return Direction::downcast(getSubfeature(axisIndices[1]));}
    const Direction& z() const {return Direction::downcast(getSubfeature(axisIndices[2]));}

    SIMTK_DOWNCAST(OrientationRep,FeatureRep);

protected:
    virtual void initializeStandardSubfeatures() {
        Direction& x = Direction::downcast(addSubfeatureLike(Direction("x"), "x"));
        Direction& y = Direction::downcast(addSubfeatureLike(Direction("y"), "y"));
        Direction& z = Direction::downcast(addSubfeatureLike(Direction("z"), "z"));

        axisIndices[0] = x.getIndexInParent();
        axisIndices[1] = y.getIndexInParent();
        axisIndices[2] = z.getIndexInParent();

        for (int i=0; i<3; ++i)
            updSubfeature(axisIndices[i]).place(Placement(getMyHandle(), i));
    }

private:
    int axisIndices[3];
};

class FrameRep : public FeatureRep {
public:
    FrameRep(Frame& f, const std::string& nm) 
      : FeatureRep(f,nm), RIndex(-1), OIndex(-1) { }
    // must call initializeStandardSubfeatures() to complete construction.

    ~FrameRep() { }

    // still overrideable for bodies.
    virtual std::string   getFeatureTypeName() const { return "Frame"; }
    virtual FeatureRep*   clone() const { return new FrameRep(*this); }

    PlacementType getRequiredPlacementType() const { return FramePlacementType; }

    PlacementRep* createFeatureReference(Placement& p, int i) const { 
        PlacementRep* prep=0;
        if (i == -1) 
            prep = new FrameFeaturePlacementRep(getMyHandle());
        else if (i == 0)
            prep = new OrientationFeaturePlacementRep(getMyHandle(), 0);
        else if (i == 1)
            prep = new StationFeaturePlacementRep(getMyHandle(), 1);

        if (prep) {
            prep->setMyHandle(p); p.setRep(prep);
            return prep;
        }

        SIMTK_THROW3(Exception::IndexOutOfRangeForFeaturePlacementReference,
            getFullName(), getFeatureTypeName(), i);
        //NOTREACHED
        return 0;
    }

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToFramePlacement();
    }

    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        PlacementRep* prep = new FrameFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getOrientation());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getOrigin());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    const Orientation& getOrientation() const {return Orientation::downcast(getSubfeature(RIndex));}
    const Station&     getOrigin()      const {return Station::downcast(getSubfeature(OIndex)); }

    SIMTK_DOWNCAST(FrameRep,FeatureRep);

protected:
    virtual void initializeStandardSubfeatures() {
        Orientation& R = Orientation::downcast(addSubfeatureLike(Orientation("R"), "orientation"));
        Station&     O = Station::downcast    (addSubfeatureLike(Station("O"),     "origin"));

        RIndex = R.getIndexInParent();
        OIndex = O.getIndexInParent();

        updSubfeature(RIndex).place(OrientationPlacement(Mat33(1)));
        updSubfeature(OIndex).place(StationPlacement(Vec3(0)));
    }

private:
    int RIndex, OIndex; // feature indices
};

inline SubFeature::SubFeature(const SubFeature& sf) : Feature() {
    if (sf.rep) { rep = sf.rep->clone(); rep->setMyHandle(*this); }
}
inline SubFeature& SubFeature::operator=(const SubFeature& sf) {
    if (&sf != this) {
        delete rep; rep=0;
        if (sf.rep) { rep = sf.rep->clone(); rep->setMyHandle(*this); }
    }
    return *this;
}

inline SubPlacement::SubPlacement(const SubPlacement& sp) : Placement() {
    if (sp.rep) { rep = sp.rep->clone(); rep->setMyHandle(*this); }
}
inline SubPlacement& SubPlacement::operator=(const SubPlacement& sp) {
    if (&sp != this) {
        delete rep; rep=0;
        if (sp.rep) { rep = sp.rep->clone(); rep->setMyHandle(*this); }
    }
    return *this;
}

} // namespace simtk


#endif // SIMTK_FEATURE_REP_H_
