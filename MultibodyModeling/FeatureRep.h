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
 * This is a "helper" class which is just a Subsystem with different constructor
 * and assignment behavior. It exists because there are times we must treat
 * the "top level" Subsystem differently than its children, particularly when
 * we are copying subtrees around. Copying a Subsystem results in a copy with
 * all its external placements removed and all its internal pointers tidied up.
 * Copying a SubSubsystem essentially copies all the bits leaving pointers pointing
 * at whatever old junk they were pointing at before.
 */
class SubSubsystem : public Subsystem {
public:
    SubSubsystem() : Subsystem() { }
    // Copy & assign do *not* invoke the Feature copy constructor.
    inline SubSubsystem(const SubSubsystem& sf);
    inline SubSubsystem& operator=(const SubSubsystem& sf);
    ~SubSubsystem() { }
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
 * See SubFeature for an explanation of this class, which is just a PlacementValue
 * with modified copy & assignment behavior.
 */
class SubPlacementValue : public PlacementValue {
public:
    SubPlacementValue() : PlacementValue() { }
    // Copy & assign do *not* invoke the Placement copy constructor.
    inline SubPlacementValue(const SubPlacementValue& sf);
    inline SubPlacementValue& operator=(const SubPlacementValue& sf);
    ~SubPlacementValue() { }
};

/**
 * A Subsystem and its SubsystemRep are logically part of the same object. There
 * is always a Subsystem handle (just a pointer) for every SubsystemRep, and there
 * must be exactly one where this->handle.rep == this!
 *
 * SubsystemRep is an abstract base class from which the reps for all the
 * concrete Subsystems (e.g., Feature reps like StationRep) derive.
 * The concrete subsystems themselves (e.g. Station) are dataless classes
 * which derive from the concrete Subsystem class (which contains only a single pointer).
 */
class SubsystemRep {
public:
    SubsystemRep(Subsystem& s, const std::string& nm) 
      : myHandle(&s), name(nm), parent(0), indexInParent(-1) { }
    virtual ~SubsystemRep() { }

    // Some Subsystem types need to get control after we have performed modeling
    // operations that affect them. After construction, they may need to 
    // install some standard subfeatures (or sub-subsystems); after a new subfeature has been
    // added they may need to perform some additional modeling (for example,
    // adding a subfeature with mass may result in an update to the arguments
    // of a 'totalMass' measure); and we offer control after the feature has
    // been placed (that doesn't mean you can necessarily get a *value* for
    // that placement; just the expression defining that value).
    virtual void initializeStandardSubfeatures() { }
    virtual void postProcessNewSubsystem(Subsystem&) { }
    virtual void postProcessNewPlacement() { }

    // This allows a concrete SubsystemRep to prevent the addition
    // of inappropriate sub-subsystems.
    virtual bool canAddSubsystemLike(const Subsystem&)   const 
    {return false;} //TODO: should be pure virtual

    virtual SubsystemRep* clone() const = 0;

    // Copying a Subsystem is tricky. The result should have all the child subsystems
    // and the *internal* feature placements and *internal* placement values.
    // External placements and values should evaporate. Note that the index
    // numbers for subsystems, placements, and values we own must stay the
    // same so that internal references in the copy are the same as in the original.
    // However, this is all handled in the Subsystem copy & assignment methods --
    // SubsystemRep copying is elementwise and dumb and thus dangerous. The idea is
    // to get a straight copy and then go clean up the mess afterwards.

    // default (bitwise) copy constructor and assignment -- look out!

    // This is the guts of the smart Subsystem copy constructor that knows
    // how to clean up all the bad pointers.
    void cloneWithoutParentOrExternalPlacements(Subsystem& newHandle) const;

    void realize(/*State,*/ Stage g) const;

    void             setMyHandle(Subsystem& s) {myHandle = &s;}
    const Subsystem& getMyHandle() const     {assert(myHandle); return *myHandle;}
    Subsystem&       updMyHandle()           {assert(myHandle); return *myHandle;}

    void               setName(const std::string& nm) {name = nm;}
    const std::string& getName() const                {return name;}

    void             setParentSubsystem(Subsystem& p, int ix) {parent = &p; indexInParent=ix;}
    bool             hasParentSubsystem() const {return parent != 0;}
    const Subsystem& getParentSubsystem() const {assert(hasParentSubsystem()); return *parent;}
    Subsystem&       updParentSubsystem() const {assert(hasParentSubsystem()); return *parent;}
    int              getIndexInParent()   const {assert(hasParentSubsystem()); return indexInParent;}

    int getNSubsystems()           const {return childSubsystems.size();}
    int getNPlacementExpressions() const {return placementExpressions.size();}
    int getNPlacementValues()      const {return placementValues.size();}

    const Subsystem&      getSubsystem(size_t i)           const {return childSubsystems[i];}
    Subsystem&            updSubsystem(size_t i)                 {return childSubsystems[i];}

    const Placement&      getPlacementExpression(size_t i) const {return placementExpressions[i];}
    Placement&            updPlacementExpression(size_t i)       {return placementExpressions[i];}

    // PlacementValues here are 'mutable', so the upd routine is const.
    const PlacementValue& getPlacementValue(size_t i)      const {return placementValues[i];}
    PlacementValue&       updPlacementValue(size_t i)      const {return placementValues[i];}

    const Feature& getFeature(size_t i) const {
        const Subsystem& s = getSubsystem(i);
        if (!Feature::isInstanceOf(s)) {
            SIMTK_THROW2(Exception::ExpectedFeatureButGotSubsystem, getFullName(), i);
            //NOTREACHED
        }
        return Feature::downcast(s);
    }
    Feature& updFeature(size_t i) {
        return const_cast<Feature&>(getFeature(i));
    }


    // This can accept a path name
    const Subsystem* findSubsystem(const std::string& nm) const {
        std::vector<std::string> segments;
        if (!isLegalSubsystemPathname(nm, &segments)) {
            if (segments.size())
                SIMTK_THROW2(Exception::IllegalSubsystemPathname,nm,segments.back());
            else
                SIMTK_THROW(Exception::EmptySubsystemPathname);
            //NOTREACHED
        }
        const Subsystem* found = &getMyHandle();
        for (size_t i=0; i<segments.size(); ++i) {
            size_t index;
            found = (found->getRep().findSubsystemIndex(segments[i],index) 
                     ? &found->getRep().childSubsystems[index] : 0);
            if (!found) break;
        }
        return found;
    }
    Subsystem* findUpdSubsystem(const std::string& nm) {
        return const_cast<Subsystem*>(findSubsystem(nm));
    }

    const Feature* findFeature(const std::string& nm) const {
        const Subsystem* sp = findSubsystem(nm);
        assert(sp==0 || Feature::isInstanceOf(*sp));
        return reinterpret_cast<const Feature*>(sp);
    }
    Feature* findUpdFeature(const std::string& nm) {
        return const_cast<Feature*>(findFeature(nm));
    }

    const Subsystem& getSubsystem(const std::string& nm) const {
        const Subsystem* s = findSubsystem(nm);
        if (!s) SIMTK_THROW2(Exception::SubsystemNameNotFound,nm,getFullName());
        return *s;
    }

    Subsystem& updSubsystem(const std::string& nm) {
        Subsystem* s = findUpdSubsystem(nm);
        if (!s) SIMTK_THROW2(Exception::SubsystemNameNotFound,nm,getFullName());
        return *s;
    }

    std::string getFullName() const;

    Subsystem&      addSubsystemLike     (const Subsystem& f, const std::string& nm);
    Feature&        addFeatureLike       (const Subsystem& f, const std::string& nm);
    Placement&      addPlacementLike     (const Placement& p);

    // This is const because PlacementValues are mutable.
    PlacementValue& addPlacementValueLike(const PlacementValue& v) const;

    const Subsystem& findRootSubsystem() const;
    Subsystem&       findUpdRootSubsystem();      

    // Wipe out all PlacementValues. Note that this routine is const because PlacementValues
    // are mutable.
    void deleteAllPlacementValues() const {
        // Erase all references to the placement values (these may be in child subsystems).
        for (size_t i=0; i < placementValues.size(); ++i) {
            PlacementValue& pv = placementValues[i];
            if (pv.getRep().hasClientPlacement())
                pv.getRep().getClientPlacement().getRep().clearValueSlot();
        }
        placementValues.clear();
    }

    // Run around this subsystem and its children looking for Placements which can be
    // evaluated without going any higher in the Subsystem tree than the target Subsysytem, and
    // allocate PlacementValue slots for them in the target.
    void allocatePlacementValueSlots(const Subsystem& target) const {
        // Always let the children go first.
        for (int i=0; i < getNSubsystems(); ++i)
            getSubsystem(i).getRep().allocatePlacementValueSlots(target);

        for (int i=0; i < getNPlacementExpressions(); ++i) {
            const PlacementRep& pr = getPlacementExpression(i).getRep();
            assert(pr.hasOwner());
            if (pr.hasValueSlot())
                continue;
            // can we get it a slot?
            const Subsystem* s = pr.findPlacementValueOwnerSubsystem(target);
            if (s && s->isSameSubsystem(target)) {
                PlacementValue& pv = 
                    s->getRep().addPlacementValueLike(pr.createEmptyPlacementValue());
                pv.updRep().setClientPlacement(getPlacementExpression(i));
                pr.assignValueSlot(pv);
            }
        }
    }

    // Is Subsystem f in the tree rooted at oldRoot? If so, optionally return the 
    // series of indices required to get to this Subsystem from the root.
    // Complexity is O(log n) where n is tree depth.
    static bool isSubsystemInSubsystemTree(const Subsystem& oldRoot, const Subsystem& f,
                                           std::vector<int>* trace=0);

    // Is Placement p owned by a Subsystem in the tree rooted at oldRoot?
    static bool isPlacementInSubsystemTree(const Subsystem& oldRoot, const Placement& p);

    // If Subsystem s is a member of the Subsystem tree rooted at oldRoot, find
    // the corresponding Subsystem in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Subsystem* findCorrespondingSubsystem
        (const Subsystem& oldRoot, const Subsystem& s, const Subsystem& newRoot);

    static const Feature* findCorrespondingFeature
        (const Subsystem& oldRoot, const Feature& s, const Subsystem& newRoot)
    {
        const Subsystem* newSubsys = findCorrespondingSubsystem(oldRoot,s,newRoot);
        assert(newSubsys==0 || Feature::isInstanceOf(*newSubsys));
        return reinterpret_cast<const Feature*>(newSubsys);
    }

    // If Placement p's owner Feature is a member of the Feature tree rooted at oldRoot,
    // find the corresponding Placement in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Placement* findCorrespondingPlacement
        (const Subsystem& oldRoot, const Placement& p, const Subsystem& newRoot);

    // If PlacementValue v's owner Feature is a member of the Feature tree rooted at oldRoot,
    // find the corresponding PlacementValue in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const PlacementValue* findCorrespondingPlacementValue
        (const Subsystem& oldRoot, const PlacementValue& v, const Subsystem& newRoot);

    // Given two subsystems, run up the tree towards the root to find
    // their "least common denominator", i.e. the first shared node
    // on the path back to the root. Return a pointer to that node
    // if found, otherwise NULL meaning that the features aren't on
    // the same tree. If the subsystems are the same, then
    // that feature is the answer.
    // Complexity is O(log n) where n is depth of Subsystem tree.
    static const Subsystem* findYoungestCommonAncestor(const Subsystem& s1, const Subsystem& s2);
    static Subsystem*       findUpdYoungestCommonAncestor(Subsystem& s1, const Subsystem& s2);

    // name utilities
    static bool isLegalSubsystemName(const std::string&);
    static bool isLegalSubsystemPathname(const std::string&, 
                                         std::vector<std::string>* segments=0);

    // For debugging
    void checkSubsystemConsistency(const Subsystem* expParent,
                                   int              expIndexInParent,
                                   const Subsystem& expRoot) const;
private:
    // Return true and ix==subsystem index if a subsystem of the given name is found.
    // Otherwise return false and ix==subsystems.size().
    bool findSubsystemIndex(const std::string& nm, size_t& ix) const;

    // We have just copied a Subsystem subtree so all the parent pointers are
    // still pointing to the old tree. Recursively repair them to point into
    // the new tree.
    void reparentMyChildren();

    // We have just created at newRoot a copy of the tree rooted at oldRoot, and the
    // current Subsystem (for which this is the Rep) is a node in the newRoot tree
    // (with correct myHandle). However, the 'parent' and 'owner' pointers 
    // still retain the values they had in the oldRoot tree; they must be 
    // changed to point to the corresponding entities in the newRoot tree.
    // If these pointers point outside the oldRoot tree, however, we'll just
    // set them to 0 in the newRoot copy.
    void fixPlacements(const Subsystem& oldRoot, const Subsystem& newRoot);
    
    static const Subsystem* getParentPtr(const Subsystem& s) {
        return s.rep ? s.rep->parent : 0;
    }
private:
    Subsystem*                  myHandle; // the Feature whose rep this is
    std::string                 name;

    // Owner information. If parent is 0 then this Subsystem is an unowned
    // prototype. Otherwise this Subsystem is contained in the parent's
    // childSubsystems list, with the given index.
    Subsystem*                  parent;
    int                         indexInParent;

    // Subsystems wholly owned by this Subsystem.
    StableArray<SubSubsystem>   childSubsystems;

    // Placement expressions wholly owned by this Feature. These expressions
    // can involve only Subfeatures of this Feature (or further descendents). But
    // note that we stop at the Subfeature references -- we don't care where they 
    // are placed or even *whether* they are placed. Only when Features are realized
    // do we chase through Placements to calculate values.
    StableArray<SubPlacement>   placementExpressions;

    // This is like a State cache except it holds values for Placement expressions
    // where the highest placement dependency is resolved at this Feature (i.e., is
    // one of the placementExpressions stored above. But note that these values
    // do not in general correspond to the placementExpression; they can be the
    // values of lower-level placement expressions.
    mutable StableArray<SubPlacementValue> placementValues;
};

/**
 * FeatureRep is a still-abstract SubsystemRep which adds handling of the Feature's
 * placement to the basic SubsystemRep capabilities.
 */
class FeatureRep : public SubsystemRep {
public:
    FeatureRep(Feature& p, const std::string& nm)
        : SubsystemRep(p,nm), placement(0) { }
    virtual ~FeatureRep() { }

    // let's be more precise
    const Feature& getMyHandle() const
      { return reinterpret_cast<const Feature&>(SubsystemRep::getMyHandle()); }
    Feature&       updMyHandle()
      { return reinterpret_cast<Feature&>(SubsystemRep::updMyHandle()); }

    // This routine offers control after the feature has
    // been placed (that doesn't mean you can necessarily get a *value* for
    // that placement; just the expression defining that value).
    virtual void postProcessNewPlacement() { }

    // These allow the feature to weigh in on the suitability of a proposed
    // placement for the feature.
    virtual bool canPlaceOnFeatureLike(const Feature&) const
    {return false;} //TODO: should be pure virtual
    virtual bool isRequiredPlacementType(const Placement&) const
    {return false;}
    virtual bool canConvertToRequiredPlacementType(const Placement&) const
    {return false;} //TODO: should be pure virtual

    // Given a proposed placement for this feature, alter it if necessary
    // and return either (1) a Placement that is acceptable, or (2) a
    // Placement with a null rep indicating that the proposed one was no good.
    virtual Placement convertToRequiredPlacementType(const Placement&) const = 0;

    virtual PlacementType getRequiredPlacementType()      const = 0;
    virtual std::string   getFeatureTypeName()            const = 0;

    // Create the appropriate concrete PlacementRep for a reference to the 
    // Placement of this kind of Feature, or to one of its Placement elements
    // if we're given an index (-1 means the whole Placement).
    virtual PlacementRep* createFeatureReference(Placement&, int i = -1) const = 0;

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

    bool             hasPlacement() const {return placement != 0;}

    const Placement& getPlacement() const {
        if (!placement) 
            SIMTK_THROW1(Exception::RepLevelException, "Feature has no placement");
        return *placement;
    }

    void place(const Placement& p);

    // Does the *placement* of this feature depend on the indicated one?
    // Note that we don't care about our child features' placements.
    bool dependsOn(const Feature& f) const 
        { return placement && placement->dependsOn(f); }

    // This is for use by SubsystemRep after a copy to fix the placement pointer.
    void fixFeaturePlacement(const Subsystem& oldRoot, const Subsystem& newRoot);

    SIMTK_DOWNCAST(FeatureRep, SubsystemRep);
private:
    // If this Feature has been placed, this is the placement information.
    // If present, this Placement must be owned by this Feature, its parent
    // Subsystem or one of its ancestors.
    const Placement* placement;
};

class RealParameterRep : public FeatureRep {
public:
    RealParameterRep(RealParameter& p, const std::string& nm) : FeatureRep(p,nm) { }
    // no standard Subfeatures

    ~RealParameterRep() { }

    std::string getFeatureTypeName() const { return "RealParameter"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    SubsystemRep* clone() const { return new RealParameterRep(*this); }
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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToRealPlacement();
    }

    PlacementRep* useFeatureAsRealPlacement(RealPlacement& handle) const {
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(RealParameterRep,SubsystemRep);
};

class Vec3ParameterRep : public FeatureRep {
public:
    Vec3ParameterRep(Vec3Parameter& p, const std::string& nm) : FeatureRep(p,nm) { }
    // no standard Subfeatures

    ~Vec3ParameterRep() { }

    std::string getFeatureTypeName() const { return "Vec3Parameter"; }
    PlacementType getRequiredPlacementType() const { return Vec3PlacementType; }
    SubsystemRep* clone() const { return new Vec3ParameterRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToVec3Placement();
    }

    PlacementRep* useFeatureAsVec3Placement(Vec3Placement& handle) const {
        PlacementRep* prep = new Vec3FeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(Vec3ParameterRep,SubsystemRep);
};

class StationParameterRep : public FeatureRep {
public:
    StationParameterRep(StationParameter& p, const std::string& nm) : FeatureRep(p,nm) { }
    // no standard Subfeatures

    ~StationParameterRep() { }

    std::string getFeatureTypeName() const { return "StationParameter"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    SubsystemRep* clone() const { return new StationParameterRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationParameterRep,SubsystemRep);
};

class RealMeasureRep : public FeatureRep {
public:
    RealMeasureRep(RealMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~RealMeasureRep() { }

    std::string getFeatureTypeName() const { return "RealMeasure"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    SubsystemRep* clone() const { return new RealMeasureRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToRealPlacement();
    }

    PlacementRep* useFeatureAsRealPlacement(RealPlacement& handle) const {
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(RealMeasureRep,SubsystemRep);
};

class Vec3MeasureRep : public FeatureRep {
public:
    Vec3MeasureRep(Vec3Measure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~Vec3MeasureRep() { }

    std::string getFeatureTypeName() const { return "Vec3Measure"; }
    PlacementType getRequiredPlacementType() const { return Vec3PlacementType; }
    SubsystemRep* clone() const { return new Vec3MeasureRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToVec3Placement();
    }

    PlacementRep* useFeatureAsVec3Placement(Vec3Placement& handle) const {
        PlacementRep* prep = new Vec3FeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(Vec3MeasureRep,SubsystemRep);
};

class StationMeasureRep : public FeatureRep {
public:
    StationMeasureRep(StationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~StationMeasureRep() { }

    std::string getFeatureTypeName() const { return "StationMeasure"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    SubsystemRep* clone() const { return new StationMeasureRep(*this); }
 
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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationMeasureRep,SubsystemRep);
};

class StationRep : public FeatureRep {
public:
    StationRep(Station& s, const std::string& nm) : FeatureRep(s,nm) { }
    // no standard Subfeatures

    ~StationRep() { }

    std::string getFeatureTypeName() const { return "Station"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    SubsystemRep* clone() const { return new StationRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentSubsystem() && Frame::isInstanceOf(getParentSubsystem()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Station", "Orientation");
            //NOTREACHED
        }
        const Frame& parentFrame = Frame::downcast(getParentSubsystem());
        PlacementRep* prep = new FrameExprPlacementRep(parentFrame.getOrientation(), 
                                                       Station::downcast(getMyHandle()));
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationRep,SubsystemRep);
};

class DirectionMeasureRep : public FeatureRep {
public:
    DirectionMeasureRep(DirectionMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~DirectionMeasureRep() { }

    std::string getFeatureTypeName() const { return "DirectionMeasure"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    SubsystemRep* clone() const { return new DirectionMeasureRep(*this); }
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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToDirectionPlacement();
    }

    PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement& handle) const {
        PlacementRep* prep = new DirectionFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(DirectionMeasureRep,SubsystemRep);
};

class DirectionRep : public FeatureRep {
public:
    DirectionRep(Direction& s, const std::string& nm) : FeatureRep(s,nm) { }
    // no standard Subfeatures

    ~DirectionRep() { }

    std::string getFeatureTypeName() const { return "Direction"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    SubsystemRep* clone() const { return new DirectionRep(*this); }
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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToDirectionPlacement();
    }

    PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement& handle) const {
        PlacementRep* prep = new DirectionFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(DirectionRep,SubsystemRep);
};


class OrientationMeasureRep : public FeatureRep {
public:
    OrientationMeasureRep(OrientationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }
    // no standard Subfeatures

    ~OrientationMeasureRep() { }

    std::string getFeatureTypeName() const { return "OrientationMeasure"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    SubsystemRep* clone() const { return new OrientationMeasureRep(*this); }
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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToOrientationPlacement();
    }

    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(OrientationMeasureRep,SubsystemRep);
};

class OrientationRep : public FeatureRep {
public:
    OrientationRep(Orientation& o, const std::string& nm) : FeatureRep(o,nm)
      { axisIndices[0]=axisIndices[1]=axisIndices[2] = -1; }
    // must call initializeStandardSubfeatures() to complete construction.

    ~OrientationRep() { }

    std::string   getFeatureTypeName() const { return "Orientation"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    SubsystemRep*   clone() const { return new OrientationRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
        return p.getRep().castToOrientationPlacement();
    }

    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentSubsystem() && Frame::isInstanceOf(getParentSubsystem()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Orientation", "Station");
            //NOTREACHED
        }
        const Frame& parentFrame = Frame::downcast(getParentSubsystem());
        PlacementRep* prep = new FrameExprPlacementRep(Orientation::downcast(getMyHandle()),
                                                       parentFrame.getOrigin());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    const Direction& getAxis(int i) const
      { assert(0<=i&&i<=2); return Direction::downcast(getFeature(axisIndices[i])); }
    const Direction& x() const {return Direction::downcast(getFeature(axisIndices[0]));}
    const Direction& y() const {return Direction::downcast(getFeature(axisIndices[1]));}
    const Direction& z() const {return Direction::downcast(getFeature(axisIndices[2]));}

    SIMTK_DOWNCAST(OrientationRep,SubsystemRep);

protected:
    virtual void initializeStandardSubfeatures() {
        Direction& x = Direction::downcast(addFeatureLike(Direction("x"), "x"));
        Direction& y = Direction::downcast(addFeatureLike(Direction("y"), "y"));
        Direction& z = Direction::downcast(addFeatureLike(Direction("z"), "z"));

        axisIndices[0] = x.getIndexInParent();
        axisIndices[1] = y.getIndexInParent();
        axisIndices[2] = z.getIndexInParent();

        for (int i=0; i<3; ++i)
            updFeature(axisIndices[i]).place(Placement(getMyHandle(), i));
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
    virtual SubsystemRep*   clone() const { return new FrameRep(*this); }

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

    Placement convertToRequiredPlacementType(const Placement& p) const {
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
    const Orientation& getOrientation() const {return Orientation::downcast(getFeature(RIndex));}
    const Station&     getOrigin()      const {return Station::downcast(getFeature(OIndex)); }

    SIMTK_DOWNCAST(FrameRep,SubsystemRep);

protected:
    virtual void initializeStandardSubfeatures() {
        Orientation& R = Orientation::downcast(addFeatureLike(Orientation("R"), "orientation"));
        Station&     O = Station::downcast    (addFeatureLike(Station("O"),     "origin"));

        RIndex = R.getIndexInParent();
        OIndex = O.getIndexInParent();

        updFeature(RIndex).place(OrientationPlacement(Mat33(1)));
        updFeature(OIndex).place(StationPlacement(Vec3(0)));
    }

private:
    int RIndex, OIndex; // feature indices
};

inline SubSubsystem::SubSubsystem(const SubSubsystem& sf) : Subsystem() {
    if (sf.rep) { rep = sf.rep->clone(); rep->setMyHandle(*this); }
}
inline SubSubsystem& SubSubsystem::operator=(const SubSubsystem& sf) {
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

inline SubPlacementValue::SubPlacementValue(const SubPlacementValue& sp) : PlacementValue() {
    if (sp.rep) { rep = sp.rep->clone(); rep->setMyHandle(*this); }
}
inline SubPlacementValue& SubPlacementValue::operator=(const SubPlacementValue& sp) {
    if (&sp != this) {
        delete rep; rep=0;
        if (sp.rep) { rep = sp.rep->clone(); rep->setMyHandle(*this); }
    }
    return *this;
}

} // namespace simtk


#endif // SIMTK_FEATURE_REP_H_
