#ifndef SIMTK_SUBSYSTEM_REP_H_
#define SIMTK_SUBSYSTEM_REP_H_

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
 * Declaration of the library side (rep) of the client Subsystem class. Rebuilding this
 * code does not require recompiling client-side code, and if the library is a DLL
 * (shared library) then the client doesn't require relinking either.
 */

#include "SimbodyCommon.h"
#include "Subsystem.h"
#include "Feature.h"
#include "Placement.h"
#include "PlacementValue.h"
#include "PlacementRep.h"
#include "PlacementValueRep.h"

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
    virtual void finalizeStandardSubfeatures() { }

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

    int getNSubsystems()       const {return childSubsystems.size();}
    int getNPlacements()       const {return placementSlots.size();}
    int getNPlacementValues()  const {return placementValueSlots.size();}

    const Subsystem&      getSubsystem(size_t i)           const {return childSubsystems[i];}
    Subsystem&            updSubsystem(size_t i)                 {return childSubsystems[i];}

    const Placement&  getPlacementSlot(size_t i) const {return placementSlots[i];}
    Placement&        updPlacementSlot(size_t i)       {return placementSlots[i];}

    // PlacementValues here are 'mutable', so the upd routine is const.
    const PlacementValueSlot& getPlacementValueSlot(size_t i) const {return placementValueSlots[i];}
    PlacementValueSlot&       updPlacementValueSlot(size_t i) const {return placementValueSlots[i];}

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
    PlacementValueSlot& addPlacementValueLike(const PlacementValue& v) const;

    const Subsystem& findRootSubsystem() const;
    Subsystem&       findUpdRootSubsystem();      

    // Wipe out all PlacementValues. Note that this routine is const because PlacementValues
    // are mutable.
    void deleteAllPlacementValues() const {
        // Erase all references to the placement values (these may be in child subsystems).
        for (size_t i=0; i < placementValueSlots.size(); ++i) {
            PlacementValueSlot& pv = placementValueSlots[i];
            if (pv.hasClientPlacement())
                pv.getClientPlacement().getRep().clearValueSlot();
        }
        placementValueSlots.clear();
    }

    // Run around this subsystem and its children looking for Placements which can be
    // evaluated without going any higher in the Subsystem tree than the target Subsysytem, and
    // allocate PlacementValueSlots for them in the target.
    void allocatePlacementValueSlots(const Subsystem& target) const {
        // Always let the children go first.
        for (int i=0; i < getNSubsystems(); ++i)
            getSubsystem(i).getRep().allocatePlacementValueSlots(target);

        for (int i=0; i < getNPlacements(); ++i) {
            const Placement& ps = getPlacementSlot(i);
            assert(ps.hasOwner());
            if (ps.getRep().hasValueSlot())
                continue;
            // can we get it a slot?
            const Subsystem* s = ps.getRep().findPlacementValueOwnerSubsystem(target);
            if (s && s->isSameSubsystem(target)) {
                PlacementValueSlot& pvs = 
                    s->getRep().addPlacementValueLike(ps.getRep().createEmptyPlacementValue());
                pvs.setClientPlacement(ps);
                ps.getRep().assignValueSlot(pvs);
            }
        }
    }

    // Is Subsystem f in the tree rooted at oldRoot? If so, optionally return the 
    // series of indices required to get to this Subsystem from the root.
    // Complexity is O(log n) where n is tree depth.
    static bool isSubsystemInSubsystemTree(const Subsystem& oldRoot, const Subsystem& f,
                                           std::vector<int>* trace=0);

    // Is PlacementSlot p owned by a Subsystem in the tree rooted at oldRoot?
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

    // If PlacementSlot p's owner Subsystem is a member of the Subsystem tree rooted at oldRoot,
    // find the corresponding PlacementSlot in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Placement* findCorrespondingPlacementSlot
        (const Subsystem& oldRoot, const Placement& p, const Subsystem& newRoot);

    // If PlacementValueSlot v's owner Subsystem is a member of the Subsystem tree rooted at oldRoot,
    // find the corresponding PlacementValueSlot in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const PlacementValueSlot* findCorrespondingPlacementValueSlot
        (const Subsystem& oldRoot, const PlacementValueSlot& v, const Subsystem& newRoot);

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
    StableArray<SubPlacement>   placementSlots;

    // This is like a State cache except it holds values for Placement expressions
    // where the highest placement dependency is resolved at this Feature (i.e., is
    // one of the placementExpressions stored above. But note that these values
    // do not in general correspond to the placementExpression; they can be the
    // values of lower-level placement expressions.
    mutable StableArray<PlacementValueSlot> placementValueSlots;
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


#endif // SIMTK_SUBSYSTEM_REP_H_
