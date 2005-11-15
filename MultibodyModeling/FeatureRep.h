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
      { f.setRep(this); }

    // Copying a Feature is tricky. The result should have all the child features
    // and the *internal* placements. External placements should evaporate. Note
    // that the index numbers for features & placements we own must stay the
    // same so that internal references in the copy are the same as in the original.
    // However, this is all handled in the Feature copy & assignment methods --
    // FetureRep copying is elementwise and dumb and thus dangerous. The idea is
    // to get a straight copy and then go clean up the mess afterwards.

    
    virtual ~FeatureRep() { }

    void           setMyHandle(Feature& f) {myHandle = &f;}
    const Feature& getMyHandle() const     {assert(myHandle); return *myHandle;}
    Feature&       updMyHandle()           {assert(myHandle); return *myHandle;}

    virtual PlacementType getRequiredPlacementType()      const = 0;
    virtual std::string   getFeatureTypeName()            const = 0;
    virtual FeatureRep*   clone()                         const = 0;

    // If this Feature can be used as the indicated placement type, return
    // a new, unowned placement of the right type. Most commonly, the returned
    // Placement will just be a recast FeaturePlacement referencing this Feature.
    // For composite Features, this can be a recast FeaturePlacement
    // referencing one of the Features's subfeatures.
    // For example, if a Frame is used as a StationPlacement, we return a
    // reference to the Frame's origin feature. Otherwise we return invalid
    // placements (handles with no rep).
    virtual void useAsRealPlacement(RealPlacement&) const {
        SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                     getFullName(), getFeatureTypeName(), "Real");
    }
    virtual void useAsVec3Placement(Vec3Placement&) const {
        SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                     getFullName(), getFeatureTypeName(), "Vec3");
    }
    virtual void useAsStationPlacement(StationPlacement&) const {
        SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                     getFullName(), getFeatureTypeName(), "Station");
    }
    virtual void useAsDirectionPlacement(DirectionPlacement&) const {
        SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                     getFullName(), getFeatureTypeName(), "Direction");
    }    
    virtual void useAsOrientationPlacement(OrientationPlacement&) const {
        SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                     getFullName(), getFeatureTypeName(), "Orientation");
    } 
    virtual void useAsFramePlacement(FramePlacement&) const {
        SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                     getFullName(), getFeatureTypeName(), "Frame");
    } 

    void cloneWithoutExternalPlacements(Feature& newHandle) const;

    void               setName(const std::string& nm) {name = nm;}
    const std::string& getName() const                {return name;}

    void           setParentFeature(Feature& p, int ix) {parent = &p; indexInParent=ix;}
    bool           hasParentFeature() const {return parent != 0;}
    const Feature& getParentFeature() const {assert(hasParentFeature()); return *parent;}
    Feature&       updParentFeature() const {assert(hasParentFeature()); return *parent;}
    int            getIndexInParent() const {assert(hasParentFeature()); return indexInParent;}

    bool hasPlacement() const { return placement != 0; }
    const Placement& getPlacement() const { return *placement; }
    void place(const Placement& p);


    int getNSubfeatures()          const {return subfeatures.size();}
    int getNPlacementExpressions() const {return placementExpressions.size();}

    const Feature&   getSubfeature(size_t i)          const {return subfeatures[i];}
    Feature&         updSubfeature(size_t i)                {return subfeatures[i];}
    const Placement& getPlacementExpression(size_t i) const {return placementExpressions[i];}

    const Feature* findSubfeature(const std::string& nm) const {
        size_t index;
        return findSubfeatureIndex(nm,index) ? &subfeatures[index] : 0;
    }
    Feature* findUpdSubfeature(const std::string& nm) {
        size_t index;
        return findSubfeatureIndex(nm,index) ? &subfeatures[index] : 0;
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

    Feature&   addSubfeatureLike(const Feature& f, const std::string& nm);
    Placement& addPlacementLike(const Placement& p);

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

    // Is Placement p owned by a Feature in the tree rooted at oldRoot? If so, 
    // optionally return the series of indices required to get to this Placement's
    // owner Feature from the root.
    static bool isPlacementInFeatureTree(const Feature& oldRoot, const Placement& p);

    // If Feature f is a member of the Feature tree rooted at oldRoot, find
    // the corresponding Feature in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Feature* findCorrespondingFeature
        (const Feature& oldRoot, const Feature& f, const Feature& newRoot);

    // Given two features, run up the tree towards the root to find
    // their "least common denominator", i.e. the first shared node
    // on the path back to the root. Return a pointer to that node
    // if found, otherwise NULL meaning that the features aren't on
    // the same tree. If the features are the same, then
    // that feature is the answer.
    // Complexity is O(log n) (3 passes) where n is depth of Feature tree.
    static const Feature* findYoungestCommonAncestor(const Feature& f1, const Feature& f2);
    static Feature* findUpdYoungestCommonAncestor(Feature& f1, const Feature& f2);

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

    // If Placement p's owner Feature is a member of the Feature tree rooted at oldRoot,
    // find the corresponding Placement in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Placement* findCorrespondingPlacement
        (const Feature& oldRoot, const Placement& p, const Feature& newRoot);
    
    static const Feature* getParentPtr(const Feature& f) {
        return f.rep ? f.rep->parent : 0;
    }
private:
    Feature*                myHandle; // the Feature whose rep this is
    std::string             name;

    // Owner information. If parent is 0 then this Feature is an unplaced
    // prototype. Otherwise this Feature is contained in the parent's
    // childFeatures list, with the given index.
    Feature*                parent;
    int                     indexInParent;

    // If this Feature has been placed, this is the placement information.
    // If present, this Placement must be owned by this Feature, its parent,
    // or one of its ancestors.
    const Placement*        placement;

    // Subfeatures wholly owned by this Feature.
    StableArray<SubFeature> subfeatures;

    // Placement expressions wholly owned by this Feature. These compute values
    // for the childFeatures of this Feature.
    StableArray<Placement>  placementExpressions;
};

class RealParameterRep : public FeatureRep {
public:
    RealParameterRep(RealParameter& p, const std::string& nm) : FeatureRep(p,nm) { }

    std::string getFeatureTypeName() const { return "RealParameter"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* clone() const { return new RealParameterRep(*this); }

    void useAsRealPlacement(RealPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }

    SIMTK_DOWNCAST(RealParameterRep,FeatureRep);
};

class StationParameterRep : public FeatureRep {
public:
    StationParameterRep(StationParameter& p, const std::string& nm) : FeatureRep(p,nm) { }

    std::string getFeatureTypeName() const { return "StationParameter"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationParameterRep(*this); }

    void useAsStationPlacement(StationPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }

    void useAsVec3Placement(Vec3Placement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }

    SIMTK_DOWNCAST(StationParameterRep,FeatureRep);
};

class RealMeasureRep : public FeatureRep {
public:
    RealMeasureRep(RealMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }

    std::string getFeatureTypeName() const { return "RealMeasure"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* clone() const { return new RealMeasureRep(*this); }

    void useAsRealPlacement(RealPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }

    SIMTK_DOWNCAST(RealMeasureRep,FeatureRep);
};

class StationMeasureRep : public FeatureRep {
public:
    StationMeasureRep(StationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }

    std::string getFeatureTypeName() const { return "StationMeasure"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationMeasureRep(*this); }

    void useAsStationPlacement(StationPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }

    void useAsVec3Placement(Vec3Placement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }

    SIMTK_DOWNCAST(StationMeasureRep,FeatureRep);
};

class StationRep : public FeatureRep {
public:
    StationRep(Station& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Station"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationRep(*this); }

    void useAsStationPlacement(StationPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }
    void useAsVec3Placement(Vec3Placement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }
    void useAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentFeature() && Frame::isInstanceOf(getParentFeature()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Station", "Orientation");
            //NOTREACHED
        }
        const Frame& parentFrame = Frame::downcast(getParentFeature());
        (void)new FramePlacementRep(handle, parentFrame.getOrientation(), 
                                    Station::downcast(getMyHandle()));
    }

    SIMTK_DOWNCAST(StationRep,FeatureRep);
};

class DirectionRep : public FeatureRep {
public:
    DirectionRep(Direction& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Direction"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    FeatureRep* clone() const { return new DirectionRep(*this); }

    void useAsDirectionPlacement(DirectionPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }
    void useAsVec3Placement(Vec3Placement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }
    SIMTK_DOWNCAST(DirectionRep,FeatureRep);
};

class OrientationRep : public FeatureRep {
public:
    OrientationRep(Orientation& o, const std::string& nm) : FeatureRep(o,nm)
        { initializeSubfeatures(); }

    std::string   getFeatureTypeName() const { return "Orientation"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    FeatureRep*   clone() const { return new OrientationRep(*this); }

    void useAsOrientationPlacement(OrientationPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }
    void useAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentFeature() && Frame::isInstanceOf(getParentFeature()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Orientation", "Station");
            //NOTREACHED
        }
        const Frame& parentFrame = Frame::downcast(getParentFeature());
        (void)new FramePlacementRep(handle, Orientation::downcast(getMyHandle()),
                                    parentFrame.getOrigin());
    }
    const Direction& getAxis(int i) const
      { assert(0<=i&&i<=2); return Direction::downcast(getSubfeature(axisIndices[i])); }
    const Direction& x() const {return Direction::downcast(getSubfeature(axisIndices[0]));}
    const Direction& y() const {return Direction::downcast(getSubfeature(axisIndices[1]));}
    const Direction& z() const {return Direction::downcast(getSubfeature(axisIndices[2]));}

    SIMTK_DOWNCAST(OrientationRep,FeatureRep);
private:
    void initializeSubfeatures() {
        Direction& x = Direction::downcast(addSubfeatureLike(Direction("x"), "x"));
        Direction& y = Direction::downcast(addSubfeatureLike(Direction("y"), "y"));
        Direction& z = Direction::downcast(addSubfeatureLike(Direction("z"), "z"));

        axisIndices[0] = x.getIndexInParent();
        axisIndices[1] = y.getIndexInParent();
        axisIndices[2] = z.getIndexInParent();

        for (int i=0; i<3; ++i)
            updSubfeature(axisIndices[i]).place(FeaturePlacement(getMyHandle(), i));
    }

    int axisIndices[3];
};

class FrameRep : public FeatureRep {
public:
    FrameRep(Frame& f, const std::string& nm) 
      : FeatureRep(f,nm), RIndex(-1), OIndex(-1) {
        initializeSubfeatures();
    }

    std::string   getFeatureTypeName() const { return "Frame"; }
    PlacementType getRequiredPlacementType() const { return FramePlacementType; }
    FeatureRep*   clone() const { return new FrameRep(*this); }

    void useAsFramePlacement(FramePlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getMyHandle());
    }
    void useAsOrientationPlacement(OrientationPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getOrientation());
    }
    void useAsStationPlacement(StationPlacement& handle) const {
        (void)new FeaturePlacementRep(
                    reinterpret_cast<FeaturePlacement&>(handle),
                    getOrigin());
    }
    const Orientation& getOrientation() const {return Orientation::downcast(getSubfeature(RIndex));}
    const Station&     getOrigin()      const {return Station::downcast(getSubfeature(OIndex)); }

    SIMTK_DOWNCAST(FrameRep,FeatureRep);
private:
    void initializeSubfeatures() {
        Orientation& R = Orientation::downcast(addSubfeatureLike(Orientation("R"), "orientation"));
        Station&     O = Station::downcast    (addSubfeatureLike(Station("O"),     "origin"));

        RIndex = R.getIndexInParent();
        OIndex = O.getIndexInParent();

        updSubfeature(RIndex).place(FeaturePlacement(getMyHandle(), 0));
        updSubfeature(OIndex).place(FeaturePlacement(getMyHandle(), 1));
    }

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


} // namespace simtk


#endif // SIMTK_FEATURE_REP_H_
