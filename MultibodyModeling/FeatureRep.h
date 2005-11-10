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

// Returns -1, 0, 1 according to key {<,==,>} test ignoring case.
static int caseInsensitiveCompare(const std::string& key, const std::string& test);


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
    { f.rep = this; }

    // Copying a Feature is tricky. The result should have all the child features
    // and the *internal* placements. External placements should evaporate. Note
    // that the index numbers for features & placements we own must stay the
    // same so that internal references in the copy are the same as in the original.
    // However, this is all handled in the Feature copy & assignment methods --
    // FetureRep copying is elementwise and dumb and thus dangerous. The idea is
    // to get a straight copy and then go clean up the mess afterwards.

    
    virtual ~FeatureRep() { }

    void           setMyHandle(Feature& f) {myHandle = &f;}
    const Feature& getMyHandle() const     {return *myHandle;}

    virtual PlacementType getRequiredPlacementType()      const = 0;
    virtual std::string   getFeatureTypeName()            const = 0;
    virtual FeatureRep*   clone()                         const = 0;

    void cloneWithoutExternalPlacements(Feature& newHandle) const {
        FeatureRep* copy = clone();
        copy->setMyHandle(newHandle);
        FeatureRep::setRep(newHandle, copy);

        // Re-parent all the copied child Features to their new parent,
        // and fix the owned Placements to acknoledge their new owner.
        copy->reparentMyChildren();

        // Fix up all the internal placement references and delete the
        // external ones.
        copy->fixPlacements(this->getMyHandle(), copy->getMyHandle());
    }

    void               setName(const std::string& nm) {name = nm;}
    const std::string& getName() const                {return name;}

    void           setParentFeature(const Feature& p, int ix) {parent = &p; indexInParent=ix;}
    bool           hasParentFeature() const {return parent != 0;}
    const Feature& getParentFeature() const {assert(hasParentFeature()); return *parent;}
    int            getIndexInParent() const {assert(hasParentFeature()); return indexInParent;}

    bool hasPlacement() const { return placement != 0; }
    void setPlacement(const Placement& p) { 
        assert(PlacementRep::getRep(p)->getPlacementType() == getRequiredPlacementType());
        placement = &p; 
    }
    const Placement& getPlacement() const { return *placement; }

    int getNChildFeatures()        const {return childFeatures.size();}
    int getNPlacementExpressions() const {return placementExpressions.size();}

    const Feature&   getChildFeature(size_t i)        const {return childFeatures[i];}
    const Placement& getPlacementExpression(size_t i) const {return placementExpressions[i];}

    const Feature* getChildFeature(const std::string& nm) const {
        size_t index;
        return findChildFeatureIndex(nm,index) ? &childFeatures[index] : 0;
    }
    Feature* updChildFeature(const std::string& nm) {
        size_t index;
        return findChildFeatureIndex(nm,index) ? &childFeatures[index] : 0;
    }

    const std::string getFullName() const { 
        std::string s;
        if (hasParentFeature())
            s = getParentFeature().getFullName() + ".";
        return s + getName(); 
    }

    Feature& addFeatureLike(const Feature& f, const std::string& nm) {
        assert(getRep(f) && nm.size() > 0);
        const int index = (int)childFeatures.size();
        childFeatures.push_back(SubFeature()); // an empty handle
        Feature& newFeature = childFeatures[index];
        FeatureRep::getRep(f)->cloneWithoutExternalPlacements(newFeature);
        FeatureRep::updRep(newFeature)->setParentFeature(getMyHandle(), index);
        FeatureRep::updRep(newFeature)->setName(nm);
        return newFeature;
    }

    // TODO: this should only allow placements involving this feature, its children,
    // grandchildren, etc.
    Placement& addPlacementLike(const Placement& p) {
        assert(PlacementRep::getRep(p));

        Feature* offender;
        if (!PlacementRep::getRep(p)->isLimitedToSubtree(getMyHandle(),offender)) {
            SIMTK_THROW3(Exception::PlacementMustBeLocal,"FeatureRep::addPlacementLike",
                this->getFullName(),offender->getFullName());
        }

        const int index = (int)placementExpressions.size();
        placementExpressions.push_back(Placement());
        Placement& newPlacement = placementExpressions[index];
        PlacementRep::getRep(p)->clone(newPlacement);
        PlacementRep::updRep(newPlacement)->setOwner(getMyHandle(), index);
        return newPlacement;
    }

    
    // Is Feature f in the tree rooted at oldRoot? If so, optionally return the 
    // series of indices required to get to this Feature from the root.
    static bool isFeatureInFeatureTree(const Feature& oldRoot, const Feature& f,
                                       std::vector<int>* trace=0)
    {
        if (trace) trace->clear();
        const Feature* const oldp = &oldRoot;
        const Feature*       fp   = &f;
        for(;;) {
            if (fp == oldp) return true;
            if (!fp->rep->hasParentFeature())
                break;
            if (trace) trace->push_back(fp->rep->getIndexInParent());
            fp = &fp->rep->getParentFeature();
        }
        if (trace) trace->clear(); // never mind ...
        return false;
    }

    // Is Placement p owned by a Feature in the tree rooted at oldRoot? If so, 
    // optionally return the series of indices required to get to this Placement's
    // owner Feature from the root.
    static bool isPlacementInFeatureTree(const Feature& oldRoot, const Placement& p)
    {
        if (!p.hasOwner())
            return false;   // a disembodied Placement

        return isFeatureInFeatureTree(oldRoot, p.getOwner());
    }

    SIMTK_REP_HELPERS(Feature,FeatureRep)
private:
    // Return true and ix==feature index if a feature of the given name is found.
    // Otherwise return false and ix==childFeatures.size().
    bool findChildFeatureIndex(const std::string& nm, size_t& ix) const {
        for (ix=0; ix < (size_t)getNChildFeatures(); ++ix)
            if (caseInsensitiveCompare(nm, childFeatures[ix].getName())==0)
                return true;
        return false;   // not found
    }

    // We have just copied a Feature subtree so all the parent pointers are
    // still pointing to the old tree. Recursively repair them to point into
    // the new tree.
    void reparentMyChildren() {
        for (size_t i=0; i < (size_t)getNChildFeatures(); ++i) {
            assert(childFeatures[i].rep->getIndexInParent() == i); // shouldn't change
            childFeatures[i].rep->setParentFeature(getMyHandle(), i);
            childFeatures[i].rep->reparentMyChildren();            // recurse
        }
        for (size_t i=0; i < (size_t)getNPlacementExpressions(); ++i) {
            assert(placementExpressions[i].getIndexInOwner() == i); // shouldn't change
            PlacementRep::updRep(placementExpressions[i])->setOwner(getMyHandle(), i);
        }
    }

    // We have just created at newRoot a copy of the tree rooted at oldRoot, and the
    // current Feature (for which this is the Rep) is a node in the newRoot tree
    // (with correct myHandle). However, the 'parent' and 'placement' pointers 
    // still retain the values they had in the oldRoot tree; they must be 
    // changed to point to the corresponding entities in the newRoot tree.
    // If these pointers point outside the oldRoot tree, however, we'll just
    // set them to 0 in the newRoot copy.
    void fixPlacements(const Feature& oldRoot, const Feature& newRoot) {
        for (size_t i=0; i < (size_t)getNChildFeatures(); ++i)
            childFeatures[i].rep->fixPlacements(oldRoot, newRoot);

        if (placement)
            placement = findCorrespondingPlacement(oldRoot,*placement,newRoot);
    }



    // If Feature f is a member of the Feature tree rooted at oldRoot, find
    // the corresponding Feature in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Feature* findCorrespondingFeature
        (const Feature& oldRoot, const Feature& f, const Feature& newRoot)
    {
        std::vector<int> trace;
        if (!isFeatureInFeatureTree(oldRoot,f,&trace))
            return 0;

        // Trace holds the indices needed to step from newRoot down to
        // the corresponding Feature (in reverse order).
        const Feature* newTreeRef = &newRoot;
        for (size_t i=trace.size(); i >=1; --i)
            newTreeRef = &newTreeRef->rep->getChildFeature(trace[i-1]);
        return newTreeRef;
    }

    
    // If Placement p's owner Feature is a member of the Feature tree rooted at oldRoot,
    // find the corresponding Placement in the tree rooted at newRoot (which is expected
    // to be a copy of oldRoot). Return NULL if not found for any reason.
    static const Placement* findCorrespondingPlacement
        (const Feature& oldRoot, const Placement& p, const Feature& newRoot)
    {
        if (!p.hasOwner()) return 0;
        const Feature* corrOwner = findCorrespondingFeature(oldRoot,p.getOwner(),newRoot);
        if (!corrOwner) return 0;
        assert(corrOwner->rep);

        const Placement* newTreeRef = 
            &corrOwner->rep->getPlacementExpression(p.getIndexInOwner());
        assert(newTreeRef);
        assert(&newTreeRef->getOwner() == corrOwner);
        assert(newTreeRef->getIndexInOwner() == p.getIndexInOwner());
        return newTreeRef;
    }
private:
    Feature*                myHandle; // the Feature whose rep this is
    std::string             name;

    // Owner information. If parent is 0 then this Feature is an unplaced
    // prototype. Otherwise this Feature is contained in the parent's
    // childFeatures list, with the given index.
    const Feature*          parent;
    int                     indexInParent;

    // If this Feature has been placed, this is the placement information.
    // If present, this Placement must be owned by this Feature, its parent,
    // or one of its ancestors.
    const Placement*        placement;

    // Subfeatures wholly owned by this Feature.
    StableArray<SubFeature> childFeatures;

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

    SIMTK_DOWNCAST(RealParameterRep,FeatureRep);
};

class StationParameterRep : public FeatureRep {
public:
    StationParameterRep(StationParameter& p, const std::string& nm) : FeatureRep(p,nm) { }

    std::string getFeatureTypeName() const { return "StationParameter"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationParameterRep(*this); }

    SIMTK_DOWNCAST(StationParameterRep,FeatureRep);
};

class RealMeasureRep : public FeatureRep {
public:
    RealMeasureRep(RealMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }

    std::string getFeatureTypeName() const { return "RealMeasure"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* clone() const { return new RealMeasureRep(*this); }

    SIMTK_DOWNCAST(RealMeasureRep,FeatureRep);
};

class StationMeasureRep : public FeatureRep {
public:
    StationMeasureRep(StationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }

    std::string getFeatureTypeName() const { return "StationMeasure"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationMeasureRep(*this); }

    SIMTK_DOWNCAST(StationMeasureRep,FeatureRep);
};

class StationRep : public FeatureRep {
public:
    StationRep(Station& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Station"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new StationRep(*this); }

    SIMTK_DOWNCAST(StationRep,FeatureRep);
};

class DirectionRep : public FeatureRep {
public:
    DirectionRep(Direction& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Direction"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    FeatureRep* clone() const { return new DirectionRep(*this); }

    SIMTK_DOWNCAST(DirectionRep,FeatureRep);
};

class OrientationRep : public FeatureRep {
public:
    OrientationRep(Orientation& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Orientation"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    FeatureRep* clone() const { return new OrientationRep(*this); }

    SIMTK_DOWNCAST(OrientationRep,FeatureRep);
};

class FrameRep : public FeatureRep {
public:
    FrameRep(Frame& f, const std::string& nm) : FeatureRep(f,nm) {
        initializeFeatures();
    }

    std::string getFeatureTypeName() const { return "Frame"; }
    PlacementType getRequiredPlacementType() const { return FramePlacementType; }
    FeatureRep* clone() const { return new FrameRep(*this); }

    const Orientation& getOrientation() const {return Orientation::downcast(getChildFeature(0));}
    const Station&     getOrigin()      const {return Station::downcast(getChildFeature(1)); }

    SIMTK_DOWNCAST(FrameRep,FeatureRep);
private:
    void initializeFeatures() {
        Station&     O = Station::downcast    (addFeatureLike(Station("O"),     "origin"));
        Orientation& R = Orientation::downcast(addFeatureLike(Orientation("R"), "orientation"));

        O.setPlacement(addPlacementLike(StationPlacement  (Vec3(0))));
        R.setPlacement(addPlacementLike(OrientationPlacement(Mat33(1))));
    }
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

static int caseInsensitiveCompare(const std::string& key, const std::string& test) {
    const size_t minlen = std::min(key.size(), test.size());
    for (size_t i=0; i < minlen; ++i) {
        const int k = tolower(key[i]), t = tolower(test[i]);
        if (k < t) return -1;
        else if (k > t) return 1;
    }
    // caution -- size() is unsigned, don't get clever here
    if (key.size() > minlen) return 1;
    else if (test.size() > minlen) return -1;
    return 0;
}

#endif // SIMTK_FEATURE_REP_H_
