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

#include "Feature.h"
#include "Placement.h"
#include "PlacementRep.h"

#include <string>
#include <vector>
#include <cassert>
#include <sstream>
#include <cctype>

// Returns -1, 0, 1 according to key {<,==,>} test ignoring case.
static int caseInsensitiveCompare(const std::string& key, const std::string& test);


namespace simtk {

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
      : handle(&f), name(nm), parent(0), indexInParent(-1), placement(0){ }
    
    virtual ~FeatureRep() { 
        for (size_t i=0; i < childFeatures.size(); ++i)
            delete childFeatures[i];
        for (size_t i=0; i < placementExpressions.size(); ++i)
            delete placementExpressions[i];
    }

    virtual PlacementType getRequiredPlacementType() const = 0;
    virtual std::string getFeatureTypeName() const = 0;
    virtual FeatureRep* cloneWithoutPlacement(Feature&) const = 0;

    void setName(const std::string& nm)         {name=nm;}
    void setParent(const Feature& p, int ix)    {parent = &p; indexInParent=ix;}

    const std::string& getName()   const { return name; }
    const Feature&     getHandle() const { return *handle; }

    bool           hasParentFeature() const {return parent != 0;}
    const Feature& getParentFeature() const {assert(hasParentFeature()); return *parent;}
    int            getIndexInParent() const {assert(hasParentFeature()); return indexInParent;}

    bool hasPlacement() const { return placement != 0; }
    void setPlacement(const Placement& p) { placement = &p; }
    const Placement& getPlacement() const { return *placement; }

    const std::vector<Feature*>&   getChildFeatures() const {return childFeatures; }
    const std::vector<Placement*>& getPlacementExpressions() const {return placementExpressions;}

    const Feature& getChildFeature(size_t i) const {
        assert(childFeatures[i]);
        return *childFeatures[i];
    }

    const Placement& getPlacementExpression(size_t i) const {
        assert(placementExpressions[i]);
        return *placementExpressions[i];
    }

    const Feature* getChildFeature(const std::string& nm) const {
        size_t index;
        return findChildFeatureIndex(nm,index) ? childFeatures[index] : 0;
    }
    Feature* updChildFeature(const std::string& nm) {
        size_t index;
        return findChildFeatureIndex(nm,index) ? childFeatures[index] : 0;
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
        childFeatures.push_back(new Feature(f)); // a copy of f but without external placements
        Feature& newFeature = *childFeatures[index];
        updRep(newFeature)->setParent(*handle, index);
        updRep(newFeature)->setName(nm);
        return newFeature;
    }

    // TODO: this should only allow placements involving this feature, its children,
    // grandchildren, etc.
    Placement& addPlacementLike(const Placement& p) {
        assert(PlacementRep::getRep(p));
        const int index = (int)placementExpressions.size();
        placementExpressions.push_back(new Placement(p)); // a copy of p without an owner
        Placement& newPlacement = *placementExpressions[index];
        PlacementRep::updRep(newPlacement)->setOwner(*handle, index);
        return newPlacement;
    }

    SIMTK_REP_HELPERS(Feature,FeatureRep)
protected:
    // Clones of FeatureReps correspond to a new handle Feature, and they
    // lose their associations with the upstream world (parent & ancestors).
    // That means all Placement references which are not owned by this 
    // feature or its children are removed.
    void cleanUpAfterClone(Feature& f) {
        handle = &f;
        parent=0;
        indexInParent = -1;
        removeExternalPlacements();
    }
private:
    // Return true and ix==feature index if a feature of the given name is found.
    // Otherwise return false and ix==childFeatures.size().
    bool findChildFeatureIndex(const std::string& nm, size_t& ix) const {
        for (ix=0; ix < getChildFeatures().size(); ++ix)
            if (caseInsensitiveCompare(nm, childFeatures[ix]->getName())==0)
                return true;
        return false;   // not found
    }

    void removeExternalPlacements() {
        std::cout << "removeExternalPlacement() needed but not written yet!" << std::endl;
    }
private:
    Feature*                handle; // the Feature whose rep this is
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
    std::vector<Feature*>    childFeatures;

    // Placement expressions wholly owned by this Feature. These compute values
    // for the childFeatures of this Feature.
    std::vector<Placement*>  placementExpressions;
};

class RealParameterRep : public FeatureRep {
public:
    RealParameterRep(RealParameter& p, const std::string& nm) : FeatureRep(p,nm) { }

    std::string getFeatureTypeName() const { return "RealParameter"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        RealParameterRep* copy = new RealParameterRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(RealParameterRep,FeatureRep);
};

class StationParameterRep : public FeatureRep {
public:
    StationParameterRep(StationParameter& p, const std::string& nm) : FeatureRep(p,nm) { }

    std::string getFeatureTypeName() const { return "StationParameter"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        StationParameterRep* copy = new StationParameterRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(StationParameterRep,FeatureRep);
};

class RealMeasureRep : public FeatureRep {
public:
    RealMeasureRep(RealMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }

    std::string getFeatureTypeName() const { return "RealMeasure"; }
    PlacementType getRequiredPlacementType() const { return RealPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        RealMeasureRep* copy = new RealMeasureRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(RealMeasureRep,FeatureRep);
};

class StationMeasureRep : public FeatureRep {
public:
    StationMeasureRep(StationMeasure& m, const std::string& nm) : FeatureRep(m,nm) { }

    std::string getFeatureTypeName() const { return "StationMeasure"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        StationMeasureRep* copy = new StationMeasureRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(StationMeasureRep,FeatureRep);
};

class StationRep : public FeatureRep {
public:
    StationRep(Station& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Station"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        StationRep* copy = new StationRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(StationRep,FeatureRep);
};

class DirectionRep : public FeatureRep {
public:
    DirectionRep(Direction& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Direction"; }
    PlacementType getRequiredPlacementType() const { return DirectionPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        DirectionRep* copy = new DirectionRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(DirectionRep,FeatureRep);
};

class OrientationRep : public FeatureRep {
public:
    OrientationRep(Orientation& s, const std::string& nm) : FeatureRep(s,nm) { }

    std::string getFeatureTypeName() const { return "Orientation"; }
    PlacementType getRequiredPlacementType() const { return OrientationPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        OrientationRep* copy = new OrientationRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST(OrientationRep,FeatureRep);
};

class FrameRep : public FeatureRep {
public:
    FrameRep(Frame& f, const std::string& nm) : FeatureRep(f,nm) {
        initializeFeatures();
    }

    std::string getFeatureTypeName() const { return "Frame"; }
    PlacementType getRequiredPlacementType() const { return FramePlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        FrameRep* copy = new FrameRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    const Station& getOrigin() const { return Station::downcast(*getChildFeatures()[3]); }
    const Direction& getAxis(int i) const 
      { assert(0<=i && i<=2); return Direction::downcast(*getChildFeatures()[i]); }
    const Direction& x() const { return Direction::downcast(*getChildFeatures()[0]); }
    const Direction& y() const { return Direction::downcast(*getChildFeatures()[1]); }
    const Direction& z() const { return Direction::downcast(*getChildFeatures()[2]); }

    SIMTK_DOWNCAST(FrameRep,FeatureRep);
private:
    void initializeFeatures() {
        Station&   O = Station::downcast  (addFeatureLike(Station("O"), "O"));
        Direction& x = Direction::downcast(addFeatureLike(Direction("x"), "x"));
        Direction& y = Direction::downcast(addFeatureLike(Direction("y"), "y"));
        Direction& z = Direction::downcast(addFeatureLike(Direction("z"), "z"));

        O.setPlacement(addPlacementLike(StationPlacement  (Vec3(0,0,0))));
        x.setPlacement(addPlacementLike(DirectionPlacement(Vec3(1,0,0))));
        y.setPlacement(addPlacementLike(DirectionPlacement(Vec3(0,1,0))));
        z.setPlacement(addPlacementLike(DirectionPlacement(Vec3(0,0,1))));
    }
};



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
