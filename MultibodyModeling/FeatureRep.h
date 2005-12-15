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
#include "SubsystemRep.h"
#include "Feature.h"
#include "Placement.h"
#include "PlacementRep.h"

#include <string>
#include <cassert>
#include <sstream>
#include <cctype>

namespace simtk {

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

} // namespace simtk


#endif // SIMTK_FEATURE_REP_H_
