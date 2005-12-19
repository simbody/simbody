#ifndef SIMTK_BASIC_FEATURES_REP_H_
#define SIMTK_BASIC_FEATURES_REP_H_

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
 * Definitions of the BasicFeature Rep methods.
 */

#include "SimbodyCommon.h"
#include "FeatureRep.h"
#include "Feature.h"
#include "Placement.h"
#include "PlacementRep.h"

#include <string>
#include <cassert>
#include <sstream>
#include <cctype>

namespace simtk {

class RealParameterRep : public FeatureRep {
public:
    RealParameterRep(RealParameter& p, const std::string& nm) 
        : FeatureRep(p,nm,RealPlacement(NTraits<Real>::getNaN())) { }
    // no standard Subfeatures

    ~RealParameterRep() { }

    std::string getFeatureTypeName() const { return "RealParameter"; }

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


    PlacementRep* useFeatureAsRealPlacement(RealPlacement& handle) const {
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(RealParameterRep,SubsystemRep);
};

class Vec3ParameterRep : public FeatureRep {
public:
    Vec3ParameterRep(Vec3Parameter& p, const std::string& nm) 
        : FeatureRep(p,nm,Vec3Placement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~Vec3ParameterRep() { }

    std::string getFeatureTypeName() const { return "Vec3Parameter"; }
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

    PlacementRep* useFeatureAsVec3Placement(Vec3Placement& handle) const {
        PlacementRep* prep = new Vec3FeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(Vec3ParameterRep,SubsystemRep);
};

class StationParameterRep : public FeatureRep {
public:
    StationParameterRep(StationParameter& p, const std::string& nm) 
        : FeatureRep(p,nm,StationPlacement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~StationParameterRep() { }

    std::string getFeatureTypeName() const { return "StationParameter"; }
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

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationParameterRep,SubsystemRep);
};

class RealMeasureRep : public FeatureRep {
public:
    RealMeasureRep(RealMeasure& m, const std::string& nm) 
        : FeatureRep(m,nm,RealPlacement(NTraits<Real>::getNaN())) { }
    // no standard Subfeatures

    ~RealMeasureRep() { }

    std::string getFeatureTypeName() const { return "RealMeasure"; }
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

    PlacementRep* useFeatureAsRealPlacement(RealPlacement& handle) const {
        PlacementRep* prep = new RealFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(RealMeasureRep,SubsystemRep);
};

class Vec3MeasureRep : public FeatureRep {
public:
    Vec3MeasureRep(Vec3Measure& m, const std::string& nm) 
        : FeatureRep(m,nm,Vec3Placement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~Vec3MeasureRep() { }

    std::string getFeatureTypeName() const { return "Vec3Measure"; }
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

    PlacementRep* useFeatureAsVec3Placement(Vec3Placement& handle) const {
        PlacementRep* prep = new Vec3FeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(Vec3MeasureRep,SubsystemRep);
};

class StationMeasureRep : public FeatureRep {
public:
    StationMeasureRep(StationMeasure& m, const std::string& nm) 
        : FeatureRep(m,nm,StationPlacement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~StationMeasureRep() { }

    std::string getFeatureTypeName() const { return "StationMeasure"; }
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

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationMeasureRep,SubsystemRep);
};

class StationRep : public FeatureRep {
public:
    StationRep(Station& s, const std::string& nm) 
        : FeatureRep(s,nm,StationPlacement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~StationRep() { }

    std::string getFeatureTypeName() const { return "Station"; }
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

    PlacementRep* useFeatureAsStationPlacement(StationPlacement& handle) const {
        PlacementRep* prep = new StationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentSubsystem() && FrameFeature::isInstanceOf(getParentSubsystem()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Station", "Orientation");
            //NOTREACHED
        }
        const FrameFeature& parentFrame = FrameFeature::downcast(getParentSubsystem());
        PlacementRep* prep = new FrameExprPlacementRep(parentFrame.getOrientation(), 
                                                       Station::downcast(getMyHandle()));
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(StationRep,SubsystemRep);
};

class DirectionMeasureRep : public FeatureRep {
public:
    DirectionMeasureRep(DirectionMeasure& m, const std::string& nm) 
        : FeatureRep(m,nm,DirectionPlacement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~DirectionMeasureRep() { }

    std::string getFeatureTypeName() const { return "DirectionMeasure"; }
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

    PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement& handle) const {
        PlacementRep* prep = new DirectionFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(DirectionMeasureRep,SubsystemRep);
};

class DirectionRep : public FeatureRep {
public:
    DirectionRep(Direction& s, const std::string& nm) 
        : FeatureRep(s,nm,DirectionPlacement(Vec3(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~DirectionRep() { }

    std::string getFeatureTypeName() const { return "Direction"; }
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

    PlacementRep* useFeatureAsDirectionPlacement(DirectionPlacement& handle) const {
        PlacementRep* prep = new DirectionFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(DirectionRep,SubsystemRep);
};


class OrientationMeasureRep : public FeatureRep {
public:
    OrientationMeasureRep(OrientationMeasure& m, const std::string& nm) 
        : FeatureRep(m,nm,OrientationPlacement(Mat33(NTraits<Real>::getNaN()))) { }
    // no standard Subfeatures

    ~OrientationMeasureRep() { }

    std::string getFeatureTypeName() const { return "OrientationMeasure"; }
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

    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }

    SIMTK_DOWNCAST(OrientationMeasureRep,SubsystemRep);
};

class OrientationRep : public FeatureRep {
public:
    OrientationRep(Orientation& o, const std::string& nm) 
        : FeatureRep(o,nm,OrientationPlacement(Mat33(NTraits<Real>::getNaN())))
      { axisIndices[0]=axisIndices[1]=axisIndices[2] = -1; }
    // must call initializeStandardSubfeatures() to complete construction.

    ~OrientationRep() { }

    std::string   getFeatureTypeName() const { return "Orientation"; }
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

    PlacementRep* useFeatureAsOrientationPlacement(OrientationPlacement& handle) const {
        PlacementRep* prep = new OrientationFeaturePlacementRep(getMyHandle());
        prep->setMyHandle(handle); handle.setRep(prep);
        return prep;
    }
    PlacementRep* useFeatureAsFramePlacement(FramePlacement& handle) const {
        if (!(hasParentSubsystem() && FrameFeature::isInstanceOf(getParentSubsystem()))) {
            SIMTK_THROW3(Exception::FeatureUsedAsFramePlacementMustBeOnFrame,
                     getFullName(), "Orientation", "Station");
            //NOTREACHED
        }
        const FrameFeature& parentFrame = FrameFeature::downcast(getParentSubsystem());
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
    FrameRep(FrameFeature& f, const std::string& nm) 
      : FeatureRep(f,nm,FramePlacement(Mat34(NTraits<Real>::getNaN()))), RIndex(-1), OIndex(-1) { }
    // must call initializeStandardSubfeatures() to complete construction.

    ~FrameRep() { }

    // still overrideable for bodies.
    virtual std::string   getFeatureTypeName() const { return "FrameFeature"; }
    virtual SubsystemRep*   clone() const { return new FrameRep(*this); }

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


#endif // SIMTK_BASIC_FEATURES_REP_H_
