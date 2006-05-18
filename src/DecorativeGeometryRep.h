#ifndef SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_
#define SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/AnalyticGeometry.h"

#include <cmath>
#include <vector>

#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkProperty.h"
#include "vtkObject.h"

namespace SimTK {

static const Real Pi = std::acos(Real(-1));

class DecorativeGeometryRep {
public:
    DecorativeGeometryRep() 
      : myHandle(0),  resolution(-1), scale(-1), placement(), 
        colorRGB(-1,-1,-1), opacity(-1), representation(-1)
    { 
    }
    virtual ~DecorativeGeometryRep() {
        deleteVTKGeometry();
        clearMyHandle();
    }

    void generateVTKPipeline() {
        deleteVTKGeometry();
        rememberVTKObject(createVTKPolyData()); // push this on last
    }

    vtkPolyData* getVTKPolyData() {
        assert(vtkObjects.size());
        return vtkPolyData::SafeDownCast(vtkObjects.back());
    }


    // Caller must be sure to call VTK's Delete() methods on these objects
    // when done with them.
    virtual vtkPolyData* createVTKPolyData() = 0;

    // Combine a SimTK Transform with a scale factor and return the equivalent
    // vtkTransform that will translate, rotate, and scale the same way.
    vtkTransform* createVTKTransform(const Transform&, const Real&);

    // Take some polygon data as input, copy it and transform it, returning the
    // transformed polygon data.
    vtkPolyData* transformVTKPolyData(const Transform&, const Real&, vtkPolyData*);


    void setPlacement(const Transform& X_BG) {placement = X_BG;}
    const Transform& getPlacement() const    {return placement;}

    // This sets resolution to some factor times the object-specific default.
    // Anything 0 or less becomes -1 and means "use default".
    void setResolution(Real r) {
        resolution = r > 0 ? r : -1.;
    }
    Real getResolution() const {return resolution;}

    // This sets the scale to some factor times the default size of the object,
    // which will be somewhere around 1 length unit. Set to 0 or less to mean
    // "use default".
    void setScale(Real s) {
        scale = s > 0 ? s : -1.;
    }
    Real getScale() const {return scale;}


    void setColor(const Vec3& rgb) {
        assert(0<=rgb[0]&&rgb[0]<=1); // TODO
        assert(0<=rgb[1]&&rgb[1]<=1);
        assert(0<=rgb[2]&&rgb[2]<=1);
        colorRGB=rgb;
    }
    const Vec3& getColor() const {return colorRGB;}

    // Opacity should be greater than zero (invisible) and less than 
    // or equal to 1. The default will generally be 1, which is opaque,
    // but we use -1 to mean "use default" and let the client decide.
    void setOpacity(Real o) {
        opacity = o > 0 ? o : -1.;
    }
    Real getOpacity() const {return opacity;}

    void setRepresentationToPoints()     {representation=VTK_POINTS;}
    void setRepresentationToWireframe()  {representation=VTK_WIREFRAME;}
    void setRepresentationToSurface()    {representation=VTK_SURFACE;}
    void setRepresentationToUseDefault() {representation=-1;}

    DecorativeGeometryRep* clone() const {
        DecorativeGeometryRep* dup = cloneDecorativeGeometryRep();
        dup->clearMyHandle();
        dup->vtkObjects.resize(0);  // remove pointers to source geometry
        dup->generateVTKPipeline(); // make private VTK objects
        return dup;
    }
    virtual DecorativeGeometryRep* cloneDecorativeGeometryRep() const = 0;

    void setMyHandle(DecorativeGeometry& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

protected:
    void rememberVTKObject(vtkObject* o) {
        vtkObjects.push_back(o);
    }

private:
    friend class DecorativeGeometry;
    DecorativeGeometry* myHandle;     // the owner of this rep

    void deleteVTKGeometry() {
        // Delete in reverse order of allocation
        for (int i=(int)vtkObjects.size()-1; i >= 0; --i)
            vtkObjects[i]->Delete(), vtkObjects[i]=0;
        vtkObjects.resize(0);
    }

    // These will be handled as we generate the PolyData.
    Real      resolution;   // -1 means use default
    Real      scale;        // -1 means use default
    Transform placement;    // default is identity

    // These must wait until we are associated with an actor.
    Vec3 colorRGB;          // set R to -1 for "use default"
    Real opacity;           // -1 means "use default"
    int  representation;    // -1, VTK_POINTS, VTK_WIREFRAME, VTK_SURFACE

    // As we build the pipeline, we accumulate VTK objects which must
    // have their Delete() methods called in the desctructor.
    std::vector<vtkObject*> vtkObjects;
};


class DecorativeLineRep : public DecorativeGeometryRep {
public:
    DecorativeLineRep() : length(1) { }
    DecorativeLineRep(const Real& l) : length(l) {
        assert(l > 0); // TODO
    }

    // virtuals
    vtkPolyData* createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeLineRep(*this);
    }

    SimTK_DOWNCAST(DecorativeLineRep, DecorativeGeometryRep);
private:
    Real length;
};

class DecorativeCircleRep : public DecorativeGeometryRep {
public:
    DecorativeCircleRep() : r(1) { }
    DecorativeCircleRep(const Real& rad) : r(rad) {
        assert(r > 0); // TODO
    }

    const Real& getRadius() const {return r;}

    // virtuals
    vtkPolyData* createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeCircleRep(*this);
    }

    SimTK_DOWNCAST(DecorativeCircleRep, DecorativeGeometryRep);
private:
    Real r;
};

class DecorativeSphereRep : public DecorativeGeometryRep {
    static const int DefaultResolution = 15;
public:
    DecorativeSphereRep() : r(1) { }
    DecorativeSphereRep(const Real& rad) : r(rad) {
        assert(r > 0); // TODO
    }

    const Real& getRadius() const {return r;}

    // virtuals
    vtkPolyData* createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeSphereRep(*this);
    }

    SimTK_DOWNCAST(DecorativeSphereRep, DecorativeGeometryRep);
private:
    Real r;
};


class DecorativeBrickRep : public DecorativeGeometryRep {
public:
    DecorativeBrickRep() : halfLengths(0.5) { }
    DecorativeBrickRep(const Vec3& xyzLengths) : halfLengths(xyzLengths) {
        assert(halfLengths[0]>0&&halfLengths[1]>0&&halfLengths[2]>0); // TODO
    }

    const Vec3& getXYZHalfLengths() const {return halfLengths;}

    // virtuals
    vtkPolyData* createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeBrickRep(*this);
    }

    SimTK_DOWNCAST(DecorativeBrickRep, DecorativeGeometryRep);
private:
    Vec3 halfLengths;
};


class DecorativeCylinderRep : public DecorativeGeometryRep {
    static const int DefaultResolution = 10;
public:
    DecorativeCylinderRep() : radius(0.5), halfLength(0.5) { }
    DecorativeCylinderRep(Real r, Real h) 
      : radius(r), halfLength(h) {
        assert(radius>0&&halfLength>0); // TODO
    }

    Real getRadius()     const {return radius;}
    Real getHalfLength() const {return halfLength;}

    // virtuals
    vtkPolyData* createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeCylinderRep(*this);
    }

    SimTK_DOWNCAST(DecorativeCylinderRep, DecorativeGeometryRep);
private:
    Real radius, halfLength;
};

class DecorativeFrameRep : public DecorativeGeometryRep {
public:
    DecorativeFrameRep() : axisLength(1) { }
    DecorativeFrameRep(const Real& axisLen) : axisLength(axisLen) {
        assert(axisLen > 0); // TODO
    }

    const Real& getAxisLength() const {return axisLength;}

    // virtuals
    vtkPolyData* createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeFrameRep(*this);
    }

    SimTK_DOWNCAST(DecorativeFrameRep, DecorativeGeometryRep);
private:
    Real axisLength;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_
