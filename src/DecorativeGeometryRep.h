#ifndef SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_
#define SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/AnalyticGeometry.h"

#include <cmath>
#include <vector>

#include "vtkPolyData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkTransform.h"
#include "vtkProperty.h"
#include "vtkObject.h"

namespace SimTK {

static const Real Pi = std::acos(Real(-1));

class DecorativeGeometryRep {
public:
    DecorativeGeometryRep() 
      : myHandle(0),  resolution(-1), scale(-1), placement(), 
        colorRGB(-1,-1,-1), opacity(-1), lineThickness(-1), representation(-1)
    { 
    }
    virtual ~DecorativeGeometryRep() {
        deleteVTKGeometry();
        clearMyHandle();
    }

    void generateVTKPipeline() {
        deleteVTKGeometry();
        createVTKPolyData();
    }

    void updateVTKPipeline() {
        // TODO
        generateVTKPipeline(); // start from scratch (bad)
    }

    vtkPolyData* updVTKPolyData() {
        assert(vtkObjects.size());
        return vtkPolyDataAlgorithm::SafeDownCast(vtkObjects.back())->GetOutput();
    }


    // Caller must be sure to call VTK's Delete() methods on these objects
    // when done with them.
    virtual void createVTKPolyData() = 0;

    // Combine a SimTK Transform with a scale factor and return the equivalent
    // vtkTransform that will translate, rotate, and scale the same way.
    vtkTransform* createVTKTransform(const Transform&, const Real&);

    // Take some polygon data as input, copy it and transform it, returning the
    // transformed polygon data.
    vtkPolyData* transformVTKPolyData(const Transform&, const Real&, vtkPolyData*);


    void setPlacement(const Transform& X_BG) {
        placement = X_BG;
        updateVTKPipeline();
    }
    const Transform& getPlacement() const    {return placement;}

    // This sets resolution to some factor times the object-specific default.
    // Anything 0 or less becomes -1 and means "use default".
    void setResolution(Real r) {
        resolution = r > 0 ? r : -1.;
        updateVTKPipeline();
    }
    Real getResolution() const {return resolution;}

    // This sets the scale to some factor times the default size of the object,
    // which will be somewhere around 1 length unit. Set to 0 or less to mean
    // "use default".
    void setScale(Real s) {
        scale = s > 0 ? s : -1.;
        updateVTKPipeline();
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

    void setLineThickness(Real t) {
        lineThickness = t > 0 ? t : -1.;
    }
    Real getLineThickness() const {return lineThickness;}

    void setRepresentationToPoints()     {representation=VTK_POINTS;}
    void setRepresentationToWireframe()  {representation=VTK_WIREFRAME;}
    void setRepresentationToSurface()    {representation=VTK_SURFACE;}
    void setRepresentationToUseDefault() {representation=-1;}

    int getRepresentation() const {return representation;}

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
        for (int i=(int)vtkObjects.size()-1; i >= 0; --i) {
            vtkObject* obj = vtkObjects[i];
            //std::cout << "ABOUT TO DELETE\n";
            //obj->Print(std::cout);
            obj->Delete();
            vtkObjects[i]=0;
        }
        vtkObjects.resize(0);
    }

    // These will be handled as we generate the PolyData.
    Real      resolution;   // -1 means use default
    Real      scale;        // -1 means use default
    Transform placement;    // default is identity

    // These must wait until we are associated with an actor.
    Vec3 colorRGB;          // set R to -1 for "use default"
    Real opacity;           // -1 means "use default"
    Real lineThickness;     // -1 means "use default"
    int  representation;    // -1, VTK_POINTS, VTK_WIREFRAME, VTK_SURFACE

    // As we build the pipeline, we accumulate VTK objects which must
    // have their Delete() methods called in the desctructor.
    std::vector<vtkObject*> vtkObjects;
};

    ///////////////////////
    // DecorativeLineRep //
    ///////////////////////

class DecorativeLineRep : public DecorativeGeometryRep {
public:
    // no default constructor
    DecorativeLineRep(const Vec3& p1, const Vec3& p2) : point1(p1), point2(p2) {
        generateVTKPipeline();
    }

    void setPoint1(const Vec3& p) {point1=p; updateVTKPipeline();}
    void setPoint2(const Vec3& p) {point2=p; updateVTKPipeline();}
    void setEndpoints(const Vec3& p1, const Vec3& p2) {
        point1=p1; point2=p2; updateVTKPipeline();
    }

    const Vec3& getPoint1() const {return point1;}
    const Vec3& getPoint2() const {return point2;}

    // virtuals
    void createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeLineRep(*this);
    }

    SimTK_DOWNCAST(DecorativeLineRep, DecorativeGeometryRep);
private:
    Vec3 point1, point2;
};

    /////////////////////////
    // DecorativeCircleRep //
    /////////////////////////

class DecorativeCircleRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeCircleRep(const Real& rad) : r(rad) {
        assert(r > 0); // TODO
        generateVTKPipeline();
    }

    void setRadius(const Real& rad) {
        assert(rad > 0); // TODO;
        r = rad;
        updateVTKPipeline();
    }
    const Real& getRadius() const {return r;}

    // virtuals
    void createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeCircleRep(*this);
    }

    SimTK_DOWNCAST(DecorativeCircleRep, DecorativeGeometryRep);
private:
    Real r;
};

    /////////////////////////
    // DecorativeSphereRep //
    /////////////////////////

class DecorativeSphereRep : public DecorativeGeometryRep {
    static const int DefaultResolution = 15;
public:
    // no default constructor
    explicit DecorativeSphereRep(const Real& rad) : r(rad) {
        assert(r > 0); // TODO
        generateVTKPipeline();
    }

    void setRadius(const Real& rad) {
        assert(rad > 0); // TODO;
        r = rad;
        updateVTKPipeline();
    }
    const Real& getRadius() const {return r;}

    // virtuals
    void createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeSphereRep(*this);
    }

    SimTK_DOWNCAST(DecorativeSphereRep, DecorativeGeometryRep);
private:
    Real r;
};

    ////////////////////////
    // DecorativeBrickRep //
    ////////////////////////

class DecorativeBrickRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeBrickRep(const Vec3& xyzHalfLengths) : halfLengths(xyzHalfLengths) {
        assert(halfLengths[0]>0&&halfLengths[1]>0&&halfLengths[2]>0); // TODO
        generateVTKPipeline();
    }

    void setHalfLengths(const Vec3& hl) {
        assert(hl[0]>0&&hl[1]>0&&hl[2]>0); // TODO;
        halfLengths = hl;
        updateVTKPipeline();
    }
    const Vec3& getHalfLengths() const {return halfLengths;}

    // virtuals
    void createVTKPolyData();
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
    // no default constructor
    DecorativeCylinderRep(Real r, Real hh) 
      : radius(r), halfHeight(hh) {
        assert(radius>0&&halfHeight>0); // TODO
        generateVTKPipeline();
    }

    void setRadius(const Real& rad) {
        assert(rad > 0); // TODO;
        radius = rad;
        updateVTKPipeline();
    }
    void setHalfHeight(const Real& hh) {
        assert(hh > 0); // TODO;
        halfHeight = hh;
        updateVTKPipeline();
    }
    Real getRadius()     const {return radius;}
    Real getHalfHeight() const {return halfHeight;}

    // virtuals
    void createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeCylinderRep(*this);
    }

    SimTK_DOWNCAST(DecorativeCylinderRep, DecorativeGeometryRep);
private:
    Real radius, halfHeight;
};

class DecorativeFrameRep : public DecorativeGeometryRep {
public:
    // no default constructor
    DecorativeFrameRep(const Real& len) : axisLength(len) {
        assert(len > 0); // TODO
        generateVTKPipeline();
    }

    void setAxisLength(const Real& len) {
        assert(len > 0); // TODO;
        axisLength = len;
        updateVTKPipeline();
    }
    const Real& getAxisLength() const {return axisLength;}

    // virtuals
    void createVTKPolyData();
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeFrameRep(*this);
    }

    SimTK_DOWNCAST(DecorativeFrameRep, DecorativeGeometryRep);
private:
    Real axisLength;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_
