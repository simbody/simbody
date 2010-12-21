#ifndef SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_REP_H_
#define SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Jack Middleton                                               *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/DecorativeGeometry.h"

#include <cmath>

namespace SimTK {

class DecorativeGeometryRep {
public:
    DecorativeGeometryRep() 
      : myHandle(0), body(0), placement(), resolution(-1), scale(-1),
      colorRGB(-1,-1,-1), opacity(-1), lineThickness(-1), faceCamera(false), representation(DecorativeGeometry::DrawDefault)
    { 
    }

    ~DecorativeGeometryRep() {
        clearMyHandle();
    }

    void setBodyId(int b) {
        body = b;
    }
    int getBodyId() const {return body;}

    void setTransform(const Transform& X_BD) {
        placement = X_BD;
    }
    const Transform& getTransform() const    {return placement;}

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

    void setLineThickness(Real t) {
        lineThickness = t > 0 ? t : -1.;
    }
    Real getLineThickness() const {return lineThickness;}
    
    void setFaceCamera(bool face) {
        faceCamera = face;
    }
    
    bool getFaceCamera() const {return faceCamera;}

    void setRepresentation(const DecorativeGeometry::Representation& r) {representation=r;}
    DecorativeGeometry::Representation getRepresentation() const {return representation;}

    DecorativeGeometryRep* clone() const {
        DecorativeGeometryRep* dup = cloneDecorativeGeometryRep();
        dup->clearMyHandle();
        return dup;
    }

    virtual DecorativeGeometryRep* cloneDecorativeGeometryRep() const = 0;

    virtual void implementGeometry(DecorativeGeometryImplementation&) const = 0;

    void setMyHandle(DecorativeGeometry& h) {
        myHandle = &h;
    }
    void clearMyHandle() {
        myHandle=0;
    }

private:
    friend class DecorativeGeometry;

    int       body;
    Transform placement;    // default is identity
    Real      resolution;   // -1 means use default
    Real      scale;        // -1 means use default

    Vec3 colorRGB;          // set R to -1 for "use default"
    Real opacity;           // -1 means "use default"
    Real lineThickness;     // -1 means "use default"
    bool faceCamera;
    DecorativeGeometry::Representation  representation; // e.g. points, wireframe, surface

protected:
    DecorativeGeometry* myHandle;         // the owner of this rep
};

    ///////////////////////
    // DecorativeLineRep //
    ///////////////////////

class DecorativeLineRep : public DecorativeGeometryRep {
public:
    // no default constructor
    DecorativeLineRep( const Vec3& p1, const Vec3& p2) : point1(p1), point2(p2) {
    }

    void setPoint1(const Vec3& p) {point1=p; }
    void setPoint2(const Vec3& p) {point2=p; }
    void setEndpoints(const Vec3& p1, const Vec3& p2) { point1=p1; point2=p2; }

    const Vec3& getPoint1() const {return point1;}
    const Vec3& getPoint2() const {return point2;}

    // virtuals
    DecorativeLineRep* cloneDecorativeGeometryRep() const {
        DecorativeLineRep* DGRep = new DecorativeLineRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementLineGeometry(getMyLineHandle());
    }

    SimTK_DOWNCAST(DecorativeLineRep, DecorativeGeometryRep);
private:
    Vec3 point1, point2;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeLine& getMyLineHandle() const {
        return *reinterpret_cast<const DecorativeLine*>(myHandle);
    }
};

    /////////////////////////
    // DecorativeCircleRep //
    /////////////////////////

class DecorativeCircleRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeCircleRep( const Real& rad) : r(rad) {
        assert(r > 0); // TODO
    }

    void setRadius(const Real& rad) {
        assert(rad > 0); // TODO;
        r = rad;
    }
    const Real& getRadius() const {return r;}

    // virtuals
    DecorativeCircleRep* cloneDecorativeGeometryRep() const {
        DecorativeCircleRep* DGRep = new DecorativeCircleRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementCircleGeometry(getMyCircleHandle());
    }

    SimTK_DOWNCAST(DecorativeCircleRep, DecorativeGeometryRep);
private:
    Real r;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeCircle& getMyCircleHandle() const {
        return *reinterpret_cast<const DecorativeCircle*>(myHandle);
    }
};

    /////////////////////////
    // DecorativeSphereRep //
    /////////////////////////

class DecorativeSphereRep : public DecorativeGeometryRep {
    static const int DefaultResolution = 15;
public:
    // no default constructor
    explicit DecorativeSphereRep( const Real& rad) : r(rad) {
        assert(r > 0); // TODO
    }

    void setRadius(const Real& rad) {
        assert(rad > 0); // TODO;
        r = rad;
    }
    const Real& getRadius() const {return r;}

    // virtuals
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        DecorativeSphereRep* DGRep = new DecorativeSphereRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementSphereGeometry(getMySphereHandle());
    }

    SimTK_DOWNCAST(DecorativeSphereRep, DecorativeGeometryRep);
private:
    Real r;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeSphere& getMySphereHandle() const {
        return *reinterpret_cast<const DecorativeSphere*>(myHandle);
    }
};


    ////////////////////////////
    // DecorativeEllipsoidRep //
    ////////////////////////////

class DecorativeEllipsoidRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeEllipsoidRep( const Vec3& xyzRadii) : radii(xyzRadii) {
        assert(radii[0]>0&&radii[1]>0&&radii[2]>0); // TODO
    }

    void setRadii(const Vec3& r) {
        assert(r[0]>0&&r[1]>0&&r[2]>0); // TODO;
        radii = r;
    }
    const Vec3& getRadii() const {return radii;}

    // virtuals
    DecorativeEllipsoidRep* cloneDecorativeGeometryRep() const {
        DecorativeEllipsoidRep* DGRep = new DecorativeEllipsoidRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementEllipsoidGeometry(getMyEllipsoidHandle());
    }

    SimTK_DOWNCAST(DecorativeEllipsoidRep, DecorativeGeometryRep);
private:
    Vec3 radii;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeEllipsoid& getMyEllipsoidHandle() const {
        return *reinterpret_cast<const DecorativeEllipsoid*>(myHandle);
    }
};


    ////////////////////////
    // DecorativeBrickRep //
    ////////////////////////

class DecorativeBrickRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeBrickRep( const Vec3& xyzHalfLengths) : halfLengths(xyzHalfLengths) {
        assert(halfLengths[0]>0&&halfLengths[1]>0&&halfLengths[2]>0); // TODO
    }

    void setHalfLengths(const Vec3& hl) {
        assert(hl[0]>0&&hl[1]>0&&hl[2]>0); // TODO;
        halfLengths = hl;
    }
    const Vec3& getHalfLengths() const {return halfLengths;}

    // virtuals
    DecorativeBrickRep* cloneDecorativeGeometryRep() const {
        DecorativeBrickRep* DGRep = new DecorativeBrickRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementBrickGeometry(getMyBrickHandle());
    }

    SimTK_DOWNCAST(DecorativeBrickRep, DecorativeGeometryRep);
private:
    Vec3 halfLengths;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeBrick& getMyBrickHandle() const {
        return *reinterpret_cast<const DecorativeBrick*>(myHandle);
    }
};


class DecorativeCylinderRep : public DecorativeGeometryRep {
    static const int DefaultResolution = 10;
public:
    // no default constructor
    DecorativeCylinderRep( Real r, Real hh) 
      : radius(r), halfHeight(hh) {
        assert(radius>0&&halfHeight>0); // TODO
    }

    void setRadius(const Real& rad) {
        assert(rad > 0); // TODO;
        radius = rad;
    }
    void setHalfHeight(const Real& hh) {
        assert(hh > 0); // TODO;
        halfHeight = hh;
    }
    Real getRadius()     const {return radius;}
    Real getHalfHeight() const {return halfHeight;}

    // virtuals
    DecorativeCylinderRep* cloneDecorativeGeometryRep() const {
        DecorativeCylinderRep* DGRep = new DecorativeCylinderRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementCylinderGeometry(getMyCylinderHandle());
    }

    SimTK_DOWNCAST(DecorativeCylinderRep, DecorativeGeometryRep);
private:
    Real radius, halfHeight;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeCylinder& getMyCylinderHandle() const {
        return *reinterpret_cast<const DecorativeCylinder*>(myHandle);
    }
};

class DecorativeFrameRep : public DecorativeGeometryRep {
public:
    // no default constructor
    DecorativeFrameRep(const Real& len) : axisLength(len) {
        assert(len > 0); // TODO
    }

    void setAxisLength(const Real& len) {
        assert(len > 0); // TODO;
        axisLength = len;
    }
    const Real& getAxisLength() const {return axisLength;}

    // virtuals
    DecorativeFrameRep* cloneDecorativeGeometryRep() const {
        DecorativeFrameRep* DGRep = new DecorativeFrameRep(*this);
        return( DGRep ); 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementFrameGeometry(getMyFrameHandle());
    }

    SimTK_DOWNCAST(DecorativeFrameRep, DecorativeGeometryRep);
private:
    Real axisLength;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.
    const DecorativeFrame& getMyFrameHandle() const {
        return *reinterpret_cast<const DecorativeFrame*>(myHandle);
    }
};

///////////////////////
// DecorativeTextRep //
///////////////////////

class DecorativeTextRep : public DecorativeGeometryRep {
static const int DefaultResolution = 15;
public:
// no default constructor
explicit DecorativeTextRep(const std::string& label) : text(label) {
}

void setText(const std::string& label) {
    text = label;
}
const std::string& getText() const {
    return text;
}

// virtuals
DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
    DecorativeTextRep* DGRep = new DecorativeTextRep(*this);
    return( DGRep ); 
}

void implementGeometry(DecorativeGeometryImplementation& geometry) const {
    geometry.implementTextGeometry(getMyTextHandle());
}

SimTK_DOWNCAST(DecorativeTextRep, DecorativeGeometryRep);
private:
std::string text;

// This is just a static downcast since the DecorativeGeometry handle class is not virtual.
const DecorativeText& getMyTextHandle() const {
    return *reinterpret_cast<const DecorativeText*>(myHandle);
}
};

///////////////////////
// DecorativeMeshRep //
///////////////////////

class DecorativeMeshRep : public DecorativeGeometryRep {
public:
// no default constructor
explicit DecorativeMeshRep(const PolygonalMesh& mesh) : mesh(mesh) {
}

const PolygonalMesh& getMesh() const {
    return  mesh;
}

// virtuals
DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
    DecorativeMeshRep* DGRep = new DecorativeMeshRep(*this);
    return( DGRep ); 
}

void implementGeometry(DecorativeGeometryImplementation& geometry) const {
    geometry.implementMeshGeometry(getMyMeshHandle());
}

SimTK_DOWNCAST(DecorativeMeshRep, DecorativeGeometryRep);
private:
PolygonalMesh mesh;

// This is just a static downcast since the DecorativeGeometry handle class is not virtual.
const DecorativeMesh& getMyMeshHandle() const {
    return *reinterpret_cast<const DecorativeMesh*>(myHandle);
}
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_REP_H_
