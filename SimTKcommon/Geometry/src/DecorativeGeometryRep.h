#ifndef SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_REP_H_
#define SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_REP_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Jack Middleton, Peter Eastman, Ayman Habib                   *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/DecorativeGeometry.h"

#include <cmath>

namespace SimTK {

class DecorativeGeometryRep {
public:
    DecorativeGeometryRep() 
    :   body(0), indexOnBody(-1), userRef(0), 
        placement(), scaleFactors(-1,-1,-1), resolution(-1),
        colorRGB(-1,-1,-1), opacity(-1), lineThickness(-1), faceCamera(-1),
        representation(DecorativeGeometry::DrawDefault), myHandle(0)
    {}

    virtual ~DecorativeGeometryRep() {
        clearMyHandle();
    }

    virtual DecorativeGeometryRep* cloneDecorativeGeometryRep() const = 0;
    virtual void implementGeometry(DecorativeGeometryImplementation&) const = 0;

    void setBodyId(int b) {body = b;}
    int getBodyId() const {return body;}
    void setIndexOnBody(int x) {indexOnBody = x;}
    int getIndexOnBody() const {return indexOnBody;}
    void setUserRef(void* p) {userRef = p;}
    void* getUserRef() const {return userRef;}

    void setTransform(const Transform& X_BD) {placement = X_BD;}
    const Transform& getTransform() const    {return placement;}

    // This sets resolution to some factor times the object-specific default.
    // Anything 0 or less becomes -1 and means "use default".
    void setResolution(Real r) {
        resolution = r > 0 ? r : Real(-1);
    }

    Real getResolution() const {return resolution;}

    // This sets the scale to some factor times the default size of the object,
    // which will be somewhere around 1 length unit. Set to 0 or less to mean
    // "use default".
    void setScaleFactors(const Vec3& s) {
        for (int i=0; i<3; ++i)
            scaleFactors[i] = s[i] > 0 ? s[i] : Real(-1);
    }
    const Vec3& getScaleFactors() const {return scaleFactors;}


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
        opacity = o > 0 ? o : Real(-1);
    }
    Real getOpacity() const {return opacity;}

    void setLineThickness(Real t) {
        lineThickness = t > 0 ? t : Real(-1);
    }
    Real getLineThickness() const {return lineThickness;}
    
    void setFaceCamera(int shouldFace) {
        faceCamera = shouldFace==0 ? 0 : (shouldFace>0 ? 1 : -1);
    }
    
    int getFaceCamera() const {return faceCamera;}

    void setRepresentation(const DecorativeGeometry::Representation& r) 
    {   representation=r; }
    DecorativeGeometry::Representation getRepresentation() const 
    {   return representation; }

    DecorativeGeometryRep* clone() const {
        DecorativeGeometryRep* dup = cloneDecorativeGeometryRep();
        dup->clearMyHandle();
        return dup;
    }

    void setMyHandle(DecorativeGeometry& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

    // Used by the composite Decorations object:
    // Set any properties that still have their default values to the ones
    // from src. Source body, index, and userRef are transferred 
    // unconditionally; source placement is composed with this one 
    // unconditionally. Scale factors, resolution, opacity, and line thickness 
    // are composed, with a default (-1) value treated as 1.
    void inheritPropertiesFrom(const DecorativeGeometryRep& srep) {
        body        = srep.body;
        indexOnBody = srep.indexOnBody;
        userRef     = srep.userRef;

        placement = srep.placement * placement;

        // These are multiplied together if both are valid, otherwise we
        // take the valid one, or set to -1 if neither is valid.
        for (int i=0; i<3; ++i)
            scaleFactors[i] = compose(scaleFactors[i], srep.scaleFactors[i]);
        resolution    = compose(resolution,    srep.resolution);
        opacity       = compose(opacity,       srep.opacity);
        lineThickness = compose(lineThickness, srep.lineThickness);

        // These are left alone if already specified in the destination,
        // otherwise they are given the source value.
        if (colorRGB[0] == -1) colorRGB = srep.colorRGB;
        if (faceCamera == -1) faceCamera = srep.faceCamera;
        if (representation == DecorativeGeometry::DrawDefault)
            representation = srep.representation;
    }
private:
    // Given arguments where ai < 0 means it hasn't been specified,
    // compose them if they are both valid, otherwise return the valid
    // one if there is one, otherwise return -1. Note that we're treating
    // 0 as specified; composing it with anything will return 0.
    static Real compose(Real a1, Real a2) {
        if (a1 >= 0) return a2 >= 0 ? a1*a2 : a1;
        if (a2 >= 0) return a2;
        return -1;
    }

protected:
    friend class DecorativeGeometry;

    int         body;
    int         indexOnBody;
    void*       userRef;

    Transform   placement;          // default is identity
    Vec3        scaleFactors;       // -1 means use default in that direction

    Real        resolution;         // -1 means use default
    Vec3        colorRGB;           // set R to -1 for "use default"
    Real        opacity;            // -1 means "use default"
    Real        lineThickness;      // -1 means "use default"
    int         faceCamera;
    DecorativeGeometry::Representation  
                representation;     // e.g. points, wireframe, surface

    DecorativeGeometry* myHandle;   // the owner of this rep
};

    ////////////////////////
    // DecorativePointRep //
    ////////////////////////

class DecorativePointRep : public DecorativeGeometryRep {
public:
    // no default constructor
    DecorativePointRep(const Vec3& p) : point(p) {}

    void setPoint(const Vec3& p) {point=p; }

    const Vec3& getPoint() const {return point;}

    // virtuals
    DecorativePointRep* cloneDecorativeGeometryRep() const {
        DecorativePointRep* DGRep = new DecorativePointRep(*this);
        return DGRep; 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementPointGeometry(getMyPointHandle());
    }

    SimTK_DOWNCAST(DecorativePointRep, DecorativeGeometryRep);
private:
    Vec3 point;

    // This is just a static downcast since the DecorativeGeometry handle class
    // is not virtual.
    const DecorativePoint& getMyPointHandle() const {
        return *static_cast<const DecorativePoint*>(myHandle);
    }
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
        return *static_cast<const DecorativeLine*>(myHandle);
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
        return *static_cast<const DecorativeCircle*>(myHandle);
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
        return *static_cast<const DecorativeSphere*>(myHandle);
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
        return *static_cast<const DecorativeEllipsoid*>(myHandle);
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
        return *static_cast<const DecorativeBrick*>(myHandle);
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
        return *static_cast<const DecorativeCylinder*>(myHandle);
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
        return *static_cast<const DecorativeFrame*>(myHandle);
    }
};

///////////////////////
// DecorativeTextRep //
///////////////////////

class DecorativeTextRep : public DecorativeGeometryRep {
static const int DefaultResolution = 15;
public:
// no default constructor
explicit DecorativeTextRep(const std::string& label) 
: text(label), isScreenText(false) {}

void setText(const std::string& label) {
    text = label;
}
const std::string& getText() const {
    return text;
}

void setIsScreenText(bool isScreen) {isScreenText=isScreen;}
bool getIsScreenText() const {return isScreenText;}

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
bool        isScreenText; // in screen coordinates

// This is just a static downcast since the DecorativeGeometry handle class is not virtual.
const DecorativeText& getMyTextHandle() const {
    return *static_cast<const DecorativeText*>(myHandle);
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
    return *static_cast<const DecorativeMesh*>(myHandle);
}
};

///////////////////////////
// DecorativeMeshFileRep //
///////////////////////

class DecorativeMeshFileRep : public DecorativeGeometryRep {
public:
// no default constructor
explicit DecorativeMeshFileRep(const std::string& meshFileName) : meshFile(meshFileName) {
}

const std::string& getMeshFile() const {
    return  meshFile;
}

// virtuals
DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
    DecorativeMeshFileRep* DGRep = new DecorativeMeshFileRep(*this);
    return DGRep ; 
}

void implementGeometry(DecorativeGeometryImplementation& geometry) const {
    geometry.implementMeshFileGeometry(getMyMeshFileHandle());
}

SimTK_DOWNCAST(DecorativeMeshFileRep, DecorativeGeometryRep);
private:
std::string meshFile;

// This is just a static downcast since the DecorativeGeometry handle class is not virtual.

const DecorativeMeshFile& getMyMeshFileHandle() const {
    return *static_cast<const DecorativeMeshFile*>(myHandle);
}
};


////////////////////////
// DecorativeTorusRep //
////////////////////////

class DecorativeTorusRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeTorusRep(Real inner, Real outer) : innerRadius(inner), outerRadius(outer) {
    }

    const Real& getInnerRadius() const {
        return  innerRadius;
    }

    const Real& getOuterRadius() const {
        return  outerRadius;
    }

    void setInnerRadius(Real ir) {
        innerRadius = ir;
    }

    void setOuterRadius(Real orad) {
        outerRadius = orad;
    }
    // virtuals
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        DecorativeTorusRep* DGRep = new DecorativeTorusRep(*this);
        return DGRep;
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementTorusGeometry(getMyTorusHandle());
    }

    SimTK_DOWNCAST(DecorativeTorusRep, DecorativeGeometryRep);
private:
    Real innerRadius;
    Real outerRadius;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.

    const DecorativeTorus& getMyTorusHandle() const {
        return *static_cast<const DecorativeTorus*>(myHandle);
    }
};


///////////////////////
// DecorativArrowRep //
///////////////////////

class DecorativeArrowRep : public DecorativeGeometryRep {
public:
    // no default constructor
    explicit DecorativeArrowRep(const Vec3& orig, const Vec3& dir, Real& len) : origin(orig), direction(dir), length(len) {
    }

    const Vec3& getOrigin() const {
        return  origin;
    }

    const Vec3& getDirection() const {
        return  direction;
    }

    const Real& getLength() const {
        return  length;
    }

    void setOrigin(const Vec3& orig) {
        origin = orig;
    }
    void setDirection(const Vec3& dir) {
        direction = dir;
    }
    void setLength(Real& len) {
        length = len;
    }

    // virtuals
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        DecorativeArrowRep* DGRep = new DecorativeArrowRep(*this);
        return DGRep;
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        geometry.implementArrowGeometry(getMyArrowHandle());
    }

    SimTK_DOWNCAST(DecorativeArrowRep, DecorativeGeometryRep);
private:
    Vec3 origin;
    Vec3 direction;
    Real length;

    // This is just a static downcast since the DecorativeGeometry handle class is not virtual.

    const DecorativeArrow& getMyArrowHandle() const {
        return *static_cast<const DecorativeArrow*>(myHandle);
    }
};

    ////////////////////
    // DecorationsRep //
    ////////////////////

class DecorationsRep : public DecorativeGeometryRep {
public:
    DecorationsRep() {}

    int addDecoration(const DecorativeGeometry& decoration) {
        geom.push_back(decoration);
        return (int)(geom.size()-1);
    }
    int addDecoration(const Transform& X_DE,
                      const DecorativeGeometry& element) {
        geom.push_back(element);
        // The current transform goes from the element's frame E to the 
        // frame G of the actual geometry. We would normally put E at the
        // Decorations frame D, but now we want to relocate it.
        const Transform& X_EG = geom.back().getTransform();
        geom.back().setTransform(X_DE*X_EG);
        return (int)(geom.size()-1);
    }

    int getNumDecorations() const {return (int)geom.size();}
    const DecorativeGeometry& getDecoration(int i) const {return geom[i];}

    // virtuals

    DecorationsRep* cloneDecorativeGeometryRep() const {
        DecorationsRep* DGRep = new DecorationsRep(*this);
        return DGRep; 
    }

    void implementGeometry(DecorativeGeometryImplementation& geometry) const {
        for (unsigned i=0; i < geom.size(); ++i) {
            DecorativeGeometry copy = geom[i];
            copy.updRep().inheritPropertiesFrom(*this);
            copy.implementGeometry(geometry);
        }
    }

    SimTK_DOWNCAST(DecorationsRep, DecorativeGeometryRep);
private:
    Array_<DecorativeGeometry> geom;

    // This is just a static downcast since the DecorativeGeometry handle class
    // is not virtual.
    const Decorations& getMyHandle() const {
        return *static_cast<const Decorations*>(myHandle);
    }
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_REP_H_
