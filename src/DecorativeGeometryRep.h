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

class vtkPolyData;

namespace SimTK {

static const Real Pi = std::acos(Real(-1));

class DecorativeGeometryRep {
public:
    DecorativeGeometryRep() 
      : myHandle(0), colorRGB(0), opacity(1), resolution(0), scale(1) 
    { 
    }
    virtual ~DecorativeGeometryRep() {clearMyHandle();}

    virtual vtkPolyData* createVTKPolyData() const = 0;

    void setPlacement(const Transform& X_BG) {placement = X_BG;}
    const Transform& getPlacement() const    {return placement;}


    void setColor(const Vec3& rgb) {
        assert(0<=rgb[0]&&rgb[0]<=1); // TODO
        assert(0<=rgb[1]&&rgb[1]<=1);
        assert(0<=rgb[2]&&rgb[2]<=1);
        colorRGB=rgb;
    }
    const Vec3& getColor() const {return colorRGB;}

    void setOpacity(Real o) {
        assert(0<=o && o <= 1); // TODO
        opacity=o;
    }
    Real getOpacity() const {return opacity;}

    void setResolution(int r) {
        assert(0<=r);
        resolution=r;   // 0 means use default
    }
    int getResolution() const {return resolution;}

    void setScale(Real s) {
        assert(0<s); // TODO
        scale=s;
    }
    Real getScale() const {return scale;}


    DecorativeGeometryRep* clone() const {
        DecorativeGeometryRep* dup = cloneDecorativeGeometryRep();
        dup->clearMyHandle();
        return dup;
    }
    virtual DecorativeGeometryRep* cloneDecorativeGeometryRep() const = 0;

    void setMyHandle(DecorativeGeometry& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class DecorativeGeometry;
    DecorativeGeometry* myHandle;     // the owner of this rep

    Vec3 colorRGB;   // default is 000 (black)
    Real opacity;    // default is 1
    int  resolution; // 0 means use default
    Real scale;      // default is 1

    Transform placement;    // default is identity
};


class DecorativeLineRep : public DecorativeGeometryRep {
public:
    DecorativeLineRep() : length(1) { }
    DecorativeLineRep(const Real& l) : length(l) {
        assert(l > 0); // TODO
    }

    // virtuals
    vtkPolyData* createVTKPolyData() const;
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
    vtkPolyData* createVTKPolyData() const;
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeCircleRep(*this);
    }

    SimTK_DOWNCAST(DecorativeCircleRep, DecorativeGeometryRep);
private:
    Real r;
};

class DecorativeSphereRep : public DecorativeGeometryRep {
public:
    DecorativeSphereRep() : r(1) { }
    DecorativeSphereRep(const Real& rad) : r(rad) {
        assert(r > 0); // TODO
    }

    const Real& getRadius() const {return r;}

    // virtuals
    vtkPolyData* createVTKPolyData() const;
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeSphereRep(*this);
    }

    SimTK_DOWNCAST(DecorativeSphereRep, DecorativeGeometryRep);
private:
    Real r;
};


class DecorativeBrickRep : public DecorativeGeometryRep {
public:
    DecorativeBrickRep() : lengths(1) { }
    DecorativeBrickRep(const Vec3& xyzLengths) : lengths(xyzLengths) {
        assert(lengths[0]>0&&lengths[1]>0&&lengths[2]>0); // TODO
    }

    const Vec3& getXYZLengths() const {return lengths;}

    // virtuals
    vtkPolyData* createVTKPolyData() const;
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeBrickRep(*this);
    }

    SimTK_DOWNCAST(DecorativeBrickRep, DecorativeGeometryRep);
private:
    Vec3 lengths;
};

class DecorativeFrameRep : public DecorativeGeometryRep {
public:
    DecorativeFrameRep() : axisLength(1) { }
    DecorativeFrameRep(const Real& axisLen) : axisLength(axisLen) {
        assert(axisLen > 0); // TODO
    }

    const Real& getAxisLength() const {return axisLength;}

    // virtuals
    vtkPolyData* createVTKPolyData() const;
    DecorativeGeometryRep* cloneDecorativeGeometryRep() const {
        return new DecorativeFrameRep(*this);
    }

    SimTK_DOWNCAST(DecorativeFrameRep, DecorativeGeometryRep);
private:
    Real axisLength;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DECORATIVE_GEOMETRY_REP_H_
