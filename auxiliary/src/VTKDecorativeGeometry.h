#ifndef SimTK_SIMBODY_VTK_DECORATIVE_GEOMETRY_H_
#define SimTK_SIMBODY_VTK_DECORATIVE_GEOMETRY_H_

/* Portions copyright (c) 2007 Stanford University and Jack Middleton.
 * Contributors: Michael Sherman
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

/** @file
 * This is an implementation of the SimTK::DecorativeGeometry facility using
 * Kitware's VTK visualization system. The resulting VTK objects can be
 * used by a VTK application or reporter to display Simbody DecorativeGeometry
 * visualization objects.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/DecorativeGeometry.h"

#include <cmath>
#include <vector>

class vtkPolyData;
class vtkTransform;
class vtkObject;

using namespace SimTK;

// Each object of class VTKDecorativeGeometry implements a single SimTK::DecorativeGeometry
// object. Use it like this: (1) We are given an object DG of class DecorativeGeometry, but
// we don't know what it is specifically. (2) Instantiate an object of type VTKDecorativeGeometry,
// which implements the abstract class DecorativeGeometryImplementation. (3) Call DG's 
// implementGeometry() method, passing it the VTKDecorativeGeometry object. (4) DG knows
// its own type, so it can call the appropriate "implement" method defined below, allowing
// the VTKDecorativeGeometry object to be constructed as the right kind of object. 
//
// Here we've chosen vtkPolyData as the abstract representation of geometry, so we have our
// various implementXXXGeometry() methods build a VTK pipeline appropriate to the specific
// DecorativeGeometry object.
// 
class VTKDecorativeGeometry : public SimTK::DecorativeGeometryImplementation {
public:
    VTKDecorativeGeometry() { }

    // Implement the virtual methods required for a DecorativeGeometryImplementation.
    /*virtual*/ ~VTKDecorativeGeometry() {
        deleteVTKGeometry();
    }
    /*virtual*/ void implementLineGeometry     (const DecorativeLine&);
    /*virtual*/ void implementBrickGeometry    (const DecorativeBrick&);
    /*virtual*/ void implementCylinderGeometry (const DecorativeCylinder&);
    /*virtual*/ void implementCircleGeometry   (const DecorativeCircle&); 
    /*virtual*/ void implementSphereGeometry   (const DecorativeSphere&);
    /*virtual*/ void implementEllipsoidGeometry(const DecorativeEllipsoid&);
    /*virtual*/ void implementFrameGeometry    (const DecorativeFrame&);

    // The last vtkObject is the end of the VTK pipeline -- its output is the
    // final representation of the object.
    vtkPolyData* getVTKPolyData();

    // Combine a SimTK Transform with x,y,z scale factors and return the equivalent
    // vtkTransform that will translate, rotate, and scale the same way.
    vtkTransform* createVTKTransform(const Transform&, const Vec3& scale);

    // Take some polygon data as input, copy it and transform it, returning the
    // transformed polygon data.
    vtkPolyData* transformVTKPolyData(const Transform&, const Vec3& scale, vtkPolyData*);

protected:
    void rememberVTKObject(vtkObject* o) {
        vtkObjects.push_back(o);
    }

    void deleteVTKGeometry();

    // As we build the pipeline, we accumulate VTK objects which must
    // have their Delete() methods called in the destructor.
    std::vector<vtkObject*> vtkObjects;
};


#endif // SimTK_SIMBODY_VTK_DECORATIVE_GEOMETRY_REP_H_
