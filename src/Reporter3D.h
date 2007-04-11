#ifndef SimTK_SIMBODY_REPORTER_3D_H_
#define SimTK_SIMBODY_REPORTER_3D_H_

/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
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


/** @file
 * This is the user-visible handle class for the VTK Reporter which
 * provides a 3d animation window for viewing Mechanical System
 * simulations.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/DecorativeGeometry.h"


#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace SimTK {

enum {                            
    drawPoints     =  0, // only draw  points at polygon vertices
    drawWireFrame  =  1, // only draw lines along polygon edges
    drawSurface    =  2  // draw shaded polygons
};


class DecorativeGeometry;


class SimTK_SIMBODY_EXPORT Reporter3DGeom {
  public:

    virtual ~Reporter3DGeom() {}
    virtual void* getReporterPolyData() = 0;

};

class SimTK_SIMBODY_EXPORT Reporter3D {
public:
    virtual Reporter3DGeom* generateLineGeometry(     const DecorativeGeometry&, const Vec3&, const Vec3& ) = 0;
    virtual Reporter3DGeom* generateBrickGeometry(    const DecorativeGeometry&, const Vec3& )              = 0;
    virtual Reporter3DGeom* generateCylinderGeometry( const DecorativeGeometry&, Real, Real)                = 0;
    virtual Reporter3DGeom* generateCircleGeometry(   const DecorativeGeometry&, Real)                      = 0; 
    virtual Reporter3DGeom* generateSphereGeometry(   const DecorativeGeometry&, Real)                      = 0;
    virtual Reporter3DGeom* generateFrameGeometry(    const DecorativeGeometry&, Real)                      = 0;

};


} // namespace SimTK

#endif // SimTK_SIMBODY_REPORTER_3D_H_
