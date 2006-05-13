#ifndef SimTK_SIMBODY_VTK_REPORTER_H_
#define SimTK_SIMBODY_VTK_REPORTER_H_

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


/** @file
 * This is the user-visible handle class for the VTK Reporter which
 * provides a 3d animation window for viewing Mechanical System
 * simulations.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"
#include "simbody/internal/DecorativeGeometry.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

class vtkActor;
class vtkRenderWindow;
class vtkRenderer;
class vtkPolyDataMapper;
class vtkProp3D;


namespace SimTK {
    
typedef std::pair<vtkProp3D*, Transform> BodyActor;
typedef std::vector<BodyActor>           ActorList;

class VTKDecoration;

class VTKReporter {
public:
    void addDecoration(int bodyNum, VTKDecoration& d, Transform X_GD);
    void addActor(int bodyNum, vtkActor* a, Transform X_GA);
    ~VTKReporter();
    VTKReporter(const MultibodySystem& m);

    void report(const State& s);
private:

    const MultibodySystem& mbs;
    std::vector<ActorList> bodies;
    vtkRenderWindow* renWin;
    vtkRenderer*     renderer;

    vtkPolyDataMapper *sphereMapper;
    vtkPolyDataMapper *cubeMapper;
    vtkPolyDataMapper *lineMapper;
    vtkPolyDataMapper *cylinderMapper;
    vtkPolyDataMapper *axesMapper;

    void zeroPointers();
    void deletePointers();
    void makeShapes();
    void setConfiguration(int bodyNum, const Transform& X_GB);

};

} // namespace SimTK

#endif // SimTK_SIMBODY_VTK_REPORTER_H_
