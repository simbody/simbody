#ifndef SimTK_SIMBODY_VISUALIZATION_PROTOCOL_H_
#define SimTK_SIMBODY_VISUALIZATION_PROTOCOL_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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

#include "simbody/internal/common.h"
#include <pthread.h>
#include <utility>

/** @file
 * This file defines commands that are used for communication between the main application
 * and the visualization GUI.
 */

// Commands sent to the GUI.

static const char START_OF_SCENE = 0;
static const char END_OF_SCENE = 1;
static const char ADD_SOLID_MESH = 2;
static const char ADD_POINT_MESH = 3;
static const char ADD_WIREFRAME_MESH = 4;
static const char ADD_LINE = 5;
static const char ADD_TEXT = 6;
static const char ADD_COORDS = 7;
static const char DEFINE_MESH = 8;
static const char DEFINE_MENU = 9;
static const char DEFINE_SLIDER = 10;
static const char SET_CAMERA = 11;
static const char ZOOM_CAMERA = 12;
static const char LOOK_AT = 13;
static const char SET_FIELD_OF_VIEW = 14;
static const char SET_CLIP_PLANES = 15;
static const char SET_GROUND_POSITION = 16;

// Events sent from the GUI back to the application.

static const char KEY_PRESSED = 0;
static const char MENU_SELECTED = 1;
static const char SLIDER_MOVED = 2;

namespace SimTK {

class Visualizer;

class VisualizationProtocol {
public:
    VisualizationProtocol(Visualizer& visualizer, const String& title);
    void beginScene();
    void finishScene();
    void drawBox(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawCylinder(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawCircle(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawPolygonalMesh(const PolygonalMesh& mesh, const Transform& transform, Real scale, const Vec4& color, int representation);
    void drawLine(const Vec3& end1, const Vec3& end2, const Vec4& color, Real thickness);
    void drawText(const Vec3& position, Real scale, const Vec4& color, const std::string& string);
    void drawCoords(const Transform& transform, Real axisLength, const Vec4& color);
    void addMenu(const String& title, const Array_<std::pair<String, int> >& items);
    void addSlider(const String& title, int id, Real min, Real max, Real value);
    void setGroundPosition(const CoordinateAxis& axis, Real height);

    void setCameraTransform(const Transform& transform) const;
    void zoomCamera() const;
    void lookAt(const Vec3& point, const Vec3& upDirection) const;
    void setFieldOfView(Real fov) const;
    void setClippingPlanes(Real near, Real far) const;
private:
    void drawMesh(const Transform& transform, const Vec3& scale, const Vec4& color, 
                  short representation, short meshIndex);
    int outPipe;
    mutable std::map<const void*, int> meshes;
    mutable pthread_mutex_t sceneLock;
};

}

#endif // SimTK_SIMBODY_VISUALIZATION_PROTOCOL_H_
