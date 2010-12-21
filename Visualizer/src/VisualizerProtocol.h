#ifndef SimTK_SIMBODY_VISUALIZATION_PROTOCOL_H_
#define SimTK_SIMBODY_VISUALIZATION_PROTOCOL_H_

/* -------------------------------------------------------------------------- *
 *                             SimTK Simbody(tm)                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                                             *
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
#include "simbody/internal/Visualizer.h"
#include <pthread.h>
#include <utility>

using namespace SimTK;

/** @file
 * This file defines commands that are used for communication between the 
 * simulation application and the visualization GUI.
 */

// Increment this every time you make *any* change to the protocol;
// we insist on an exact match.
static const unsigned ProtocolVersion   = 26;

// The VisualizerGUI has several predefined cached meshes for common
// shapes so that we don't have to send them. These are the mesh 
// indices for them; they must start with zero.
static const unsigned short MeshBox              = 0;
static const unsigned short MeshEllipsoid        = 1;    // works for sphere
static const unsigned short MeshCylinder         = 2;
static const unsigned short MeshCircle           = 3;

// This serves as the first index number for unique meshes that are 
// defined during this run.
static const unsigned short NumPredefinedMeshes  = 4;

// Commands sent to the GUI.

// This should always be command #1 so we can reliably check whether
// we're talking to a compatible protocol.
static const unsigned char StartupHandshake      = 1;

static const unsigned char StartOfScene          = 2;
static const unsigned char EndOfScene            = 3;
static const unsigned char AddSolidMesh          = 4;
static const unsigned char AddPointMesh          = 5;
static const unsigned char AddWireframeMesh      = 6;
static const unsigned char AddLine               = 7;
static const unsigned char AddText               = 8;
static const unsigned char AddCoords             = 9;
static const unsigned char DefineMesh            = 10;
static const unsigned char DefineMenu            = 11;
static const unsigned char DefineSlider          = 12;
static const unsigned char SetSliderValue        = 13;
static const unsigned char SetSliderRange        = 14;
static const unsigned char SetCamera             = 15;
static const unsigned char ZoomCamera            = 16;
static const unsigned char LookAt                = 17;
static const unsigned char SetFieldOfView        = 18;
static const unsigned char SetClipPlanes         = 19;
static const unsigned char SetSystemUpDirection  = 20;
static const unsigned char SetGroundHeight       = 21;
static const unsigned char SetWindowTitle        = 22;
static const unsigned char SetMaxFrameRate       = 23;
static const unsigned char SetBackgroundColor    = 24;
static const unsigned char SetShowShadows        = 25;
static const unsigned char SetBackgroundType     = 26;



// Events sent from the GUI back to the application.

// This should always be command #1 so we can reliably check whether
// we're talking to a compatible protocol.
static const unsigned char ReturnHandshake       = 1;

static const unsigned char KeyPressed            = 2;
static const unsigned char MenuSelected          = 3;
static const unsigned char SliderMoved           = 4;

class VisualizerProtocol {
public:
    VisualizerProtocol(Visualizer& visualizer);
    void shakeHandsWithGUI(int toGUIPipe, int fromGUIPipe);
    void beginScene(Real simTime);
    void finishScene();
    void drawBox(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawCylinder(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawCircle(const Transform& transform, const Vec3& scale, const Vec4& color, int representation);
    void drawPolygonalMesh(const PolygonalMesh& mesh, const Transform& transform, Real scale, const Vec4& color, int representation);
    void drawLine(const Vec3& end1, const Vec3& end2, const Vec4& color, Real thickness);
    void drawText(const Vec3& position, Real scale, const Vec4& color, const std::string& string);
    void drawCoords(const Transform& transform, Real axisLength, const Vec4& color);
    
    void addMenu(const String& title, int id, const Array_<std::pair<String, int> >& items);
    void addSlider(const String& title, int id, Real min, Real max, Real value);
    void setSliderValue(int id, Real newValue) const;
    void setSliderRange(int id, Real newMin, Real newMax) const;
    
    void setSystemUpDirection(const CoordinateDirection& upDir);
    void setGroundHeight(Real height);
    void setWindowTitle(const String& title) const;
    void setMaxFrameRate(Real rateInFPS) const;
    void setBackgroundColor(const Vec3& color) const;
    void setShowShadows(bool showShadows) const;
    void setBackgroundType(Visualizer::BackgroundType type) const;
    void setCameraTransform(const Transform& transform) const;
    void zoomCamera() const;
    void lookAt(const Vec3& point, const Vec3& upDirection) const;
    void setFieldOfView(Real fov) const;
    void setClippingPlanes(Real near, Real far) const;
private:
    void drawMesh(const Transform& transform, const Vec3& scale, const Vec4& color, 
                  short representation, unsigned short meshIndex);
    int outPipe;

    // For user-defined meshes, map their unique memory addresses to the 
    // assigned VisualizerGUI cache index.
    mutable std::map<const void*, unsigned short> meshes;
    mutable pthread_mutex_t sceneLock;
};


#endif // SimTK_SIMBODY_VISUALIZATION_PROTOCOL_H_
