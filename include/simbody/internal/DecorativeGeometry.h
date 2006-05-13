#ifndef SimTK_SIMBODY_DECORATIVE_GEOMETRY_H_
#define SimTK_SIMBODY_DECORATIVE_GEOMETRY_H_

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
 * This is the client-side interface to decorative geometry suitable for
 * visualization. DO NOT confuse this with AnalyticGeometry which can represent
 * physically meaningful objects. However, any AnalyticGeometry object can
 * generate some decorative geometry for visualization.
 *
 * Each object has its own local coordinate system and is defined self-consistently
 * but independent of anything else. Clients can associate these with a Transform
 * to produce the effect of placing DecorativeGeometry objects in a scene. We support
 * both 3D objects which are attached to actors in the scene, and 2D "screen"
 * objects like titles which are attached to the display rather than the actors.
 * The classes here deal only with the local-frame definitions of the geometric
 * objects, not their placement in the scene.
 */

#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

class AnalyticGeometry;

/// This is a handle class using the PIMPL design pattern to hide the private
/// implementation. This is effectively an abstract class although the virtual
/// function table is hidden in the private part.
class DecorativeGeometry {
public:
    DecorativeGeometry() : rep(0) { }
    ~DecorativeGeometry();
    DecorativeGeometry(const DecorativeGeometry&);
    DecorativeGeometry& operator=(const DecorativeGeometry&);

    /// Implicit conversion
    DecorativeGeometry(const AnalyticGeometry&);

    void setColor(const Vec3& rgb); // 0-1 for each color; default is 0,0,0 (black)
    void setOpacity(Real);          // 0-1; default is 1 (opaque)
    void setResolution(int);        // some kind of level of detail knob
    void setScale(Real);            // default is 1
 
    /// Is this handle the owner of this rep?
    bool isOwnerHandle() const;
protected:
    class DecorativeGeometryRep* rep;
};

class DecorativeLine : public DecorativeGeometry {
public:
    DecorativeLine() { }
    DecorativeLine(Real length);

    SimTK_PIMPL_DOWNCAST(DecorativeLine, DecorativeGeometry);
};

class DecorativeCircle : public DecorativeGeometry {
public:
    DecorativeCircle() { }
    DecorativeCircle(Real radius);

    SimTK_PIMPL_DOWNCAST(DecorativeCircle, DecorativeGeometry);
};

class DecorativeSphere : public DecorativeGeometry {
public:
    DecorativeSphere() { }
    DecorativeSphere(Real radius);

    SimTK_PIMPL_DOWNCAST(DecorativeSphere, DecorativeGeometry);
};

class DecorativeBrick : public DecorativeGeometry {
public:
    DecorativeBrick() { }
    DecorativeBrick(const Vec3& xyzLengths);

    SimTK_PIMPL_DOWNCAST(DecorativeBrick, DecorativeGeometry);
};

class DecorativeFrame : public DecorativeGeometry {
public:
    DecorativeFrame() { }
    DecorativeFrame(Real axisLength);

    SimTK_PIMPL_DOWNCAST(DecorativeFrame, DecorativeGeometry);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DECORATIVE_GEOMETRY_H_
