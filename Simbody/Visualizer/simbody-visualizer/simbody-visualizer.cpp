/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/Visualizer_InputListener.h"
#include "../src/VisualizerProtocol.h"
#include "lodepng.h"

#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <thread>
#include <condition_variable>
#include <sys/stat.h>
#ifdef _WIN32
    #include <direct.h>
#endif

// Get gl and glut using the appropriate platform-dependent incantations.
#if defined(__APPLE__)
    // OSX comes with a glut implementation. In OSX 10.9 and 10.10, this
    // glut is deprecated and emits deprecation warnings.
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    #include <GLUT/glut.h>
#elif defined(_WIN32)
    #include "glut32/glut.h"    // we have our own private headers
    #include "glut32/glext.h"

    // A Windows-only extension for disabling vsync, allowing unreasonably
    // high frame rates.
    PFNWGLSWAPINTERVALFARPROC wglSwapIntervalEXT;

    // These will hold the dynamically-determined function addresses.

    // These functions are needed for basic Visualizer functionality.
    PFNGLGENBUFFERSPROC glGenBuffers;
    PFNGLBINDBUFFERPROC glBindBuffer;
    PFNGLBUFFERDATAPROC glBufferData;
    PFNGLACTIVETEXTUREPROC glActiveTexture;

    // These are needed only for saving images and movies.
    // Use old EXT names for these so we only require OpenGL 2.0.
    PFNGLGENFRAMEBUFFERSEXTPROC glGenFramebuffersEXT;
    PFNGLGENRENDERBUFFERSEXTPROC glGenRenderbuffersEXT;
    PFNGLBINDFRAMEBUFFEREXTPROC glBindFramebufferEXT;
    PFNGLBINDRENDERBUFFEREXTPROC glBindRenderbufferEXT;
    PFNGLRENDERBUFFERSTORAGEEXTPROC glRenderbufferStorageEXT;
    PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC glFramebufferRenderbufferEXT;
    PFNGLDELETERENDERBUFFERSEXTPROC glDeleteRenderbuffersEXT;
    PFNGLDELETEFRAMEBUFFERSEXTPROC glDeleteFramebuffersEXT;

    // see initGlextFuncPointerIfNeeded() at end of this file
#else
    // Linux: assume we have a good OpenGL 2.0 and working glut or freeglut.
    #define GL_GLEXT_PROTOTYPES
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glext.h>
    #include <GL/glut.h>
#endif

// Returns true if we were able to find sufficient OpenGL functionality to 
// operate. We'll still limp along if we can't get enough to save images.
static bool initGlextFuncPointersIfNeeded(bool& canSaveImages);
static void redrawDisplay();
static void setKeepAlive(bool enable);
static void setVsync(bool enable);
static void shutdown();

// Next, get the functions necessary for reading from and writing to pipes.
#ifdef _WIN32
    #include <io.h>
    #define READ _read
    #define WRITEFUNC _write
    #define CLOSE _close
#else
    #include <unistd.h>
    #define READ read
    #define WRITEFUNC write
    #define CLOSE close
#endif

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strerror, sprintf, etc.
#endif

// gcc 4.4.3 complains bitterly if you don't check the return
// status from the write() system call. This avoids those
// warnings and maybe, someday, will catch an error.
#define WRITE(pipeno, buf, len) \
   {int status=WRITEFUNC((pipeno), (buf), (len)); \
    SimTK_ERRCHK4_ALWAYS(status!=-1, "simbody-visualizer",  \
    "An attempt to write() %d bytes to pipe %d failed with errno=%d (%s).", \
    (len),(pipeno),errno,strerror(errno));}

using namespace SimTK;
using namespace std;

// This is the transform giving the pose of the camera's local frame in the
// model's ground frame. The camera local frame has Y as the up direction,
// -Z as the "look at" direction, and X to the right. We can't know a good
// default transform for the camera until we know what the SimTK::System
// we're viewing considers to be its "up" and "look at" directions.
static fTransform X_GC;

// Communication with the simulator.
static int inPipe, outPipe;
// If the simulator told us to stop communication, then we close the outPipe
// and can no longer write to the simulator.
static std::atomic<bool> writeToSimulator{true};

static void computeBoundingSphereForVertices(const vector<float>& vertices, float& radius, fVec3& center) {
    fVec3 lower(vertices[0], vertices[1], vertices[2]);
    fVec3 upper = lower;
    for (int i = 3; i < (int) vertices.size(); i += 3) {
        for (int j = 0; j < 3; j++) {
            lower[j] = min(lower[j], vertices[i+j]);
            upper[j] = max(upper[j], vertices[i+j]);
        }
    }
    center = (lower+upper)/2;
    float rad2 = 0;
    for (int i = 0; i < (int) vertices.size(); i += 3) {
        float x = center[0]-vertices[i];
        float y = center[1]-vertices[i+1];
        float z = center[2]-vertices[i+2];
        float norm2 = x*x+y*y+z*z;
        if (norm2 > rad2)
            rad2 = norm2;
    }
    radius = sqrt(rad2);
}

class Mesh {
public:
    Mesh(vector<float>& vertices, vector<float>& normals, vector<GLushort>& faces) 
    :   numVertices((int)(vertices.size()/3)), faces(faces) {
        // Build OpenGL buffers.

        GLuint buffers[2];
        glGenBuffers(2, buffers);
        vertBuffer = buffers[0];
        normBuffer = buffers[1];
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), &vertices[0], GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glBufferData(GL_ARRAY_BUFFER, normals.size()*sizeof(float), &normals[0], GL_STATIC_DRAW);

        // Create the list of edges.

        set<pair<GLushort, GLushort> > edgeSet;
        for (int i = 0; i < (int) faces.size(); i += 3) {
            GLushort v1 = faces[i];
            GLushort v2 = faces[i+1];
            GLushort v3 = faces[i+2];
            edgeSet.insert(make_pair(min(v1, v2), max(v1, v2)));
            edgeSet.insert(make_pair(min(v2, v3), max(v2, v3)));
            edgeSet.insert(make_pair(min(v3, v1), max(v3, v1)));
        }
        for (set<pair<GLushort, GLushort> >::const_iterator iter = edgeSet.begin(); iter != edgeSet.end(); ++iter) {
            edges.push_back(iter->first);
            edges.push_back(iter->second);
        }

        // Compute the center and radius.

        computeBoundingSphereForVertices(vertices, radius, center);
    }
    void draw(short representation) const {
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glNormalPointer(GL_FLOAT, 0, 0);
        if (representation == DecorativeGeometry::DrawSurface)
            glDrawElements(GL_TRIANGLES, (GLsizei)faces.size(), GL_UNSIGNED_SHORT, &faces[0]);
        else if (representation == DecorativeGeometry::DrawPoints)
            glDrawArrays(GL_POINTS, 0, numVertices);
        else if (representation == DecorativeGeometry::DrawWireframe)
            glDrawElements(GL_LINES, (GLsizei)edges.size(), GL_UNSIGNED_SHORT, &edges[0]);
    }
    void getBoundingSphere(float& radius, fVec3& center) {
        radius = this->radius;
        center = this->center;
    }
private:
    int numVertices;
    GLuint vertBuffer, normBuffer;
    vector<GLushort> edges, faces;
    fVec3 center;
    float radius;
};

static vector<vector<Mesh*> > meshes;

class RenderedMesh {
public:
    RenderedMesh(const fTransform& transform, const fVec3& scale, const fVec4& color, short representation, unsigned short meshIndex, unsigned short resolution) :
            transform(transform), scale(scale), representation(representation), meshIndex(meshIndex), resolution(resolution) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
        this->color[3] = color[3];
    }
    void draw(bool setColor = true) {
        glPushMatrix();
        glTranslated(transform.p()[0], transform.p()[1], transform.p()[2]);
        fVec4 rot = transform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        if (setColor) {
            if (representation == DecorativeGeometry::DrawSurface)
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
            else
                glColor3fv(color);
        }
        meshes[meshIndex][resolution]->draw(representation);
        glPopMatrix();
    }
    const fTransform& getTransform() const {
        return transform;
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
        meshes[meshIndex][resolution]->getBoundingSphere(radius, center);
        center += transform.p();
        radius *= max(abs(scale[0]), max(abs(scale[1]), abs(scale[2])));
    }
private:
    fTransform transform;
    fVec3 scale;
    GLfloat color[4];
    short representation;
    unsigned short meshIndex, resolution;
};

class RenderedLine {
public:
    RenderedLine(const fVec3& color, float thickness)
    :   color(color), thickness(thickness) {}

    void draw(bool setColor = true) {
        if (setColor)
            glColor3d(color[0], color[1], color[2]);
        glLineWidth(thickness);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glVertexPointer(3, GL_FLOAT, 0, &lines[0]);
        glDrawArrays(GL_LINES, 0, (GLsizei)(lines.size()/3));
    }
    vector<GLfloat>& getLines() {
        return lines;
    }
    const fVec3& getColor() const {
        return color;
    }
    float getThickness() const {
        return thickness;
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
        computeBoundingSphereForVertices(lines, radius, center);
    }
private:
    fVec3 color;
    float thickness;
    vector<GLfloat> lines;
};

class RenderedText {
public:
    RenderedText(const fTransform& X_GT, const fVec3& scale, const fVec3& color, 
                 const string& text, bool faceCamera = true) 
    :   X_GT(X_GT), scale(scale/119), text(text),
        faceCamera(faceCamera) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
    }
    void draw() {
        glPushMatrix();
        glTranslated(X_GT.p()[0], X_GT.p()[1], X_GT.p()[2]);
        const fVec4 rot = faceCamera ? X_GC.R().convertRotationToAngleAxis()
                                     : X_GT.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        glColor3fv(color);
        for (int i = 0; i < (int) text.size(); i++)
            glutStrokeCharacter(GLUT_STROKE_ROMAN, text[i]);
        glPopMatrix();
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
        center = X_GT.p();
        radius = glutStrokeLength(GLUT_STROKE_ROMAN, 
                                  (unsigned char*)text.c_str())*scale[0];
    }
private:
    fTransform X_GT;
    fVec3 scale;
    GLfloat color[3];
    string text;
    bool faceCamera;
};

class ScreenText {
public:
    ScreenText(const string& txt) 
    :   text(txt) {}

    const string& getString() const {return text;}

private:
    string text;
};



/*==============================================================================
                                SCENE DATA
================================================================================
We keep two scenes at a time -- a "front" scene that is currently being
rendered, and a "back" scene is being filled in by the listener thread in
parallel, while the front scene is being rendered. Once the back scene is
complete, the listener thread must block until the renderer has drawn the
front scene (at least once). After the front scene is rendered, the listener
thread replaces it with the back scene. When there is no back scene, the front
scene will be rendered repeatedly whenever an event occurs that may require
that, such as a user moving the camera.

There is a mutex lock for the front scene. It is held by the renderer when
drawing, by any operations (such as camera motion) that would require redrawing
the front scene, and by the listener thread when it wants to swap in the back
scene. There is a condition variable that the listener can wait on when it is
done with the back scene but front scene rendering has not yet finished; that
is signaled by the renderer when it finishes drawing the front scene. */

// This object holds a scene. There are at most two of these around.
class Scene {
public:
    Scene() : simTime(0), sceneHasBeenDrawn(false) {}

    float simTime; // simulated time associated with this frame

    vector<RenderedMesh> drawnMeshes;
    vector<RenderedMesh> solidMeshes;
    vector<RenderedMesh> transparentMeshes;
    vector<RenderedLine> lines;
    vector<RenderedText> sceneText;
    vector<ScreenText>   screenText;

    bool sceneHasBeenDrawn;
};

static Scene* scene=0;      // This is the front scene.

// This mutex must be locked whenever the front scene is being
// modified or drawn, or when it is being swapped with the back scene.
static std::mutex sceneMutex;

// Wait on this if you want to be notified when the scene has been drawn.
static std::condition_variable  sceneHasBeenDrawn;

// This is how long it's been since we've done a redisplay for any reason.
static double lastRedisplayDone = 0; // real time

// When this is true it means something has happened that requires redisplay,
// but we hope to see it when the next scene is rendered rather than initiate
// rendering now. The idle function checks this and initiates rendering if
// no one else does.
static bool passiveRedisplayRequested = true;
// This is reset to zero every time we do an active display.
static int numPassiveRedisplaysSinceLastActive = 0;
// This is reset to zero whenever we post a passive or active display.
static int numMopUpDisplaysSinceLastRedisplay = 0;

// This is the rate that the Visualizer class (on the simulation side) is
// trying to achieve. Normally the GUI simply renders frames whenever they
// are delivered; it does not attempt to time them locally. However, when
// no frames are received, the GUI may issue its own redisplays and should
// not do so any faster than this.
static float maxFrameRate = 30; // in frames/sec

// This is set during the initial handshake with the simulator process.
// It is the file name (not full path) of the simulator executable,
// suitable for use in the title or an "about" message.
static string simulatorExecutableName;

// This is the Simbody version number from the simulator, as major/minor/patch,
// as received in the initial handshake.
static int simbodyVersion[3];
// Constructed from version numbers -- if patch is 0 will just be major.minor.
static string simbodyVersionStr;

// These are used when saving a movie.
static bool savingMovie = false, saveNextFrameToMovie = false;
static bool canSaveImages = false; // is this OpenGL version up to the job?
static string movieDir;
static int movieFrame;
static ParallelWorkQueue imageSaverQueue(5);
static void writeImage(const string& filename);

static void forceActiveRedisplay() {
    passiveRedisplayRequested = false; // cancel if pending
    numPassiveRedisplaysSinceLastActive = 0;
    numMopUpDisplaysSinceLastRedisplay = 0;
    glutPostRedisplay();
}

static void requestPassiveRedisplay() {
    passiveRedisplayRequested = true;
}

static void forcePassiveRedisplay() {
    passiveRedisplayRequested = false; // cancel if pending
    ++numPassiveRedisplaysSinceLastActive;
    numMopUpDisplaysSinceLastRedisplay = 0;
    glutPostRedisplay();
}

static void forceMopUpRedisplay() {
    passiveRedisplayRequested = false; // cancel if pending
    ++numMopUpDisplaysSinceLastRedisplay;
    glutPostRedisplay();
}

// Some commands received by the listener must be executed on the rendering
// thread. They are saved in concrete objects derived from this abstract
// class.
class PendingCommand {
public:
    virtual ~PendingCommand() {}
    virtual void execute() = 0;
};

static int viewWidth, viewHeight;
static GLfloat fieldOfView = GLfloat(SimTK_PI/4);
static GLfloat nearClip = 1;
static GLfloat farClip = 1000;
static GLfloat groundHeight = 0;
static CoordinateDirection groundNormal = YAxis; // the +Y direction
static bool showGround=true, showShadows=true;
static bool showFPS=false, showSimTime=false, showFrameNum=false;
static fVec3 backgroundColor = fVec3(1,1,1); // white is the default
static vector<PendingCommand*> pendingCommands;
static float fps = 0.0f;
static float lastSceneSimTime = 0.0f;
static int fpsCounter = 0, nextMeshIndex;
static int frameCounter = 0;
static double fpsBaseTime = 0;

static vector<string> overlayMessageLines;
static bool displayOverlayMessage = false;
static const double OverlayDisplayTimeInSec = 5; // Leave up for 5 s.
static double overlayStartTime = NaN; // set when overlay is requested
static double overlayEndTime   = NaN; // set when overlay is requested

static void setClearColorToBackgroundColor() {
    glClearColor(backgroundColor[0],backgroundColor[1],backgroundColor[2],1);
}

class PendingMesh : public PendingCommand {
public:
    PendingMesh() {
        index = nextMeshIndex++;
    }
    void execute() override {
        if ((int) meshes.size() <= index)
            meshes.resize(index+1);
        meshes[index].push_back(new Mesh(vertices, normals, faces));
    }
    vector<float> vertices;
    vector<float> normals;
    vector<GLushort> faces;
    int index;
};


static void addVec(vector<float>& data, float x, float y, float z) {
    data.push_back(x);
    data.push_back(y);
    data.push_back(z);
}

static void addVec(vector<unsigned short>& data, int x, int y, int z) {
    data.push_back((unsigned short) x);
    data.push_back((unsigned short) y);
    data.push_back((unsigned short) z);
}

static Mesh* makeBox()  {
    const float halfx = 1;
    const float halfy = 1;
    const float halfz = 1;
    vector<GLfloat> vertices;
    vector<GLfloat> normals;
    vector<GLushort> faces;

    // lower x face
    addVec(vertices, -halfx, -halfy, -halfz);
    addVec(vertices, -halfx, -halfy, halfz);
    addVec(vertices, -halfx, halfy, halfz);
    addVec(vertices, -halfx, halfy, -halfz);
    addVec(normals, -1, 0, 0);
    addVec(normals, -1, 0, 0);
    addVec(normals, -1, 0, 0);
    addVec(normals, -1, 0, 0);
    addVec(faces, 0, 1, 2);
    addVec(faces, 2, 3, 0);
    // upper x face
    addVec(vertices, halfx, halfy, halfz);
    addVec(vertices, halfx, -halfy, halfz);
    addVec(vertices, halfx, -halfy, -halfz);
    addVec(vertices, halfx, halfy, -halfz);
    addVec(normals, 1, 0, 0);
    addVec(normals, 1, 0, 0);
    addVec(normals, 1, 0, 0);
    addVec(normals, 1, 0, 0);
    addVec(faces, 4, 5, 6);
    addVec(faces, 6, 7, 4);
    // lower y face
    addVec(vertices, -halfx, -halfy, -halfz);
    addVec(vertices, halfx, -halfy, -halfz);
    addVec(vertices, halfx, -halfy, halfz);
    addVec(vertices, -halfx, -halfy, halfz);
    addVec(normals, 0, -1, 0);
    addVec(normals, 0, -1, 0);
    addVec(normals, 0, -1, 0);
    addVec(normals, 0, -1, 0);
    addVec(faces, 8, 9, 10);
    addVec(faces, 10, 11, 8);
    // upper y face
    addVec(vertices, halfx, halfy, halfz);
    addVec(vertices, halfx, halfy, -halfz);
    addVec(vertices, -halfx, halfy, -halfz);
    addVec(vertices, -halfx, halfy, halfz);
    addVec(normals, 0, 1, 0);
    addVec(normals, 0, 1, 0);
    addVec(normals, 0, 1, 0);
    addVec(normals, 0, 1, 0);
    addVec(faces, 12, 13, 14);
    addVec(faces, 14, 15, 12);
    // lower z face
    addVec(vertices, -halfx, -halfy, -halfz);
    addVec(vertices, -halfx, halfy, -halfz);
    addVec(vertices, halfx, halfy, -halfz);
    addVec(vertices, halfx, -halfy, -halfz);
    addVec(normals, 0, 0, -1);
    addVec(normals, 0, 0, -1);
    addVec(normals, 0, 0, -1);
    addVec(normals, 0, 0, -1);
    addVec(faces, 16, 17, 18);
    addVec(faces, 18, 19, 16);
    // upper z face
    addVec(vertices, halfx, halfy, halfz);
    addVec(vertices, -halfx, halfy, halfz);
    addVec(vertices, -halfx, -halfy, halfz);
    addVec(vertices, halfx, -halfy, halfz);
    addVec(normals, 0, 0, 1);
    addVec(normals, 0, 0, 1);
    addVec(normals, 0, 0, 1);
    addVec(normals, 0, 0, 1);
    addVec(faces, 20, 21, 22);
    addVec(faces, 22, 23, 20);
    return new Mesh(vertices, normals, faces);
}

static Mesh* makeSphere(unsigned short resolution) {
    const int numLatitude = 4*resolution;
    const int numLongitude = 6*resolution;
    const float radius = 1.0f;
    vector<GLfloat> vertices;
    vector<GLfloat> normals;
    vector<GLushort> faces;
    addVec(vertices, 0, radius, 0);
    addVec(normals, 0, 1, 0);
    for (int i = 0; i < numLatitude; i++) {
        float phi = (float) (((i+1)*SimTK_PI)/(numLatitude+1));
        float sphi = sin(phi);
        float cphi = cos(phi);
        float y = radius*cphi;
        float r = radius*sphi;
        for (int j = 0; j < numLongitude; j++) {
            float theta = (float) ((j*2*SimTK_PI)/numLongitude);
            float stheta = sin(theta);
            float ctheta = cos(theta);
            addVec(vertices, r*ctheta, y, r*stheta);
            addVec(normals, sphi*ctheta, cphi, sphi*stheta);
        }
    }
    addVec(vertices, 0, -radius, 0);
    addVec(normals, 0, -1, 0);
    for (int i = 1; i < numLongitude; i++)
        addVec(faces, 0, i+1, i);
    addVec(faces, 0, 1, numLongitude);
    for (int i = 1; i < numLatitude; i++) {
        int base = (i-1)*numLongitude+1;
        for (int j = 0; j < numLongitude; j++) {
            int v1 = base+j;
            int v2 = (j == numLongitude-1 ? base : v1+1);
            int v3 = v1+numLongitude;
            int v4 = v2+numLongitude;
            addVec(faces, v1, v4, v3);
            addVec(faces, v1, v2, v4);
        }
    }
    int first = (numLatitude-1)*numLongitude+1;
    int last = numLatitude*numLongitude+1;
    for (int i = first; i < last-1; i++)
        addVec(faces, i, i+1, last);
    addVec(faces, last-1, first, last);
    return new Mesh(vertices, normals, faces);
}

static Mesh* makeCylinder(unsigned short resolution) {
    const int numSides = 6*resolution;
    const float halfHeight = 1;
    const float radius = 1;
    vector<GLfloat> vertices;
    vector<GLfloat> normals;
    vector<GLushort> faces;

    // Create the top face.

    addVec(vertices, 0, halfHeight, 0);
    addVec(normals, 0, 1.0, 0);
    for (int i = 0; i < numSides; i++) {
        float theta = (float) ((i*2*SimTK_PI)/numSides);
        float stheta = sin(theta);
        float ctheta = cos(theta);
        addVec(vertices, radius*ctheta, halfHeight, radius*stheta);
        addVec(normals, 0, 1.0, 0);
    }
    for (int i = 1; i < numSides; i++)
        addVec(faces, i, 0, i+1);
    addVec(faces, numSides, 0, 1);

    // Create the bottom face.

    int bottomStart = numSides+1;
    addVec(vertices, 0, -halfHeight, 0);
    addVec(normals, 0, -1.0, 0);
    for (int i = 0; i < numSides; i++) {
        float theta = (float) ((i*2*SimTK_PI)/numSides);
        float stheta = sin(theta);
        float ctheta = cos(theta);
        addVec(vertices, radius*ctheta, -halfHeight, radius*stheta);
        addVec(normals, 0, -1.0, 0);
    }
    for (int i = 1; i < numSides; i++)
        addVec(faces, bottomStart+i, bottomStart+i+1, bottomStart);
    addVec(faces, bottomStart+numSides, bottomStart+1, bottomStart);

    // Create the sides.

    for (int i = 0; i < numSides; i++) {
        float theta = (float) ((i*2*SimTK_PI)/numSides);
        float stheta = sin(theta);
        float ctheta = cos(theta);
        float x = radius*ctheta;
        float z = radius*stheta;
        addVec(vertices, x, halfHeight, z);
        addVec(normals, ctheta, 0, stheta);
        addVec(vertices, x, -halfHeight, z);
        addVec(normals, ctheta, 0, stheta);
    }
    int sideStart = 2*numSides+2;
    for (int i = 0; i < numSides-1; i++) {
        int base = sideStart+2*i;
        addVec(faces, base, base+2, base+1);
        addVec(faces, base+1, base+2, base+3);
    }
    addVec(faces, sideStart+2*numSides-2, sideStart, sideStart+2*numSides-1);
    addVec(faces, sideStart+2*numSides-1, sideStart, sideStart+1);
    return new Mesh(vertices, normals, faces);
}

static Mesh* makeCircle(unsigned short resolution) {
    const int numSides = 6*resolution;
    const float radius = 1;
    vector<GLfloat> vertices;
    vector<GLfloat> normals;
    vector<GLushort> faces;

    // Create the front face.

    addVec(vertices, 0, 0, 0);
    addVec(normals, 0, 0, -1.0);
    for (int i = 0; i < numSides; i++) {
        float theta = (float) ((i*2*SimTK_PI)/numSides);
        float stheta = sin(theta);
        float ctheta = cos(theta);
        addVec(vertices, radius*ctheta, radius*stheta, 0);
        addVec(normals, 0, 0, -1.0);
    }
    for (int i = 1; i < numSides; i++)
        addVec(faces, i, 0, i+1);
    addVec(faces, numSides, 0, 1);

    // Create the back face.

    addVec(vertices, 0, 0, 0);
    addVec(normals, 0, 0, 1.0);
    for (int i = 0; i < numSides; i++) {
        float theta = (float) ((i*2*SimTK_PI)/numSides);
        float stheta = sin(theta);
        float ctheta = cos(theta);
        addVec(vertices, radius*ctheta, radius*stheta, 0);
        addVec(normals, 0, 0, 1.0);
    }
    int backStart = numSides+1;
    for (int i = 1; i < numSides; i++)
        addVec(faces, backStart, backStart+i, backStart+i+1);
    addVec(faces, backStart, backStart+numSides, backStart+1);
    return new Mesh(vertices, normals, faces);
}

class PendingStandardMesh : public PendingCommand {
public:
    PendingStandardMesh(unsigned short meshIndex, unsigned short resolution) : meshIndex(meshIndex), resolution(resolution) {
    }
    void execute() override {
        if ((int) meshes[meshIndex].size() <= resolution)
            meshes[meshIndex].resize(resolution+1, NULL);
        if (meshes[meshIndex][resolution] == NULL) {
            switch (meshIndex) {
                case MeshBox:
                    meshes[meshIndex][resolution] = makeBox();
                    break;
                case MeshEllipsoid:
                    meshes[meshIndex][resolution] = makeSphere(resolution);
                    break;
                case MeshCylinder:
                    meshes[meshIndex][resolution] = makeCylinder(resolution);
                    break;
                case MeshCircle:
                    meshes[meshIndex][resolution] = makeCircle(resolution);
                    break;
            }
        }
    }
    unsigned short meshIndex, resolution;
};

// Caution -- make sure scene is locked before you call this function.
static void computeSceneBounds(const Scene* scene, float& radius, fVec3& center) {
    // Record the bounding sphere of every object in the scene.

    vector<fVec3> centers;
    vector<float> radii;
    for (int i = 0; i < (int) scene->drawnMeshes.size(); i++) {
        fVec3 center;
        float radius;
        scene->drawnMeshes[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->solidMeshes.size(); i++) {
        fVec3 center;
        float radius;
        scene->solidMeshes[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->transparentMeshes.size(); i++) {
        fVec3 center;
        float radius;
        scene->transparentMeshes[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->lines.size(); i++) {
        fVec3 center;
        float radius;
        scene->lines[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->sceneText.size(); i++) {
        fVec3 center;
        float radius;
        scene->sceneText[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }

    // Find the overall bounding sphere of the scene.

    if (centers.size() == 0) {
        radius = 0;
        center = fVec3(0);
    }
    else {
        fVec3 lower = centers[0]-radii[0];
        fVec3 upper = centers[0]+radii[0];
        for (int i = 1; i < (int) centers.size(); i++) {
            for (int j = 0; j < 3; j++) {
                lower[j] = min(lower[j], centers[i][j]-radii[i]);
                upper[j] = max(upper[j], centers[i][j]+radii[i]);
            }
        }
        center = (lower+upper)/2;
        radius = 0;
        for (int i = 0; i < (int) centers.size(); i++)
            radius = max(radius, (centers[i]-center).norm()+radii[i]);
    }
}

static void zoomCameraToShowWholeScene(bool sceneAlreadyLocked=false) {
    float radius;
    fVec3 center;
    std::unique_lock<std::mutex> lock(sceneMutex, std::defer_lock);
    if (!sceneAlreadyLocked)
        lock.lock();                              //-------- LOCK SCENE --------
    computeSceneBounds(scene, radius, center);
    if (!sceneAlreadyLocked)
        lock.unlock();                            //----- UNLOCK SCENE ---------
   float viewDistance = 
        radius/tan(min(fieldOfView, fieldOfView*viewWidth/viewHeight)/2);
    // Back up 1 unit more to make sure we don't clip at this position.
    // Also add a modest offset in the (x,y) direction to avoid edge-on views.
    const float offset = std::max(1.f, viewDistance/10);
    X_GC.updP() = center+X_GC.R()*fVec3(offset, offset, viewDistance+1);
    // Now fix the aim to point at the center.
    fVec3 zdir = X_GC.p() - center;
    if (zdir.normSqr() >= square(1e-6))
        X_GC.updR().setRotationFromTwoAxes(fUnitVec3(zdir), ZAxis, 
                                           X_GC.y(),        YAxis);
}

class PendingCameraZoom : public PendingCommand {
public:
    void execute() override {
        zoomCameraToShowWholeScene(true); // scene already locked
    }
};

class PendingSetCameraTransform : public PendingCommand {
public:
    PendingSetCameraTransform(fVec3 R, fVec3 p) : Rxyz(R), p(p) { }

    void execute() override {
        X_GC.updR().setRotationToBodyFixedXYZ(Rxyz);
        X_GC.updP() = p;
    }

private:
    fVec3 Rxyz;
    fVec3 p;
};

class PendingWindowTitleChange : public PendingCommand {
public:
    PendingWindowTitleChange(const string& title) : title(title) {}
    void execute() override {glutSetWindowTitle(title.c_str());}
private:
    string title;
};


class PendingBackgroundColorChange : public PendingCommand {
public:
    PendingBackgroundColorChange(const fVec3& color) : color(color) {}
    void execute() override {
        backgroundColor=color; setClearColorToBackgroundColor();
    }
private:
    fVec3 color;
};



/*==============================================================================
                                   MENU
==============================================================================*/

// The glut callback for menu item picking.
static void menuSelected(int option);
static const int InvalidMenuId = std::numeric_limits<int>::min();
class Menu {
public:
    Menu(string title, int id, const vector<pair<string, int> >& items,
         void(*handler)(int))
    :   title(title), menuId(id), items(items), handler(handler),
        hasCreated(false) 
    {   glutId = -1; minx=miny=maxx=maxy= -1; }

    // This is called once, the first time we try to draw this menu.
    void createMenu() {
        glutId = glutCreateMenu(handler);
        vector<string> components;
        vector<int> submenuIds;
        for (int i = 0; i < (int) items.size(); i++) {
            size_t start = 0;
            string spec = items[i].first;
            for (int componentIndex = 0; ; componentIndex++) {
                size_t end = spec.find('/', start);
                while (end < spec.size()-1 && spec[end+1] == '/') {
                    // A double slash.
                    spec.erase(end, 1);
                    end = spec.find('/', end+1);
                    continue;
                }
                string substring = spec.substr(start, end-start);
                if (componentIndex < (int)components.size() && substring != components[componentIndex]) {
                    components.resize(componentIndex);
                    submenuIds.resize(componentIndex);
                }
                if (componentIndex == components.size())
                    components.push_back(substring);
                if (end == string::npos)
                    break;
                start = end+1;
            }
            int firstNewSubmenu = (int) submenuIds.size();
            for (int j = firstNewSubmenu; j < (int) components.size()-1; j++)
                submenuIds.push_back(glutCreateMenu(handler));
            glutSetMenu(firstNewSubmenu == 0 ? glutId
                                             : submenuIds[firstNewSubmenu-1]);
            for (int j = firstNewSubmenu; j < (int) components.size()-1; j++) {
                glutAddSubMenu(components[j].c_str(), submenuIds[j]);
                glutSetMenu(submenuIds[j]);
            }
            glutAddMenuEntry(components[components.size()-1].c_str(),
                             items[i].second);
        }
    }

    int draw(int x, int y) {
        if (!hasCreated) {
            createMenu();
            hasCreated = true;
        }
        minx = x;
        miny = y-18;
        int width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18,
                                     (unsigned char*) title.c_str());
        maxx = x+width+14;
        maxy = y+3;
        glColor3f(0.9f, 0.9f, 0.9f);
        glBegin(GL_POLYGON);
        glVertex2i(minx+2, miny);
        glVertex2i(minx, miny+1);
        glVertex2i(minx, maxy-1);
        glVertex2i(minx+2, maxy);
        glVertex2i(maxx-1, maxy);
        glVertex2i(maxx, maxy-1);
        glVertex2i(maxx, miny+1);
        glVertex2i(maxx-1, miny);
        glEnd();
        glColor3f(0.2f, 0.2f, 0.2f);
        glRasterPos2f(GLfloat(x+2), GLfloat(y));
        for (int i = 0; i < (int) title.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, title[i]);
        glBegin(GL_TRIANGLES);
        glVertex2i(maxx-2, y-8);
        glVertex2i(maxx-10,y-8);
        glVertex2i(maxx-6, y-4);
        glEnd();
        return width+25;
    }

    // When the mouse moves onto this menu's pulldown tab, attach the menu
    // to the mouse's left button. When it is next seen moving outside the
    // pulldown tab, unattach the menu if it is still attached. Only one menu
    // can be displayed at a time.
    void mouseMoved(int x, int y) {
        if (!hasCreated)
            return;
        if (x >= minx && y >= miny && x <= maxx && y <= maxy) {
            if (currentMenuGlutId != glutId) {
                glutSetMenu(glutId);
                glutAttachMenu(GLUT_LEFT_BUTTON);
                currentMenuGlutId = glutId;
                currentMenuId = menuId;
            }
        }
        else if (currentMenuGlutId == glutId) {
            glutSetMenu(glutId);
            glutDetachMenu(GLUT_LEFT_BUTTON);
            currentMenuGlutId = currentMenuId = InvalidMenuId;
        }
    }

    int getId() const {return menuId;}
    const string& getTitle() const {return title;}

    static const int& getCurrentMenuId() {return currentMenuId;}
    static bool isMenuAttached() {return currentMenuId != InvalidMenuId;}

private:
    string                      title;
    int                         menuId; // assigned by creator
    vector<pair<string, int> >  items;
    void                      (*handler)(int);
    int                         glutId; // assigned by glut
    int                         minx, miny, maxx, maxy;
    bool                        hasCreated;

    static int currentMenuGlutId;
    static int currentMenuId;
};

int Menu::currentMenuId     = InvalidMenuId; // user assigned id
int Menu::currentMenuGlutId = InvalidMenuId;

// This is the handler for simulator-defined menus. All it does is send back
// to the simulator the integer index that was associated with the particular
// menu entry on which the user clicked.
void menuSelected(int option) {
    if (writeToSimulator) {
        WRITE(outPipe, &MenuSelected, 1);
        WRITE(outPipe, &Menu::getCurrentMenuId(), sizeof(int));
        WRITE(outPipe, &option, sizeof(int));
    }
}

static vector<Menu>     menus;

// Look up a particular menu; returns null pointer if not there.
// Note that this is the user-assigned id, not the one from glut.
static Menu* findMenuById(int id) {
    for (unsigned i=0; i < menus.size(); ++i)
        if (menus[i].getId() == id)
            return &menus[i];
    return 0;
}



/*==============================================================================
                                 SLIDER
==============================================================================*/
class Slider {
public:
    Slider(string title, int id, float minValue, float maxValue, float value)
    :   title(title), id(id), minValue(minValue), maxValue(maxValue),
        position((value-minValue)/(maxValue-minValue))
    {
        labelWidth = glutBitmapLength(GLUT_BITMAP_HELVETICA_18,
                                      (unsigned char*) title.c_str());
        if (labelWidth > maxLabelWidth)
            maxLabelWidth = labelWidth;

        minx=miny=maxx=maxy=handlex=clickOffset= -1;
        dragging = false;
    }

    int draw(int y) {
        minx = maxLabelWidth+16;
        maxx = minx+sliderWidth+handleWidth;
        miny = y-18;
        maxy = y+3;
        int frameMinx = 10;
        int frameMaxx = maxx+2;

        // Draw the background.

        glColor3f(0.9f, 0.9f, 0.9f);
        glBegin(GL_POLYGON);
        glVertex2i(frameMinx+2, miny);
        glVertex2i(frameMinx, miny+1);
        glVertex2i(frameMinx, maxy-1);
        glVertex2i(frameMinx+2, maxy);
        glVertex2i(frameMaxx-1, maxy);
        glVertex2i(frameMaxx, maxy-1);
        glVertex2i(frameMaxx, miny+1);
        glVertex2i(frameMaxx-1, miny);
        glEnd();

        // Draw the slider.

        glColor3f(0.5f, 0.5f, 0.5f);
        glBegin(GL_QUADS);
        glVertex2i(minx, miny+12);
        glVertex2i(maxx, miny+12);
        glVertex2i(maxx, miny+9);
        glVertex2i(minx, miny+9);
        glEnd();
        glColor3f(0.2f, 0.2f, 0.2f);
        handlex = minx+(int) (sliderWidth*position);
        glBegin(GL_QUADS);
        glVertex2i(handlex+handleWidth, miny+2);
        glVertex2i(handlex, miny+2);
        glVertex2i(handlex, miny+19);
        glVertex2i(handlex+handleWidth, miny+19);
        glEnd();

        // Draw the label.

        glRasterPos2f(GLfloat(12+maxLabelWidth-labelWidth), GLfloat(y));
        for (int i = 0; i < (int) title.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, title[i]);

        // Draw the value.
        char buffer[32];
        sprintf(buffer, "%g", (double)calcValue());
        string value(buffer);
        int valueWidth = glutBitmapLength(GLUT_BITMAP_HELVETICA_10,
                            (unsigned char*) value.c_str());
        int valueMinx = maxx+4;
        int valueMaxx = valueMinx+valueWidth+4;
        int valueMiny = y-14;
        int valueMaxy = y-1;
        glColor3f(0.9f, 0.9f, 0.9f);
        glBegin(GL_POLYGON);
        glVertex2i(valueMinx+2, valueMiny);
        glVertex2i(valueMinx, valueMiny+1);
        glVertex2i(valueMinx, valueMaxy-1);
        glVertex2i(valueMinx+2, valueMaxy);
        glVertex2i(valueMaxx-1, valueMaxy);
        glVertex2i(valueMaxx, valueMaxy-1);
        glVertex2i(valueMaxx, valueMiny+1);
        glVertex2i(valueMaxx-1, valueMiny);
        glEnd();
        glColor3f(0.2f, 0.2f, 0.2f);
        glRasterPos2f(GLfloat(valueMinx+2), GLfloat(y-4));
        for (int i = 0; i < (int) value.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, value[i]);

        return y-25;
    }

    bool mousePressed(int x, int y) {
        dragging = false;
        if (x >= minx && y >= miny && x <= maxx && y <= maxy) {
            if (x < handlex) {
                position = std::max(position-0.1f, 0.0f);
                positionChanged();
            }
            else if (x > handlex+handleWidth) {
                position = std::min(position+0.1f, 1.0f);
                positionChanged();
            }
            else {
                clickOffset = x-handlex;
                dragging = true;
            }
            return true;
        }
        else
            return false;
    }

    void mouseDragged(int x) {
        if (dragging) {
            float oldPosition = position;
            position = clamp(0.f, (x-clickOffset-minx)/(float)sliderWidth, 1.f);
            if (position != oldPosition)
                positionChanged();
        }
    }

    // The position is just 0..1; map that into the actual value.
    float calcValue() const {return minValue+position*(maxValue-minValue);}

    // Use this for programmatic changes to slider value. This does not get
    // reported back to the simulation process.
    void changeValue(float newValue) {
        clampInPlace(minValue, newValue, maxValue); // fix newValue
        const float range = maxValue-minValue;
        if (range > 0) position = clamp(0.f, (newValue-minValue)/range, 1.f);
        else position = 0.5f;
        requestPassiveRedisplay();              //----- PASSIVE REDISPLAY ----
    }

    // This is for a programmatic change to the range after a slider has been
    // defined. Value is unchanged if it still fits, otherwise it moves to
    // the limit. This is not reported back to the simulation process.
    void changeRange(float newMin, float newMax) {
        SimTK_ASSERT_ALWAYS(newMax >= newMin, "Slider::changeRange(): bad range");
        const float curValue = calcValue();
        minValue = newMin; maxValue = newMax; // change bounds
        changeValue(curValue); // recalculate relative position
    }

    int getId() const {return id;}
    const string& getTitle() const {return title;}

private:
    void positionChanged() {
        if (writeToSimulator) {
            WRITE(outPipe, &SliderMoved, 1);
            WRITE(outPipe, &id, sizeof(int));
            float value = calcValue();
            WRITE(outPipe, &value, sizeof(float));
        }
        requestPassiveRedisplay();              //----- PASSIVE REDISPLAY ----
    }
    string title;
    int labelWidth;
    int id, minx, miny, maxx, maxy, handlex;
    int clickOffset;
    bool dragging;
    float minValue, maxValue;   // actual min,max values
    float position;             // always in range [0,1]
    static int maxLabelWidth;
    static const int handleWidth = 5;
    static const int sliderWidth = 100;
};

int Slider::maxLabelWidth = 0;


static vector<Slider>   sliders;

// Look up a particular slider; returns null pointer if not there.
static Slider* findSliderById(int id) {
    for (unsigned i=0; i < sliders.size(); ++i)
        if (sliders[i].getId() == id)
            return &sliders[i];
    return 0;
}



/*==============================================================================
                             GROUND AND SKY
==============================================================================*/
static void drawSkyVertex(fVec3 position, float texture) {
    glTexCoord1f(texture);
    glVertex3fv(&position[0]);
}

static void drawGroundAndSky(float farClipDistance) {
    static bool initialized = false;
    static GLuint skyTexture;
    static GLuint groundTexture;
    if (!initialized) {
        initialized = true;

        // Create a texture to use as the sky.

        glGenTextures(1, &skyTexture);
        glBindTexture(GL_TEXTURE_1D, skyTexture);
        glTexParameterf(GL_TEXTURE_1D, GL_GENERATE_MIPMAP, GL_TRUE);
        const int width = 256;
        float skyImage[3*width];
        for (int i = 0; i < width; i++) {
            float fract = pow(i/(float) width, 1.8f);
            skyImage[3*i] = fract;
            skyImage[3*i+1] = fract;
            skyImage[3*i+2] = 1;
        }
        glTexImage1D(GL_TEXTURE_1D, 0, 3, width, 0, GL_RGB, GL_FLOAT, skyImage);

        // Create a texture to use as the ground.

        srand(0);
        glGenTextures(1, &groundTexture);
        glBindTexture(GL_TEXTURE_2D, groundTexture);
        glTexParameterf(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
        float groundImage[width*width];
        for (int i = 0; i < width; i++) {
            float x = i/(float) width;
            for (int j = 0; j < width; j++) {
                float y = j/(float) width;
                double line = min(min(min(x, y), 1.0f-x), 1.0f-y);
                float noise = (rand()%255)/255.0f-0.5f;
                groundImage[i*width+j] = float(std::pow(line,.1)*(.35f+noise));
            }
        }
        glTexImage2D(GL_TEXTURE_2D, 0, 1, width, width, 0, GL_RED, GL_FLOAT, groundImage);
    }

    // Draw the box to represent the sky.

    float viewDistance = farClipDistance*0.5f;
    fVec3 center = X_GC.p();
    const float sign = (float)groundNormal.getDirection(); // 1 or -1
    const CoordinateAxis axis = groundNormal.getAxis();
    float top = center[axis]+sign*viewDistance;
    center[axis] = 0;
    fVec3 offset1 = viewDistance*(X_GC.R()(2)%fUnitVec3(axis));
    fVec3 offset2 = (offset1%fUnitVec3(axis));
    fVec3 corner1 = center+(-offset1-offset2);
    fVec3 corner2 = center+(offset1-offset2);
    fVec3 corner3 = center+(offset1+offset2);
    fVec3 corner4 = center+(-offset1+offset2);
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, skyTexture);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glBegin(GL_QUAD_STRIP);
    fVec3 offset3 = sign*fUnitVec3(axis);
    drawSkyVertex(corner1+top*offset3, 0);
    drawSkyVertex(corner1+groundHeight*offset3, 1);
    drawSkyVertex(corner2+top*offset3, 0);
    drawSkyVertex(corner2+groundHeight*offset3, 1);
    drawSkyVertex(corner3+top*offset3, 0);
    drawSkyVertex(corner3+groundHeight*offset3, 1);
    drawSkyVertex(corner4+top*offset3, 0);
    drawSkyVertex(corner4+groundHeight*offset3, 1);
    drawSkyVertex(corner1+top*offset3, 0);
    drawSkyVertex(corner1+groundHeight*offset3, 1);
    glEnd();
    glDisable(GL_TEXTURE_1D);
    glColor3f(0, 0, 1);
    glBegin(GL_QUADS);
    drawSkyVertex(corner1+0.99f*top*offset3, 0);
    drawSkyVertex(corner2+0.99f*top*offset3, 0);
    drawSkyVertex(corner3+0.99f*top*offset3, 0);
    drawSkyVertex(corner4+0.99f*top*offset3, 0);
    glEnd();
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    // Draw the ground plane.

    center = X_GC.p()-X_GC.R()*fVec3(0, 0, 0.9999f*farClipDistance);
    center[1] = 0;
    corner1 = center+fVec3(-farClipDistance, 0, -farClipDistance);
    corner2 = center+fVec3(farClipDistance, 0, -farClipDistance);
    corner3 = center+fVec3(farClipDistance, 0, farClipDistance);
    corner4 = center+fVec3(-farClipDistance, 0, farClipDistance);
    // We need to calculate the GL transform T_GP that gives the ground
    // plane's coordinate frame P in the ground frame G.
    Mat<4, 4, GLfloat> T_GP(1);
    fVec3::updAs(&T_GP(YAxis)[0]) = fVec3(fUnitVec3(groundNormal)); // signed
    fVec3::updAs(&T_GP(XAxis)[0]) = fVec3(fUnitVec3(axis.getPreviousAxis()));
    fVec3::updAs(&T_GP(ZAxis)[0]) = fVec3(fUnitVec3(axis.getNextAxis()));
    T_GP[axis][3] = sign*groundHeight;
    if (showShadows) {
        // Draw shadows on the ground. These are drawn in the shadow frame S,
        // which we'll distort slightly from the ground plane P.
        Mat<4, 4, GLfloat> T_PS(1);
        T_PS[axis.getPreviousAxis()][axis] = 0.2f;
        T_PS[axis][axis] = 0.0f;
        T_PS[axis.getNextAxis()][axis] = 0.2f;
        T_PS[axis][3] = sign*groundHeight;
        glPushMatrix();
        glMultMatrixf(&T_PS[0][0]);
        // Solid and transparent shadows are the same color (sorry). Trying to
        // mix light and dark shadows is much harder and any simple attempts
        // (e.g. put light shadows on top of dark ones) look terrible.
        glColor3f(0.3f, 0.2f, 0.0f);
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
            scene->solidMeshes[i].draw(false);
        for (int i = 0; i < (int) scene->transparentMeshes.size(); i++)
            scene->transparentMeshes[i].draw(false);
        for (int i = 0; i < (int) scene->lines.size(); i++)
            scene->lines[i].draw(false);
        glPopMatrix();
    }
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, groundTexture);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glDisable(GL_CULL_FACE);
    glPushMatrix();
    glMultMatrixf(&T_GP[0][0]);
    glDepthRange(0.01, 1.0);
    float color2[] = {1.0f, 0.8f, 0.7f};
    glBindTexture(GL_TEXTURE_2D, groundTexture);
    glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, color2);
    glBegin(GL_QUADS);
    glColor3f(0.3f, 0.2f, 0.0f);
    glTexCoord2d(2.0f*corner1[0], 2.0f*corner1[2]);
    glVertex3f(corner1[0], corner1[1], corner1[2]);
    glTexCoord2d(2.0f*corner2[0], 2.0f*corner2[2]);
    glVertex3f(corner2[0], corner2[1], corner2[2]);
    glTexCoord2d(2.0f*corner3[0], 2.0f*corner3[2]);
    glVertex3f(corner3[0], corner3[1], corner3[2]);
    glTexCoord2d(2.0f*corner4[0], 2.0f*corner4[2]);
    glVertex3f(corner4[0], corner4[1], corner4[2]);
    glEnd();
    glDepthRange(0.0, 1.0);
    glEnable(GL_CULL_FACE);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}



/*==============================================================================
                               RENDER SCENE
==============================================================================*/
static void renderScene(std::vector<std::string>* screenText = NULL) {
    static bool firstTime = true;
    static GLfloat prevNearClip; // initialize to near & farClip
    static GLfloat prevFarClip;
    if (firstTime) {
        prevNearClip = nearClip;
        prevFarClip = farClip;
        firstTime = false;
    }

    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();
    glDisable(GL_LIGHTING);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    std::lock_guard<std::mutex> lock(sceneMutex); //-------- LOCK SCENE --------

    if (scene != NULL) {
        // Remember the simulated time associated with the most recent rendered
        // scene.
        lastSceneSimTime = scene->simTime;

        // Execute any pending commands that need to be executed on the
        // rendering thread. This happens only the first time a particular
        // scene is drawn because we delete the pending commands after
        // executing them.
        for (int i = 0; i < (int) pendingCommands.size(); i++) {
            pendingCommands[i]->execute();
            delete pendingCommands[i];
        }
        pendingCommands.clear();

        // Set up the viewpoint.

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glViewport(0, 0, viewWidth, viewHeight);
        float sceneRadius;
        fVec3 sceneCenter;
        // Scene is already locked so OK to call this.
        computeSceneBounds(scene, sceneRadius, sceneCenter);
        float sceneScale = std::max(0.1f, sceneRadius);

        float centerDepth = ~(X_GC.p()-sceneCenter)*(X_GC.R().col(2));
        float nearClipDistance, farClipDistance;
        if (showGround) {
            nearClipDistance = nearClip;
            const float wantFarClipDistance = min(farClip,
                                                  max(2*sceneScale, 2*(centerDepth+sceneRadius)));
            if (std::abs(wantFarClipDistance-prevFarClip) <= 0.5*prevFarClip)
                farClipDistance = prevFarClip; // hysteresis to avoid jiggling ground
            else farClipDistance = wantFarClipDistance;
        }
        else {
            nearClipDistance = max(nearClip, centerDepth-sceneRadius);
            farClipDistance = min(farClip, centerDepth+sceneRadius);
        }
        prevNearClip = nearClipDistance;
        prevFarClip  = farClipDistance;

        gluPerspective(fieldOfView*SimTK_RADIAN_TO_DEGREE, (GLdouble) viewWidth/viewHeight, nearClipDistance, farClipDistance);
        glMatrixMode(GL_MODELVIEW);
        fVec3 cameraPos = X_GC.p();
        fVec3 centerPos = X_GC.p()+X_GC.R()*fVec3(0, 0, -1);
        fVec3 upDir = X_GC.R()*fVec3(0, 1, 0);
        gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2], centerPos[0], centerPos[1], centerPos[2], upDir[0], upDir[1], upDir[2]);

        // Render the objects in the scene.

        if (showGround)
            drawGroundAndSky(farClipDistance);
        for (int i = 0; i < (int) scene->lines.size(); i++)
            scene->lines[i].draw();
        glLineWidth(2);
        for (int i = 0; i < (int) scene->sceneText.size(); i++)
            scene->sceneText[i].draw();
        glLineWidth(1);
        glEnableClientState(GL_NORMAL_ARRAY);
        for (int i = 0; i < (int) scene->drawnMeshes.size(); i++)
            scene->drawnMeshes[i].draw();
        glEnable(GL_LIGHTING);
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
            scene->solidMeshes[i].draw();
        glEnable(GL_BLEND);
        glDepthMask(GL_FALSE);
        vector<pair<float, int> > order(scene->transparentMeshes.size());
        for (int i = 0; i < (int) order.size(); i++)
            order[i] = make_pair((float)(~X_GC.R()*scene->transparentMeshes[i].getTransform().p())[2], i);
        sort(order.begin(), order.end());
        for (int i = 0; i < (int) order.size(); i++)
            scene->transparentMeshes[order[i].second].draw();

        // Update the frame rate counter.

        const double now = realTime(); // in seconds at high resolution
        const double elapsed = now - fpsBaseTime;
        if (elapsed >= 1) {
            fps = float(fpsCounter/elapsed);
            fpsBaseTime = now;
            fpsCounter = 0;
        }

        // Extract the scene text if requested. This will get displayed
        // later but must be copied out now since the scene will get
        // overwritten once we say it has been drawn.
        if (screenText != NULL)
            for (int i = 0; i < (int)scene->screenText.size(); ++i)
                screenText->push_back(scene->screenText[i].getString());

        scene->sceneHasBeenDrawn = true;
    }

    // Signal in case the listener is waiting.
    sceneHasBeenDrawn.notify_one();               //----- SIGNAL CONDITION -----
                                                  //----- UNLOCK SCENE ---------
}


// Redraw everything currently being displayed. That means:
//  - the "front" scene
//  - menus
//  - overlay text like the frame rate display
//  - overlay messages (these are up for a while then time out)
//
// We also update the frame number and FPS counter. Note that numerous
// frames may be produced from the same scene, so you can expect the frame
// numbers to be higher than the number of scenes sent by the simulator.
static void redrawDisplay() {
    // If a movie is being generated, save frames.
    // ------------------------------------------------------------
    if (saveNextFrameToMovie) {
        saveNextFrameToMovie = false;
        stringstream filename;
        filename << movieDir;
        filename << "/Frame";
        filename << setw(4) << setfill('0') << movieFrame++;
        filename << ".png";
        writeImage(filename.str());
    }

    // Render the scene and extract the screen text.
    // ------------------------------------------------------------
    std::vector<string> screenText;
    renderScene(&screenText);

    // Draw menus.
    // ------------------------------------------------------------
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, viewWidth, 0, viewHeight);
    glScalef(1, -1, 1);
    glTranslatef(0, GLfloat(-viewHeight), 0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    int menux = 10;
    for (int i = 0; i < (int) menus.size(); i++)
        menux += menus[i].draw(menux, viewHeight-10);

    // Draw sliders.
    // ------------------------------------------------------------
    int slidery = viewHeight-35;
    for (int i = 0; i < (int) sliders.size(); i++)
        slidery = sliders[i].draw(slidery);

    // Draw the "heads-up" display.
    // ------------------------------------------------------------
    glColor3f(1.0f, 0.5f, 0.0f);
    void* font = (void*)GLUT_BITMAP_HELVETICA_18;
    GLfloat nextLine = 25, lineHeight = 18;

    // Frames per second
    if (showFPS) {
    char fpstxt[64]; sprintf(fpstxt, "FPS:   %.1f", fps); // 1 decimal place
        glRasterPos2f(10, nextLine);
        for (const char* p = fpstxt; *p; ++p)
            glutBitmapCharacter(font, *p);
        nextLine += lineHeight;
    }

    // Simulation time
    if (showSimTime) {
        char timetxt[64]; sprintf(timetxt, "Time:  %.3f", lastSceneSimTime);
        glRasterPos2f(10, nextLine);
        for (const char* p = timetxt; *p; ++p)
            glutBitmapCharacter(font, *p);
        nextLine += lineHeight;
    }

    // Frame number
    if (showFrameNum) {
        char cnttxt[64]; sprintf(cnttxt, "Frame: %d", frameCounter);
        glRasterPos2f(10, nextLine);
        for (const char* p = cnttxt; *p; ++p)
            glutBitmapCharacter(font, *p);
        nextLine += lineHeight;
    }

    // User specified screen text
    for (int i=0; i<(int)screenText.size(); ++i) {
        glRasterPos2f(10, nextLine);
        for (const char* p = screenText[i].c_str(); *p; ++p)
            glutBitmapCharacter(font, *p);
        nextLine += lineHeight;
    }

    // Draw a message overlay 
    // (center box, with text left justified in box).
    // ------------------------------------------------------------
    if (displayOverlayMessage) {
        void* font = (void*)GLUT_BITMAP_HELVETICA_18;
        const int lineSpacing = 25;
        int width = 0;
        for (unsigned i=0; i < overlayMessageLines.size(); ++i)
            width = max(width, glutBitmapLength(font,
                                    (unsigned char*)overlayMessageLines[i].c_str()));
        int height = lineSpacing * (int)overlayMessageLines.size();

        GLfloat xpos = GLfloat(max(0, (viewWidth-width)/2));
        GLfloat ypos = GLfloat(max(0, (viewHeight-height)/2));

        GLfloat xborder = 4;
        Vec<2,GLfloat> tl(xpos-xborder, ypos-19);
        Vec<2,GLfloat> br(xpos+(float)width+xborder, ypos-19+(float)height);

        // Draw an orange background.
        glColor3f(1.0f, 0.5f, 0.f);
        glBegin(GL_QUADS);
        glVertex2f(tl[0],tl[1]);
        glVertex2f(tl[0],br[1]);
        glVertex2f(br[0],br[1]);
        glVertex2f(br[0],tl[1]);
        glEnd();

        // Draw dark gray text.
        glColor3f(0.2f, 0.2f, 0.2f);
        for (unsigned line=0; line < overlayMessageLines.size(); ++line) {
            glRasterPos2f(xpos, ypos);
            for (unsigned i = 0; i < overlayMessageLines[line].size(); i++)
                glutBitmapCharacter(font, overlayMessageLines[line][i]);
            ypos += (float)lineSpacing;
        }
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glutSwapBuffers();
    ++fpsCounter;
    ++frameCounter;

    lastRedisplayDone = realTime();
}

// These are set when a mouse button is clicked and then referenced
// while the mouse is dragged with that button down.
static int clickModifiers;
static int clickButton;
static int clickX;
static int clickY;
static int clickedSlider = -1;
static bool rotateCenterValid = false;
static fVec3 rotateCenter;
static float sceneScale; // for scaling translations

// Handle the initial press of a mouse button. Standard glut does not
// deal with the mouse wheel, but freeglut and some tweaked gluts treat
// it as a "button press" with button codes 3 (up) and 4 (down). This
// function handles those properly if they are present.
static void mouseButtonPressedOrReleased(int button, int state, int x, int y) {
    const int GlutWheelUp = 3, GlutWheelDown = 4;

    if (scene == NULL)
        return; // probably mouse click during startup

    float radius;
    fVec3 sceneCenter;
    std::unique_lock<std::mutex> lock(sceneMutex); //------- LOCK SCENE -------
    computeSceneBounds(scene, radius, sceneCenter);
    lock.unlock();                                 //----- UNLOCK SCENE -------
    sceneScale = std::max(radius, 0.1f);

    // "state" (pressed/released) is irrelevant for mouse wheel. However, if
    // we're being called by freeglut we'll get called twice, while (patched)
    // glut calls just once, with state=GLUT_UP.
    if ((button == GlutWheelUp || button == GlutWheelDown)
        && (state==GLUT_UP))
    {
        // Scroll wheel.
        const float ZoomFractionPerWheelClick = 0.15f;   // 15% scene radius

        const int   direction = button==GlutWheelUp ? -1 : 1;
        const float zoomBy    = direction * (ZoomFractionPerWheelClick * sceneScale);

        std::unique_lock<std::mutex> lock(sceneMutex); //----- LOCK SCENE -----
        X_GC.updP() += X_GC.R()*fVec3(0, 0, zoomBy);
        requestPassiveRedisplay();               //------ PASSIVE REDISPLAY ---
        lock.unlock();                           //------ UNLOCK SCENE --------
        return;
    }

    // Not mouse wheel; currently ignore "button up" message except to
    // forget the rotation center. This is needed to get around a glut
    // bug described in mouseDragged().
    if (state == GLUT_UP) {
        if (clickButton == GLUT_LEFT_BUTTON)
            rotateCenterValid = false;
        return;
    }

    // Handle state == GLUT_DOWN:

    // Remember which button was pressed; we'll deal with it in mouseDragged().
    clickModifiers = glutGetModifiers();
    clickButton = button;
    clickX = x;
    clickY = y;

    // See if this click was on a slider.
    clickedSlider = -1;
    for (int i = 0; i < (int) sliders.size(); i++) {
        if (sliders[i].mousePressed(x, y)) {
            clickedSlider = i;
            return;
        }
    }

    // Left button is rotation; when it is first pressed we calcuate the center
    // of rotation; we'll do the actual rotating in mouseDragged().
    if (clickButton == GLUT_LEFT_BUTTON) {
        float distToCenter = (sceneCenter-X_GC.p()).norm();
        float defaultDistance = sceneScale;
        if (distToCenter > defaultDistance)
            rotateCenter = sceneCenter;
        else {
            fVec3 cameraDir = X_GC.R()*fVec3(0, 0, -1);
            fVec3 lookAt = X_GC.p()+defaultDistance*cameraDir;
            float fract = (defaultDistance-distToCenter)/defaultDistance;
            rotateCenter = fract*lookAt+(1-fract)*sceneCenter;
        }
        rotateCenterValid = true;
    }
}

// This function is called when the mouse is moved while a button is being held
// down. When the button was first clicked, we recorded which one it was in
// clickButton, and where the mouse was then in (clickX,clickY). We update
// (clickX,clickY) each call here to reflect where it was last seen.
//
// (sherm 20101121: on Windows glut I get a spurious "mouse dragged" call
// when I have a menu displayed and then click outside the menu but inside
// the graphics window.)
static void mouseDragged(int x, int y) {
    if (clickedSlider > -1) {
        sliders[clickedSlider].mouseDragged(x);
        return;
    }

    // 1/4 degree per pixel
    const float AnglePerPixel = 0.25*((float)SimTK_PI/180);

    // map 1 pixel move to 1% of scale
    const float TranslateFracPerPixel = 0.01f;
    const float translatePerPixel = TranslateFracPerPixel * sceneScale;
    const int dx = clickX-x, dy = clickY-y;

    // translate: right button or shift-left button
    if (  clickButton == GLUT_RIGHT_BUTTON
      || (clickButton == GLUT_LEFT_BUTTON && clickModifiers & GLUT_ACTIVE_SHIFT))
    {
        std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE --------
        X_GC.updP() += translatePerPixel*X_GC.R()*fVec3(float(dx),float(-dy),0);
                                                      //---- UNLOCK SCENE ------
    }

    // zoom: middle button or alt-left button (or mouse wheel; see above)
    else if (  clickButton == GLUT_MIDDLE_BUTTON
           || (clickButton == GLUT_LEFT_BUTTON && clickModifiers & GLUT_ACTIVE_ALT))
    {
        std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE --------
        X_GC.updP() += translatePerPixel* X_GC.R()*fVec3(0,0,float(dy));
                                                      //---- UNLOCK SCENE ------
    }

    // rotate: left button alone: rotate scene left/right or up/down
    //         ctrl-left button:  roll scene about camera direction
    else if (clickButton == GLUT_LEFT_BUTTON && rotateCenterValid) {
        fVec3 cameraPos = X_GC.p();
        fVec3 cameraDir = X_GC.R()*fVec3(0, 0, -1);
        fVec3 upDir = X_GC.R()*fVec3(0, 1, 0);
        fRotation r;
        if (clickModifiers & GLUT_ACTIVE_CTRL)
            r.setRotationFromAngleAboutAxis(AnglePerPixel*(dy-dx), ZAxis);
        else
            r.setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence,
                AnglePerPixel*dy, XAxis, AnglePerPixel*dx, YAxis);
        r = X_GC.R()*r*~X_GC.R();
        cameraPos = r*(cameraPos-rotateCenter)+rotateCenter;
        cameraDir = r*cameraDir;
        upDir = r*upDir;

        std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE --------
        X_GC.updP() = cameraPos;
        X_GC.updR().setRotationFromTwoAxes(fUnitVec3(-cameraDir), ZAxis, upDir, YAxis);
                                                      //---- UNLOCK SCENE ------
    }
    else
        return;

    // Something changed.

    clickX = x; clickY = y;
    requestPassiveRedisplay();                   //------ PASSIVE REDISPLAY ---
}

// Currently the only "passive" mouse motion we care about is if the
// mouse is in an active menu, in which case its position marks the
// menu item that would be selected upon clicking.
// TODO: pass this to the simulator otherwise.
static void mouseMoved(int x, int y) {
    for (int i = 0; i < (int) menus.size(); i++)
        menus[i].mouseMoved(x, y);
}

// Glut distinguishes ordinary ASCII key presses from special
// ones like arrows and function keys. We will map either case
// into the same "key pressed" command to send to the simulator.
static void ordinaryKeyPressed(unsigned char key, int x, int y) {
    if (writeToSimulator) {
        char command = KeyPressed;
        WRITE(outPipe, &command, 1);
        unsigned char buffer[2];
        buffer[0] = key;
        buffer[1] = 0;
        unsigned char modifiers = glutGetModifiers();
        if ((modifiers & GLUT_ACTIVE_SHIFT) != 0)
            buffer[1] += Visualizer::InputListener::ShiftIsDown;
        if ((modifiers & GLUT_ACTIVE_CTRL) != 0)
            buffer[1] += Visualizer::InputListener::ControlIsDown;
        if ((modifiers & GLUT_ACTIVE_ALT) != 0)
            buffer[1] += Visualizer::InputListener::AltIsDown;
        WRITE(outPipe, buffer, 2);
    }
}

static void specialKeyPressed(int key, int x, int y) {
    if (writeToSimulator) {
        char command = KeyPressed;
        WRITE(outPipe, &command, 1);
        unsigned char buffer[2];
        buffer[0] = (unsigned char)key; // this is the special key code
        buffer[1] = Visualizer::InputListener::IsSpecialKey;
        unsigned char modifiers = glutGetModifiers();
        if ((modifiers & GLUT_ACTIVE_SHIFT) != 0)
            buffer[1] &= Visualizer::InputListener::ShiftIsDown;
        if ((modifiers & GLUT_ACTIVE_CTRL) != 0)
            buffer[1] &= Visualizer::InputListener::ControlIsDown;
        if ((modifiers & GLUT_ACTIVE_ALT) != 0)
            buffer[1] &= Visualizer::InputListener::AltIsDown;
        WRITE(outPipe, buffer, 2);
    }
}

static void shutdownWhenKeyPressed(unsigned char, int, int) {
    exit(0);
}

// This function is called when the timer goes off after an overlay message
// has been displayed long enough. On some platforms the timer may go off
// early -- if so we'll reissue it here for the remaining time.
static void disableOverlayTimer(int value) {
    const double now = realTime();
    if (now + 0.1 < overlayEndTime) { // 100ms slop
        // went off too early
        glutTimerFunc((unsigned)((overlayEndTime-now)*1000), // ms
            disableOverlayTimer, 0);
        return;
    }
    displayOverlayMessage = false;
    overlayStartTime = overlayEndTime = NaN;
    requestPassiveRedisplay();                   //------ PASSIVE REDISPLAY ---
}

static void setOverlayMessage(const string& message) {
    if (displayOverlayMessage) return; // there is already one up
    overlayMessageLines.clear();
    string::size_type start = 0;
    while (start < message.size()) {
        string::size_type newline = message.find_first_of('\n', start);
        if (newline == string::npos) newline = message.size();
        overlayMessageLines.push_back(message.substr(start, newline-start));
        start = newline+1;
    }

    if (overlayMessageLines.empty())
        return;

    displayOverlayMessage = true;
    requestPassiveRedisplay();                   //------ PASSIVE REDISPLAY ---

    overlayStartTime = realTime();
    overlayEndTime = overlayStartTime + OverlayDisplayTimeInSec;
    glutTimerFunc((unsigned)(OverlayDisplayTimeInSec*1000), // ms
        disableOverlayTimer, 0);
}

class SaveImageTask : public ParallelWorkQueue::Task {
public:
    Array_<unsigned char> data;
    SaveImageTask(const string& filename, int width, int height) : filename(filename), width(width), height(height), data(width*height*3) {
    }
    void execute() override {
        // Flip the image vertically, since OpenGL and PNG use different row orders.

        const int rowLength = 3*width;
        for (int row = 0; row < height/2; ++row) {
            const int base1 = row*rowLength;
            const int base2 = (height-1-row)*rowLength;
            for (int i = 0; i < rowLength; i++) {
                unsigned char temp = data[base1+i];
                data[base1+i] = data[base2+i];
                data[base2+i] = temp;
            }
        }
        LodePNG::encode(filename, data.empty() ? 0 : &data[0], width, height, 2, 8);
    }
private:
    string filename;
    int width, height;
};

static void writeImage(const string& filename) {
    int width = ((viewWidth+3)/4)*4; // must be a multiple of 4 pixels
    int height = viewHeight;

    // Create offscreen buffers for rendering the image.

    GLuint frameBuffer, colorBuffer, depthBuffer;
    glGenFramebuffersEXT(1, &frameBuffer);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBuffer);
    glGenRenderbuffersEXT(1, &colorBuffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, colorBuffer);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGB8, width, height);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, colorBuffer);
    glGenRenderbuffersEXT(1, &depthBuffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depthBuffer);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, width, height);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depthBuffer);

    // Render the image and load it into memory.

    renderScene();
    SaveImageTask* task = new SaveImageTask(filename, width, height);
    Array_<unsigned char>& data = task->data;
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &data[0]);
    glDeleteRenderbuffersEXT(1, &colorBuffer);
    glDeleteRenderbuffersEXT(1, &depthBuffer);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glDeleteFramebuffersEXT(1, &frameBuffer);

    // Add it to the queue to be saved to disk.

    imageSaverQueue.addTask(task);
}

static void saveImage() {
    struct stat statInfo;
    int counter = 0;
    string filename;
    do {
        counter++;
        stringstream namestream;
        namestream << simulatorExecutableName.c_str() << "_";
        namestream << counter;
        namestream << ".png";
        filename = namestream.str();
    } while (stat(filename.c_str(), &statInfo) == 0);
    writeImage(filename);
    setOverlayMessage("Image saved as:\n"+filename);
}

static void saveMovie() {
    struct stat statInfo;
    int counter = 0;
    string dirname;
    do {
        counter++;
        stringstream namestream;
        namestream << simulatorExecutableName.c_str() << "_";
        namestream << counter;
        dirname = namestream.str();
    } while (stat(dirname.c_str(), &statInfo) == 0);
#ifdef _WIN32
    int result = mkdir(dirname.c_str());
#else
    int result = mkdir(dirname.c_str(), 0777);
#endif
    if (result == -1)
        setOverlayMessage("Failed to create directory:\n"+dirname);
    else {
        movieDir = dirname;
        movieFrame = 1;
        savingMovie = true;
        setOverlayMessage("Capturing frames in:\n"+dirname);
    }
}

class ReadingInterrupted : public std::exception {};

// Read a particular number of bytes from srcPipe to the given buffer.
// This will hang until the expected number of bytes has been received.
// Throws ReadingInterrupted if the srcPipe is closed.
static void readDataFromPipe(int srcPipe, unsigned char* buffer, int bytes) {
    int totalRead = 0;
    while (totalRead < bytes) {
        auto retval = READ(srcPipe, buffer + totalRead, bytes - totalRead);
        SimTK_ERRCHK4_ALWAYS(retval!=-1, "simbody-visualizer",
            "An attempt to read() %d bytes from pipe %d failed with errno=%d (%s).", 
            bytes - totalRead, srcPipe, errno, strerror(errno));
        // The pipe was closed, perhaps because the simulator was closed.
        // Without this check, we end up in an infinite loop.
        if (retval == 0) throw ReadingInterrupted();

        totalRead += retval;
    }
}
// Throws ReadingInterrupted if inPipe is closed.
static void readData(unsigned char* buffer, int bytes) {
    readDataFromPipe(inPipe, buffer, bytes);
}

// We have just processed a StartOfScene command. Read in all the scene
// elements until we see an EndOfScene command. We allocate a new Scene
// object to hold the scene and return a pointer to it. Don't forget to
// delete that object when you are done with it.
static Scene* readNewScene() {
    unsigned char buffer[256];
    float*          floatBuffer = (float*)          buffer;
    int*            intBuffer   = (int*)            buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;

    Scene* newScene = new Scene;

    // Simulated time for this frame comes first.
    readData(buffer, sizeof(float));
    newScene->simTime = floatBuffer[0];

    bool finished = false;
    while (!finished) {
        readData(buffer, 1);
        char command = buffer[0];

        switch (command) {

        case Shutdown:
            shutdown(); // doesn't return
            break;

        case EndOfScene:
            finished = true;
            break;

        // Add a scene element that uses an already-cached mesh.
        case AddPointMesh:
        case AddWireframeMesh:
        case AddSolidMesh: {
            readData(buffer, 13*sizeof(float)+2*sizeof(short));
            fTransform position;
            position.updR().setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
            position.updP() = fVec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 scale = fVec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            fVec4 color = fVec4(floatBuffer[9], floatBuffer[10], floatBuffer[11], floatBuffer[12]);
            short representation = (command == AddPointMesh ? DecorativeGeometry::DrawPoints : (command == AddWireframeMesh ? DecorativeGeometry::DrawWireframe : DecorativeGeometry::DrawSurface));
            unsigned short meshIndex = shortBuffer[13*sizeof(float)/sizeof(short)];
            unsigned short resolution = shortBuffer[13*sizeof(float)/sizeof(short)+1];
            RenderedMesh mesh(position, scale, color, representation, meshIndex, resolution);
            if (command != AddSolidMesh)
                newScene->drawnMeshes.push_back(mesh);
            else if (color[3] == 1)
                newScene->solidMeshes.push_back(mesh);
            else
                newScene->transparentMeshes.push_back(mesh);
            if (meshIndex < NumPredefinedMeshes && (meshes[meshIndex].size() <= resolution || meshes[meshIndex][resolution] == NULL)) {
                // A real mesh will be generated from this the next
                // time the scene is redrawn.
                std::lock_guard<std::mutex> lock(sceneMutex); //-- LOCK SCENE --
                pendingCommands.insert(pendingCommands.begin(), new PendingStandardMesh(meshIndex, resolution));
                                                              //- UNLOCK SCENE -
            }
            break;
        }

        case AddLine: {
            readData(buffer, 10*sizeof(float));
            fVec3 color = fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            float thickness = floatBuffer[3];
            int index;
            int numLines = (int)newScene->lines.size();
            for (index = 0; index < numLines && (color != newScene->lines[index].getColor() || thickness != newScene->lines[index].getThickness()); index++)
                ;
            if (index == numLines)
                newScene->lines.push_back(RenderedLine(color, thickness));
            vector<GLfloat>& line = newScene->lines[index].getLines();
            line.push_back(floatBuffer[4]);
            line.push_back(floatBuffer[5]);
            line.push_back(floatBuffer[6]);
            line.push_back(floatBuffer[7]);
            line.push_back(floatBuffer[8]);
            line.push_back(floatBuffer[9]);
            break;
        }

        case AddText: {
            readData(buffer, 12*sizeof(float)+3*sizeof(short));
            fTransform X_GT;
            X_GT.updR().setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
            X_GT.updP() = fVec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 scale = fVec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            fVec3 color = fVec3(floatBuffer[9], floatBuffer[10], floatBuffer[11]);
            unsigned short* shortp = &shortBuffer[12*sizeof(float)/sizeof(short)];
            bool faceCamera = (shortp[0] != 0);
            bool isScreenText = (shortp[1] != 0);
            short length = shortp[2];
            readData(buffer, length);

            if (isScreenText)
                newScene->screenText.push_back(
                    ScreenText(string((char*)buffer, length)));
            else
                newScene->sceneText.push_back(
                    RenderedText(X_GT, scale, color, string((char*)buffer, length), faceCamera));
            break;
        }

        case AddCoords: {
            readData(buffer, 12*sizeof(float));
            fRotation rotation;
            rotation.setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], 
                                                     floatBuffer[1], 
                                                     floatBuffer[2]));
            fVec3 position(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 axisLengths(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            fVec3 textScale = fVec3(0.2f*min(axisLengths));
            float lineThickness = 1;
            fVec3 color = fVec3(floatBuffer[9], floatBuffer[10], floatBuffer[11]);
            int index;
            int numLines = (int)newScene->lines.size();
            for (index = 0; 
                 index < numLines && (color != newScene->lines[index].getColor() 
                                      || newScene->lines[index].getThickness() 
                                                              != lineThickness); 
                 index++)
                ;
            if (index == numLines)
                newScene->lines.push_back(RenderedLine(color, lineThickness));
            vector<GLfloat>& line = newScene->lines[index].getLines();
            fVec3 end = position+rotation*fVec3(axisLengths[0], 0, 0);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->sceneText.push_back(RenderedText(end, textScale, color, "X"));
            end = position+rotation*fVec3(0, axisLengths[1], 0);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->sceneText.push_back(RenderedText(end, textScale, color, "Y"));
            end = position+rotation*fVec3(0, 0, axisLengths[2]);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->sceneText.push_back(RenderedText(end, textScale, color, "Z"));
            break;
        }

        // Define a new mesh that will be assigned the next available mesh
        // index. It will be cached here and then can be referenced in this
        // scene and others by using it mesh index.
        case DefineMesh: {
            readData(buffer, 2*sizeof(short));
            PendingMesh* mesh = new PendingMesh(); // assigns next mesh index
            int numVertices = shortBuffer[0];
            int numFaces = shortBuffer[1];
            mesh->vertices.resize(3*numVertices, 0);
            mesh->normals.resize(3*numVertices);
            mesh->faces.resize(3*numFaces);
            readData((unsigned char*)&mesh->vertices[0], (int)(mesh->vertices.size()*sizeof(float)));
            readData((unsigned char*)&mesh->faces[0], (int)(mesh->faces.size()*sizeof(short)));

            // Compute normal vectors for the mesh.

            vector<fVec3> normals(numVertices, fVec3(0));
            for (int i = 0; i < numFaces; i++) {
                int v1 = mesh->faces[3*i];
                int v2 = mesh->faces[3*i+1];
                int v3 = mesh->faces[3*i+2];
                fVec3 vert1(mesh->vertices[3*v1], mesh->vertices[3*v1+1], mesh->vertices[3*v1+2]);
                fVec3 vert2(mesh->vertices[3*v2], mesh->vertices[3*v2+1], mesh->vertices[3*v2+2]);
                fVec3 vert3(mesh->vertices[3*v3], mesh->vertices[3*v3+1], mesh->vertices[3*v3+2]);
                fVec3 norm = (vert2-vert1)%(vert3-vert1);
                float length = norm.norm();
                if (length > 0) {
                    norm /= length;
                    normals[v1] += norm;
                    normals[v2] += norm;
                    normals[v3] += norm;
                }
            }
            for (int i = 0; i < numVertices; i++) {
                normals[i] = normals[i].normalize();
                mesh->normals[3*i] = normals[i][0];
                mesh->normals[3*i+1] = normals[i][1];
                mesh->normals[3*i+2] = normals[i][2];
            }

            // A real mesh will be generated from this the next
            // time the scene is redrawn.
            std::lock_guard<std::mutex> lock(sceneMutex); //--- LOCK SCENE ----
            pendingCommands.insert(pendingCommands.begin(), mesh);
            break;                                        //--- UNLOCK SCENE --
        }

        default:
            SimTK_ASSERT_ALWAYS(false, "Unexpected scene data sent to visualizer");
        }
    }

    return newScene;
}

// This is the main program for the listener thread. It reads continuously
// from the input pipe, which contains data from the simulator's calls
// to a Visualizer object. Any changes to the scene must wait until the
// rendering thread has finished with the current frame. The listener
// thread must not make any gl or glut calls; any operations that require
// those are queued for later handing in the rendering thread.
void listenForInput() {
    unsigned char buffer[256];
    float* floatBuffer = (float*) buffer;
    int* intBuffer = (int*) buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;

    try
  { while (true) {
        bool issuedActiveRedisplay = false;
        // Read commands from the simulator.
        readData(buffer, 1);
        switch (buffer[0]) {
        case DefineMenu: {
            readData(buffer, sizeof(short));
            int titleLength = shortBuffer[0];
            vector<char> titleBuffer(titleLength);
            readData((unsigned char*)&titleBuffer[0], titleLength);
            string title(&titleBuffer[0], titleLength);
            readData(buffer, sizeof(int));
            const int menuId = intBuffer[0];
            readData(buffer, sizeof(short));
            int numItems = shortBuffer[0];
            vector<pair<string, int> > items(numItems);
            for (int index = 0; index < numItems; index++) {
                readData(buffer, 2*sizeof(int));
                items[index].second = intBuffer[0];
                vector<char> textBuffer(intBuffer[1]);
                readData((unsigned char*)&textBuffer[0], intBuffer[1]);
                items[index].first = string(&textBuffer[0], intBuffer[1]);
            }
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            menus.push_back(Menu(title, menuId, items, menuSelected));
            break;                                        //--- UNLOCK SCENE ---
        }
        case DefineSlider: {
            readData(buffer, sizeof(short));
            int titleLength = shortBuffer[0];
            vector<char> titleBuffer(titleLength);
            readData((unsigned char*)&titleBuffer[0], titleLength);
            string title(&titleBuffer[0], titleLength);
            readData(buffer, sizeof(int)+3*sizeof(float));
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            sliders.push_back(Slider(title, intBuffer[0], floatBuffer[1], floatBuffer[2], floatBuffer[3]));
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetSliderValue: {
            readData(buffer, sizeof(int));
            const int sliderId = intBuffer[0];
            readData(buffer, sizeof(float));
            const float newValue = floatBuffer[0];
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            Slider* sp = findSliderById(sliderId);
            if (sp) sp->changeValue(newValue);
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetSliderRange: {
            readData(buffer, sizeof(int));
            const int sliderId = intBuffer[0];
            readData(buffer, 2*sizeof(float));
            const float newMin = floatBuffer[0];
            const float newMax = floatBuffer[1];
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            Slider* sp = findSliderById(sliderId);
            if (sp) sp->changeRange(newMin, newMax);
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetCamera: {
            readData(buffer, 6*sizeof(float));
            fVec3 R(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            fVec3 p(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            pendingCommands.push_back(new PendingSetCameraTransform(R, p));
            break;                                        //--- UNLOCK SCENE ---
        }
        case ZoomCamera: {
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            pendingCommands.push_back(new PendingCameraZoom());
            break;                                        //--- UNLOCK SCENE ---
        }
        case LookAt: {
            readData(buffer, 6*sizeof(float));
            fVec3 point(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            fVec3 updir(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            fVec3 pt2camera = X_GC.p()-point;
            if (pt2camera.normSqr() < square(1e-6))
                pt2camera = fVec3(X_GC.z()); // leave unchanged
            X_GC.updR().setRotationFromTwoAxes(fUnitVec3(pt2camera), ZAxis, 
                                               updir, YAxis);
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetFieldOfView: {
            readData(buffer, sizeof(float));
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            fieldOfView = floatBuffer[0];
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetClipPlanes: {
            readData(buffer, 2*sizeof(float));
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            nearClip = floatBuffer[0];
            farClip = floatBuffer[1];
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetSystemUpDirection: {
            readData(buffer, 2);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            groundNormal = CoordinateDirection( CoordinateAxis((int)buffer[0]),
                                                (int)(signed char)buffer[1] );
            X_GC.updR().setRotationFromTwoAxes
               (groundNormal, YAxis, X_GC.z(), ZAxis); // attempt to keep z
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetGroundHeight: {
            readData(buffer, sizeof(float));
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            groundHeight = floatBuffer[0];
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetWindowTitle: {
            readData(buffer, sizeof(short));
            int titleLength = shortBuffer[0];
            vector<char> titleBuffer(titleLength);
            readData((unsigned char*)&titleBuffer[0], titleLength);
            string title(&titleBuffer[0], titleLength);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            pendingCommands.push_back(new PendingWindowTitleChange(title));
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetMaxFrameRate: {
            readData(buffer, sizeof(float));
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            maxFrameRate = floatBuffer[0];
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetBackgroundColor: {
            readData(buffer, 3*sizeof(float));
            fVec3 color(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            pendingCommands.push_back(new PendingBackgroundColorChange(color));
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetShowShadows: {
            readData(buffer, sizeof(short));
            const bool shouldShow = (shortBuffer[0] != 0);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            showShadows = shouldShow;
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetBackgroundType: {
            readData(buffer, sizeof(short));
            const Visualizer::BackgroundType type =
                Visualizer::BackgroundType(shortBuffer[0]);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            showGround = (type == Visualizer::GroundAndSky);
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetShowFrameRate: {
            readData(buffer, sizeof(short));
            const bool shouldShow = (shortBuffer[0] != 0);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            showFPS = shouldShow;
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetShowSimTime: {
            readData(buffer, sizeof(short));
            const bool shouldShow = (shortBuffer[0] != 0);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            showSimTime = shouldShow;
            break;                                        //--- UNLOCK SCENE ---
        }
        case SetShowFrameNumber: {
            readData(buffer, sizeof(short));
            const bool shouldShow = (shortBuffer[0] != 0);
            std::lock_guard<std::mutex> lock(sceneMutex); //---- LOCK SCENE ----
            showFrameNum = shouldShow;
            break;                                        //--- UNLOCK SCENE ---
        }
        case StartOfScene: {
            Scene* newScene = readNewScene();
            std::unique_lock<std::mutex> lock(sceneMutex); //--- LOCK SCENE ----
            if (scene != NULL) {
                // -------- WAIT FOR CONDITION --------
                // Give up the lock and wait for notice.
                // Pass a predicate function so that we continue waiting if we
                // caught a spurious wakeup.
                sceneHasBeenDrawn.wait(lock,
                        [&] { return scene->sceneHasBeenDrawn; });
                // Previous scene has been drawn.
                delete scene; scene = 0;
            }
            // Swap in the new scene.
            scene = newScene;
            saveNextFrameToMovie = savingMovie;
            lock.unlock();                                 //-- UNLOCK SCENE ---
            forceActiveRedisplay();               //------- ACTIVE REDISPLAY ---
            issuedActiveRedisplay = true;
            break;
        }

        case Shutdown:
            shutdown(); // doesn't return
            break;

        case StopCommunication: {
            // Before we stop listening to the simulator, we tell the simulator
            // to stop listening to the GUI.
            writeToSimulator = false;
            int retval = CLOSE(outPipe);
            if (retval == -1) {
                std::cout << "Warning in simbody-visualizer: "
                    << "An attempt to close() pipe " << outPipe
                    << " failed with errno=" << errno
                    << " (" << strerror(errno) << ")." << std::endl;
            }

            // Stop this listener thread, but do not kill the entire GUI; the
            // user may still want to save images, etc.

            // If the user presses a key after we have disconnected from the
            // simulator, then shut down.
            glutKeyboardFunc(shutdownWhenKeyPressed);
            // Remove the callback for special keys.
            glutSpecialFunc(nullptr);

            return;
        }

        default:
            SimTK_ERRCHK1_ALWAYS(!"unrecognized command", "listenForInput()",
                "simbody-visualizer received unexpected command %u "
                "from simulator. Can't continue.",
                (unsigned)buffer[0]);
        }

        // Do this after every received command.
        if (!issuedActiveRedisplay)
            requestPassiveRedisplay();         //------- PASSIVE REDISPLAY --
    }
  } catch (const ReadingInterrupted&) {
        // Stop listening, because the simulator was closed.
  } catch (const std::exception& e) {
        std::cout << "simbody-visualizer listenerThread: unrecoverable error:\n";
        std::cout << e.what() << std::endl;
    }
}

static void dumpAboutMessageToConsole() {
    printf("\n\n=================== ABOUT SIMBODY VISUALIZER ===================\n");
    printf("Simbody(tm) %s visualizer (protocol rev. %u)\n",
        simbodyVersionStr.c_str(), ProtocolVersion);
    printf("\nName of invoking executable: %s\n", simulatorExecutableName.c_str());
    printf(  "Current working directory:\n  %s\n",
        Pathname::getCurrentWorkingDirectory().c_str());
    printf(  "simbody-visualizer executable:\n  %s\n",
        Pathname::getThisExecutablePath().c_str());
    printf(  "Current window size: %d X %d\n", viewWidth, viewHeight);
    printf("\nGL version:   %s\n", glGetString(GL_VERSION));
    printf(  "GLSL version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
    printf(  "GL renderer:  %s\n", glGetString(GL_RENDERER));
    printf(  "GL vendor:    %s\n", glGetString(GL_VENDOR));
    printf("\nVisualizer authors: Peter Eastman, Michael Sherman\n");
    printf(  "Support: Simbios, Stanford Bioengineering, NIH U54 GM072970\n");
    printf(  "https://simtk.org/home/simbody\n");
    printf("================================================================\n\n");
}

static const int MENU_VIEW_FRONT = 0;
static const int MENU_VIEW_BACK = 1;
static const int MENU_VIEW_LEFT = 2;
static const int MENU_VIEW_RIGHT = 3;
static const int MENU_VIEW_TOP = 4;
static const int MENU_VIEW_BOTTOM = 5;

static const int MENU_BACKGROUND_BLACK = 6;
static const int MENU_BACKGROUND_WHITE = 7;
static const int MENU_BACKGROUND_SKY = 8;

static const int MENU_SHOW_SHADOWS = 9;
static const int MENU_SHOW_FPS = 10;
static const int MENU_SHOW_SIM_TIME = 11;
static const int MENU_SHOW_FRAME_NUM = 12;

static const int MENU_SAVE_IMAGE = 13;
static const int MENU_SAVE_MOVIE = 14;

static const int MENU_ABOUT = 15;

// This is the handler for our built-in "View" pull down menu.
void viewMenuSelected(int option) {
    fRotation groundRotation;
    if (groundNormal.getAxis() == XAxis)
        groundRotation.setRotationFromTwoAxes(fUnitVec3(0, 1, 0), ZAxis, fVec3(1, 0, 0), YAxis);
    else if (groundNormal.getAxis() == ZAxis)
        groundRotation.setRotationFromTwoAxes(fUnitVec3(1, 0, 0), ZAxis, fVec3(0, 0, 1), YAxis);

    // Flip the camera 180 degress around one of the other axes if in a negative direction.
    if (groundNormal.getDirection() == -1)
        groundRotation = fRotation((float)Pi, groundNormal.getAxis().getNextAxis()) * groundRotation;

    switch (option) {
    case MENU_VIEW_FRONT:
        X_GC.updR().setRotationToIdentityMatrix();
        X_GC.updR() = groundRotation*X_GC.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_BACK:
        X_GC.updR().setRotationFromAngleAboutY((float)SimTK_PI);
        X_GC.updR() = groundRotation*X_GC.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_LEFT:
        X_GC.updR().setRotationFromAngleAboutY(-(float)(SimTK_PI/2));
        X_GC.updR() = groundRotation*X_GC.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_RIGHT:
        X_GC.updR().setRotationFromAngleAboutY((float)(SimTK_PI/2));
        X_GC.updR() = groundRotation*X_GC.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_TOP:
        X_GC.updR().setRotationFromAngleAboutX(-(float)(SimTK_PI/2));
        X_GC.updR() = groundRotation*X_GC.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_BOTTOM:
        X_GC.updR().setRotationFromAngleAboutX((float)(SimTK_PI/2));
        X_GC.updR() = groundRotation*X_GC.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_BACKGROUND_BLACK:
        showGround = false;
        backgroundColor = fVec3(0,0,0);
        setClearColorToBackgroundColor();
        break;
    case MENU_BACKGROUND_WHITE:
        showGround = false;
        backgroundColor = fVec3(1,1,1);
        setClearColorToBackgroundColor();
        break;
    case MENU_BACKGROUND_SKY:
        showGround = true;
        backgroundColor = fVec3(1,1,1);
        setClearColorToBackgroundColor();
        break;
    case MENU_SHOW_SHADOWS:
        showShadows = !showShadows;
        break;
    case MENU_SHOW_FPS:
        showFPS = !showFPS;
        break;
    case MENU_SHOW_SIM_TIME:
        showSimTime = !showSimTime;
        break;
    case MENU_SHOW_FRAME_NUM:
        showFrameNum = !showFrameNum;
        break;
    case MENU_SAVE_IMAGE:
        if (canSaveImages) {
            saveImage();
        } else 
            setOverlayMessage(
            "Sorry -- image capture not available due to your\n"
            "backlevel OpenGL. At least OpenGL 2.0 is required.\n"
            "See the About message for level information.");
        break;
    case MENU_SAVE_MOVIE:
        if (canSaveImages) {
            if (savingMovie) {
                savingMovie = false;
                setOverlayMessage("Frame capture off.");
            } else
                saveMovie();
        } else 
            setOverlayMessage(
            "Sorry -- movie capture not available due to your\n"
            "backlevel OpenGL. At least OpenGL 2.0 is required.\n"
            "See the About message for level information.");
        break;
    case MENU_ABOUT:
        dumpAboutMessageToConsole();
        setOverlayMessage("About: see console window");
        break;
    }

    requestPassiveRedisplay();                  //----- PASSIVE REDISPLAY ----
}

static const int DefaultWindowWidth  = 800;
static const int DefaultWindowHeight = 600;


// The glut callback for chaning window size.
static void changeSize(int width, int height) {
    if (height == 0)
        height = 1;
    viewWidth = width;
    viewHeight = height;
}


// This seems to be necessary on Mac and Linux where the built-in glut
// idle hangs if there is no activity in the gl window. In any case we
// use it to issue a redisplay when one is needed but no scene is forthcoming,
// and for final cleanup of the FPS display after scenes stop arriving.
static void keepAliveIdleFunc() {
    static const double SleepTime   = 1./100; // 10ms
    static const double MopUpTime   = 1.; // to clean up FPS display

    const double passiveRedisplayInterval = 1/(double)maxFrameRate;

    sleepInSec(SleepTime); // take a short break

    // If it has been at least one passive frame time, redisplay if there
    // has been a passive redisplay request.
    const double idle = realTime() - lastRedisplayDone;
    if (passiveRedisplayRequested && idle > passiveRedisplayInterval) {
        forcePassiveRedisplay();                //------ FORCE REDISPLAY -----
        return;
    }

    // If it has been so long since the last frame that we need to mop up
    // (currently that just means update the FPS display), then do that.
    // We want to calculate the last real FPS value, and then show a zero on
    // the subsequent update, and then stop updating.
    if (numMopUpDisplaysSinceLastRedisplay < 2 && idle >= MopUpTime) {
        if (numMopUpDisplaysSinceLastRedisplay)
            fpsCounter = 0; // clear the 1 frame we generated after 1st update
        forceMopUpRedisplay();                  //------ FORCE REDISPLAY -----
    }

}

static void setKeepAlive(bool enable) {
    if (enable) glutIdleFunc(keepAliveIdleFunc);
    else glutIdleFunc(0);
}

static void setVsync(bool enable) {
#ifdef _WIN32
    // I don't know how to disable vsync on Mac or Linux.
    if( wglSwapIntervalEXT )
      wglSwapIntervalEXT(enable?1:0); // 1==vsync; 0 is off
#endif
}


// This is executed from the main thread at startup.
static void shakeHandsWithSimulator(int fromSimPipe, int toSimPipe) {
    unsigned char handshakeCommand;
    readDataFromPipe(fromSimPipe, &handshakeCommand, 1);
    SimTK_ERRCHK2_ALWAYS(handshakeCommand == StartupHandshake,
        "simbody-visualizer::shakeHandsWithSimulator()",
        "Expected initial handshake command %u but received %u. Can't continue.",
        (unsigned)StartupHandshake, (unsigned)handshakeCommand);

    unsigned SimVersion;
    readDataFromPipe(fromSimPipe, (unsigned char*)&SimVersion, sizeof(unsigned int));
    SimTK_ERRCHK2_ALWAYS(SimVersion == ProtocolVersion,
        "simbody-visualizer::shakeHandsWithSimulator()",
        "The Simbody Visualizer class protocol version %u is not compatible with "
        " simbody-visualizer protocol %u; this may be an installation problem."
        " Can't continue.",
        SimVersion, ProtocolVersion);

    // Get Simbody version number as major,minor,patch
    readDataFromPipe(fromSimPipe, (unsigned char*)simbodyVersion, 3*sizeof(int));
    simbodyVersionStr = String(simbodyVersion[0]) + "." + String(simbodyVersion[1]);
    if (simbodyVersion[2]) simbodyVersionStr += "." + String(simbodyVersion[2]);

    unsigned exeNameLength;
    char exeNameBuf[256]; // just a file name, not a path name
    readDataFromPipe(fromSimPipe, (unsigned char*)&exeNameLength, sizeof(unsigned));
    SimTK_ASSERT_ALWAYS(exeNameLength <= 255,
        "simbody-visualizer: executable name length violates protocol.");
    readDataFromPipe(fromSimPipe, (unsigned char*)exeNameBuf, exeNameLength);
    exeNameBuf[exeNameLength] = (char)0;

    simulatorExecutableName = std::string(exeNameBuf, exeNameLength);

    WRITE(outPipe, &ReturnHandshake, 1);
    WRITE(outPipe, &ProtocolVersion, sizeof(unsigned));
}

// Received Shutdown message from simulator. Die immediately.
static void shutdown() {
    printf("\nsimbody-visualizer: received Shutdown message. Goodbye.\n");
    exit(0);
}

int main(int argc, char** argv) {
  try
  { bool talkingToSimulator = false;
      
    if (argc >= 3) {
        stringstream(argv[1]) >> inPipe;
        stringstream(argv[2]) >> outPipe;
        talkingToSimulator = true; // presumably those were the pipes
  } else {
        printf("\n**** VISUALIZER HAS NO SIMULATOR TO TALK TO ****\n");
        printf("The simbody-visualizer was invoked directly with no simulator\n");
        printf("process to talk to. Will attempt to bring up the display anyway\n");
        printf("in case you want to look at the About message.\n");
        printf("The simbody-visualizer is intended to be invoked programmatically.\n");
    }


    // Initialize GLUT, then perform initial handshake with the parent
    // from the main thread here.
    glutInit(&argc, argv);

    if (talkingToSimulator)
        shakeHandsWithSimulator(inPipe, outPipe);
    else {
        simbodyVersionStr = "?.?.?";
        simulatorExecutableName = "No simulator";
    }

    // Construct the default initial title.
    string title = "Simbody " + simbodyVersionStr + ": " + simulatorExecutableName;

    // Put the upper left corner of the glut window near the upper right
    // corner of the screen.
    int screenW = glutGet(GLUT_SCREEN_WIDTH);
    int screenH = glutGet(GLUT_SCREEN_HEIGHT);
    int windowPosX = (screenW - DefaultWindowWidth) - 50;   // 50 pixel margin
    int windowPosY = 50;

    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(windowPosX,windowPosY);
    glutInitWindowSize(DefaultWindowWidth, DefaultWindowHeight);
    glutCreateWindow(title.c_str());


    // Set up callback functions.
    glutDisplayFunc(redrawDisplay);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mouseButtonPressedOrReleased);
    glutMotionFunc(mouseDragged);
    glutPassiveMotionFunc(mouseMoved);
    glutKeyboardFunc(ordinaryKeyPressed);
    glutSpecialFunc(specialKeyPressed);

    //dumpAboutMessageToConsole();

    // On some systems (Windows at least), some of the gl functions may
    // need to be loaded dynamically.
    bool canFunction = initGlextFuncPointersIfNeeded(canSaveImages);
    if (!canFunction) {
        printf("\n\n**** FATAL ERROR ****\n");
        dumpAboutMessageToConsole();
        printf("\n**** FATAL ERROR ****\n");
        printf("Sorry, Simbody Visualizer can't function with this ancient version\n");
        printf("of OpenGL. See the message above for version information.\n");
        printf("Are you using an emulation through a remote desktop or\n");
        printf("virtual machine?\n");
        printf("**** FATAL ERROR **** Simbody Visualizer terminating.\n");
        return 1;
    }

    setVsync(true);

    // Set up lighting.

    GLfloat light_diffuse[]  = {0.8f, 0.8f, 0.8f, 1.0f};
    GLfloat light_position[] = {1.0f, 1.0f, 1.0f, 0.0f};
    GLfloat light_ambient[]  = {0.2f, 0.2f, 0.2f, 1.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambient);
    glClearColor(1, 1, 1, 1);
    glEnable(GL_LIGHT0);

    // Initialize miscellaneous OpenGL state.

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    // Make room for the predefined meshes.
    meshes.resize(NumPredefinedMeshes);
    // Note the first mesh index available for unique meshes
    // that are sent to the GUI for caching.
    nextMeshIndex = NumPredefinedMeshes;

    scene = NULL;
    pendingCommands.push_back(new PendingCameraZoom());

    vector<pair<string, int> > items;
    items.push_back(make_pair("View Direction/Front", MENU_VIEW_FRONT));
    items.push_back(make_pair("View Direction/Back", MENU_VIEW_BACK));
    items.push_back(make_pair("View Direction/Left", MENU_VIEW_LEFT));
    items.push_back(make_pair("View Direction/Right", MENU_VIEW_RIGHT));
    items.push_back(make_pair("View Direction/Top", MENU_VIEW_TOP));
    items.push_back(make_pair("View Direction/Bottom", MENU_VIEW_BOTTOM));
    items.push_back(make_pair("Background/Black", MENU_BACKGROUND_BLACK));
    items.push_back(make_pair("Background/White", MENU_BACKGROUND_WHITE));
    items.push_back(make_pair("Background/Ground and Sky", MENU_BACKGROUND_SKY));
    items.push_back(make_pair("Show//Hide/Shadows", MENU_SHOW_SHADOWS));
    items.push_back(make_pair("Show//Hide/Frame Rate", MENU_SHOW_FPS));
    items.push_back(make_pair("Show//Hide/Sim Time", MENU_SHOW_SIM_TIME));
    items.push_back(make_pair("Show//Hide/Frame #", MENU_SHOW_FRAME_NUM));
    items.push_back(make_pair("Save Image", MENU_SAVE_IMAGE));
    items.push_back(make_pair("Save Movie", MENU_SAVE_MOVIE));
    items.push_back(make_pair("About (to console)", MENU_ABOUT));
    menus.push_back(Menu("View", Visualizer::ViewMenuId, items, viewMenuSelected));

    // Spawn the listener thread. After this it runs independently.
    std::thread listenerThread;
    if (talkingToSimulator) {
        listenerThread = std::thread(listenForInput);
    } else {
        scene = new Scene;
    }

    // Avoid hangs on Mac & Linux; posts orphan redisplays on all platforms.
    setKeepAlive(true);

    // Enter the main loop. If there is nothing else to do we'll check if
    // someone was hoping for a redisplay and issue one.
    requestPassiveRedisplay();                   //------ PASSIVE REDISPLAY ---
    fpsBaseTime = realTime();
    glutMainLoop();
  } catch(const std::exception& e) {
      std::cout << "simbody-visualizer failed with exception:\n"
                << e.what() << std::endl;
      return 1;
    }

    return 0;
}


// Initialize function pointers for Windows GL extensions.
static bool initGlextFuncPointersIfNeeded(bool& glCanSaveImages) {
    glCanSaveImages = true;
#ifdef _WIN32
    wglSwapIntervalEXT = (PFNWGLSWAPINTERVALFARPROC)wglGetProcAddress( "wglSwapIntervalEXT" );

    // These are necessary for basic functioning.
    glGenBuffers    = (PFNGLGENBUFFERSPROC) glutGetProcAddress("glGenBuffers");
    glBindBuffer    = (PFNGLBINDBUFFERPROC) glutGetProcAddress("glBindBuffer");
    glBufferData    = (PFNGLBUFFERDATAPROC) glutGetProcAddress("glBufferData");
    glActiveTexture = (PFNGLACTIVETEXTUREPROC) glutGetProcAddress("glActiveTexture");

    if (!(glGenBuffers && glBindBuffer && glBufferData && glActiveTexture))
        return false; // fatal error

    // These are needed only when saving images or movies so the Visualizer can 
    // function without them.

    // Using the "EXT" names here means we only require OpenGL 2.0.
    glGenFramebuffersEXT    = (PFNGLGENFRAMEBUFFERSEXTPROC) glutGetProcAddress("glGenFramebuffersEXT");
    glGenRenderbuffersEXT   = (PFNGLGENRENDERBUFFERSEXTPROC) glutGetProcAddress("glGenRenderbuffersEXT");
    glBindFramebufferEXT    = (PFNGLBINDFRAMEBUFFEREXTPROC) glutGetProcAddress("glBindFramebufferEXT");
    glBindRenderbufferEXT   = (PFNGLBINDRENDERBUFFEREXTPROC) glutGetProcAddress("glBindRenderbufferEXT");
    glRenderbufferStorageEXT = (PFNGLRENDERBUFFERSTORAGEEXTPROC) glutGetProcAddress("glRenderbufferStorageEXT");
    glFramebufferRenderbufferEXT = (PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC) glutGetProcAddress("glFramebufferRenderbufferEXT");
    glDeleteRenderbuffersEXT = (PFNGLDELETERENDERBUFFERSEXTPROC) glutGetProcAddress("glDeleteRenderbuffersEXT");
    glDeleteFramebuffersEXT  = (PFNGLDELETEFRAMEBUFFERSEXTPROC) glutGetProcAddress("glDeleteFramebuffersEXT");

    if (! (glGenFramebuffersEXT  && glGenRenderbuffersEXT    && glBindFramebufferEXT
        && glBindRenderbufferEXT && glRenderbufferStorageEXT && glFramebufferRenderbufferEXT
        && glDeleteRenderbuffersEXT && glDeleteFramebuffersEXT))
        glCanSaveImages = false;
#endif

    return true;
}

