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
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon.h"
#include "simbody/internal/VisualizationEventListener.h"
#include "../src/VisualizationProtocol.h"
#include "lodepng.h"

#include <string>
#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <pthread.h>
#include <sys/stat.h>

// Get gl and glut using the appropriate platform-dependent incantations.
#if defined(__APPLE__)
    // OSX comes with a good glut implementation.
    #include <GLUT/glut.h>
#elif defined(_WIN32)
    #include "glut32/glut.h"    // we have our own private headers
    #include "glut32/glext.h" 

    // A Windows-only extension for disabling vsync, allowing unreasonably
    // high frame rates.
    PFNWGLSWAPINTERVALFARPROC wglSwapIntervalEXT;


    // These will hold the dynamically-determined function addresses.
    PFNGLGENBUFFERSPROC glGenBuffers;
    PFNGLBINDBUFFERPROC glBindBuffer;
    PFNGLBUFFERDATAPROC glBufferData;
    PFNGLCREATEPROGRAMPROC glCreateProgram;
    PFNGLCREATESHADERPROC glCreateShader;
    PFNGLSHADERSOURCEPROC glShaderSource;
    PFNGLCOMPILESHADERPROC glCompileShader;
    PFNGLATTACHSHADERPROC glAttachShader;
    PFNGLLINKPROGRAMPROC glLinkProgram;
    PFNGLUSEPROGRAMPROC glUseProgram;
    PFNGLUNIFORM3FPROC glUniform3f; 
    PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation; 
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
    // Linux: assume we have a good OpenGL 2.0 and working glut or free glut.
    #define GL_GLEXT_PROTOTYPES
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glext.h>
    #include <GL/glut.h>
#endif

static void initGlextFuncPointersIfNeeded();
static void redrawDisplay();
static void setKeepAlive(bool enable);
static void setVsync(bool enable);

// Next, get the functions necessary for reading from and writing to pipes.
#ifdef _WIN32
    #include <io.h>
    #define READ _read
#else
    #include <unistd.h>
    #define READ read
#endif

// gcc 4.4.3 complains bitterly if you don't check the return
// status from the write() system call. This avoids those 
// warnings and maybe, someday, will catch an error.
#define WRITE(pipeno, buf, len) \
   {int status=write((pipeno), (buf), (len)); \
    SimTK_ERRCHK4_ALWAYS(status!=-1, "VisualizationGUI",  \
    "An attempt to write() %d bytes to pipe %d failed with errno=%d (%s).", \
    (len),(pipeno),errno,strerror(errno));}

using namespace SimTK;
using namespace std;


static fTransform cameraTransform(fRotation(), fVec3(0, 0, 10));
static int inPipe, outPipe;

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
    Mesh(vector<float>& vertices, vector<float>& normals, vector<GLushort>& faces) : numVertices(vertices.size()/3), faces(faces) {
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
            glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_SHORT, &faces[0]);
        else if (representation == DecorativeGeometry::DrawPoints)
            glDrawArrays(GL_POINTS, 0, numVertices*3);
        else if (representation == DecorativeGeometry::DrawWireframe)
            glDrawElements(GL_LINES, edges.size(), GL_UNSIGNED_SHORT, &edges[0]);
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

static vector<Mesh*> meshes;

class RenderedMesh {
public:
    RenderedMesh(const fTransform& transform, const fVec3& scale, const fVec4& color, short representation, unsigned short meshIndex) :
            transform(transform), scale(scale), representation(representation), meshIndex(meshIndex) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
        this->color[3] = color[3];
    }
    void draw() {
        glPushMatrix();
        glTranslated(transform.T()[0], transform.T()[1], transform.T()[2]);
        fVec4 rot = transform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        if (representation == DecorativeGeometry::DrawSurface)
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
        else
            glColor3fv(color);
        meshes[meshIndex]->draw(representation);
        glPopMatrix();
    }
    const fTransform& getTransform() const {
        return transform;
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
        meshes[meshIndex]->getBoundingSphere(radius, center);
        center += transform.T();
        radius *= max(abs(scale[0]), max(abs(scale[1]), abs(scale[2])));
    }
private:
    fTransform transform;
    fVec3 scale;
    GLfloat color[4];
    short representation;
    unsigned short meshIndex;
};

class RenderedLine {
public:
    RenderedLine(const fVec3& color, float thickness) 
    :   color(color), thickness(thickness) {}

    void draw() {
        glColor3d(color[0], color[1], color[2]);
        glLineWidth(thickness);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glVertexPointer(3, GL_FLOAT, 0, &lines[0]);
        glDrawArrays(GL_LINES, 0, lines.size()/3);
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
    RenderedText(const fVec3& position, float scale, const fVec3& color, const string& text) :
            position(position), scale(scale/119), text(text) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
    }
    void draw() {
        glPushMatrix();
        glTranslated(position[0], position[1], position[2]);
        fVec4 rot = cameraTransform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale, scale, scale);
        glColor3fv(color);
        for (int i = 0; i < (int) text.size(); i++)
            glutStrokeCharacter(GLUT_STROKE_ROMAN, text[i]);
        glPopMatrix();
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
        center = position;
        radius = glutStrokeLength(GLUT_STROKE_ROMAN, (unsigned char*) text.c_str())*scale;
    }
private:
    fVec3 position;
    float scale;
    GLfloat color[3];
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
    Scene() : sceneHasBeenDrawn(false) {}

    vector<RenderedMesh> drawnMeshes;
    vector<RenderedMesh> solidMeshes;
    vector<RenderedMesh> transparentMeshes;
    vector<RenderedLine> lines;
    vector<RenderedText> strings;

    bool sceneHasBeenDrawn;
};

static Scene* scene=0;      // This is the front scene.

// This lock must be held whenever the front scene is being
// modified or drawn, or when it is being swapped with the back scene.
static pthread_mutex_t sceneLock;

// Wait on this if you want to be notified when the scene has been drawn.
static pthread_cond_t  sceneHasBeenDrawn;


// Some commands received by the listener must be executed on the rendering
// thread. They are saved in concrete objects derived from this abstract
// class.
class PendingCommand {
public:
    virtual void execute() = 0;
};

static int viewWidth, viewHeight;
static GLfloat fieldOfView = GLfloat(SimTK_PI/4);
static GLfloat nearClip = 1;
static GLfloat farClip = 1000;
static GLfloat groundHeight = 0;
static int groundAxis = 1;
static bool showGround = false, showShadows = true, showFPS = true;
static vector<PendingCommand*> pendingCommands;
static float fps = 0.0f;
static int fpsBaseTime = 0, fpsCounter = 0, nextMeshIndex;
static int frameCounter = 0;

static string overlayMessage;
static bool displayOverlayMessage = false;

class PendingMesh : public PendingCommand {
public:
    PendingMesh() {
        index = nextMeshIndex++;
    }
    void execute() {
        if ((int)meshes.size() <= index)
            meshes.resize(index+1);
        meshes[index] = new Mesh(vertices, normals, faces);
    }
    vector<float> vertices;
    vector<float> normals;
    vector<GLushort> faces;
    int index;
};

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
    for (int i = 0; i < (int) scene->strings.size(); i++) {
        fVec3 center;
        float radius;
        scene->strings[i].computeBoundingSphere(radius, center);
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

static void zoomCameraToShowWholeScene() {
    float radius;
    fVec3 center;
    computeSceneBounds(scene, radius, center);
    float viewDistance = radius/tan(min(fieldOfView, fieldOfView*viewWidth/viewHeight)/2);
    cameraTransform.updT() = center+cameraTransform.R()*fVec3(0, 0, viewDistance);
}

class PendingCameraZoom : public PendingCommand {
public:
    void execute() {
        zoomCameraToShowWholeScene();
    }
};

// This is the handler for simulator-defined menus. All it does is send back
// to the simulator the integer index that was associated with the particular
// menu entry on which the user clicked.
void menuSelected(int option) {
    char command = MENU_SELECTED;
    WRITE(outPipe, &command, 1);
    WRITE(outPipe, &option, sizeof(int));
}

class Menu {
public:
    Menu(string title, vector<pair<string, int> > items, void(*handler)(int)) 
    :   title(title), items(items), handler(handler), hasCreated(false) {}

    int draw(int x, int y) {
        if (!hasCreated) {
            id = glutCreateMenu(handler);
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
                glutSetMenu(firstNewSubmenu == 0 ? id : submenuIds[firstNewSubmenu-1]);
                for (int j = firstNewSubmenu; j < (int) components.size()-1; j++) {
                    glutAddSubMenu(components[j].c_str(), submenuIds[j]);
                    glutSetMenu(submenuIds[j]);
                }
                glutAddMenuEntry(components[components.size()-1].c_str(), items[i].second);
            }
            hasCreated = true;
        }
        minx = x;
        miny = y-18;
        int width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char*) title.c_str());
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
            if (currentMenu != id) {
                glutSetMenu(id);
                glutAttachMenu(GLUT_LEFT_BUTTON);
                currentMenu = id;
            }
        }
        else if (currentMenu == id) {
            glutSetMenu(id);
            glutDetachMenu(GLUT_LEFT_BUTTON);
            currentMenu = -1;
        }
    }
private:
    string title;
    vector<pair<string, int> > items;
    void(*handler)(int);
    int id, minx, miny, maxx, maxy;
    bool hasCreated;
    static int currentMenu;
};

int Menu::currentMenu = -1;

static void drawGroundAndSky(float farClipDistance) {
    static GLuint skyProgram = 0;
    if (skyProgram == 0) {
        const GLchar* vertexShaderSource =
        "varying vec3 position;\n"
        "uniform vec3 cameraPosition;\n"
        "void main() {\n"
        "gl_Position = ftransform();\n"
        "position = gl_Vertex.xyz-cameraPosition;\n"
        "}";
        const GLchar* fragmentShaderSource =
        "varying vec3 position;\n"
        "uniform vec3 upDirection;\n"
        "void main() {\n"
        "float gradient = 1.0-dot(normalize(position), upDirection);\n"
        "gradient *= gradient*gradient;\n"
        "gl_FragColor = clamp(vec4(gradient, 0.97*gradient, 1.0, 1.0), 0.0, 1.0);\n"
        "}";
        skyProgram = glCreateProgram();
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(vertexShader);
        glCompileShader(fragmentShader);
        glAttachShader(skyProgram, vertexShader);
        glAttachShader(skyProgram, fragmentShader);
        glLinkProgram(skyProgram);
    }
    static GLuint groundProgram = 0;
    if (groundProgram == 0) {
        const GLchar* vertexShaderSource =
        "varying vec3 position;\n"
        "uniform vec3 sdirection, tdirection, cameraPosition;\n"
        "void main() {\n"
        "gl_Position = ftransform();\n"
        "vec3 pos = (gl_ModelViewMatrix*gl_Vertex).xyz+cameraPosition;\n"
        "position = vec3(dot(sdirection, pos), 0.0, dot(tdirection, pos));\n"
        "}";
        const GLchar* fragmentShaderSource =
        "varying vec3 position;\n"
        "uniform vec3 color1, color2;\n"
        "void main() {\n"
        "vec2 square = floor(0.2*position.xz);\n"
        "vec2 delta = 0.2*position.xz-square.xy;\n"
        "float line = min(min(min(delta.x, delta.y), 1.0-delta.x), 1.0-delta.y);\n"
        "float blur = max(fwidth(position.x), fwidth(position.z));\n"
        "float pattern = 0.35;\n"
        "if (blur < 1.0)\n"
        "pattern += 0.5*noise1(6.0*position);\n"
        "gl_FragColor = vec4(mix(color1, color2, sqrt(line)*pattern), 1.0);\n"
        "}";
        groundProgram = glCreateProgram();
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
        glCompileShader(vertexShader);
        glCompileShader(fragmentShader);
        glAttachShader(groundProgram, vertexShader);
        glAttachShader(groundProgram, fragmentShader);
        glLinkProgram(groundProgram);
    }

    // Draw the rectangle to represent the sky.

    float viewDistance = 0.9999f*farClipDistance;
    float xwidth = viewDistance*tan(fieldOfView*viewWidth/viewHeight/2);
    float ywidth = viewDistance*tan(fieldOfView/2);
    fVec3 center = cameraTransform.T()-cameraTransform.R()*fVec3(0, 0, viewDistance);
    fVec3 corner1 = center+cameraTransform.R()*fVec3(-xwidth, -ywidth, 0);
    fVec3 corner2 = center+cameraTransform.R()*fVec3(xwidth, -ywidth, 0);
    fVec3 corner3 = center+cameraTransform.R()*fVec3(xwidth, ywidth, 0);
    fVec3 corner4 = center+cameraTransform.R()*fVec3(-xwidth, ywidth, 0);
    fVec3 cameraPosition, upDirection;
    if (groundAxis == 0) {
        cameraPosition = fVec3(groundHeight, cameraTransform.T()[1], cameraTransform.T()[2]);
        upDirection = fVec3(1, 0, 0);
    }
    else if (groundAxis == 1) {
        cameraPosition = fVec3(cameraTransform.T()[0], groundHeight, cameraTransform.T()[2]);
        upDirection = fVec3(0, 1, 0);
    }
    else {
        cameraPosition = fVec3(cameraTransform.T()[0], cameraTransform.T()[1], groundHeight);
        upDirection = fVec3(0, 0, 1);
    }
    glUseProgram(skyProgram);
    glUniform3f(glGetUniformLocation(skyProgram, "cameraPosition"), (GLfloat) cameraPosition[0], (GLfloat) cameraPosition[1], (GLfloat) cameraPosition[2]);
    glUniform3f(glGetUniformLocation(skyProgram, "upDirection"), (GLfloat) upDirection[0], (GLfloat) upDirection[1], (GLfloat) upDirection[2]);
    glDepthMask(GL_FALSE);
    glBegin(GL_QUADS);
    glVertex3d(corner1[0], corner1[1], corner1[2]);
    glVertex3d(corner2[0], corner2[1], corner2[2]);
    glVertex3d(corner3[0], corner3[1], corner3[2]);
    glVertex3d(corner4[0], corner4[1], corner4[2]);
    glEnd();
    glDepthMask(GL_TRUE);

    // Draw the ground plane.

    center[1] = 0;
    corner1 = center+fVec3(-farClipDistance, 0, -farClipDistance);
    corner2 = center+fVec3(farClipDistance, 0, -farClipDistance);
    corner3 = center+fVec3(farClipDistance, 0, farClipDistance);
    corner4 = center+fVec3(-farClipDistance, 0, farClipDistance);
    glUseProgram(groundProgram);
    Mat<4, 4, GLfloat> transform(1.0f);
    fVec3 sdir, tdir;
    if (groundAxis == 0) {
        transform[0][0] = transform[1][1] = 0.0f;
        transform[0][1] = transform[1][0] = 1.0f;
        sdir = fVec3(0, 1, 0);
        tdir = fVec3(0, 0, 1);
    }
    else if (groundAxis == 1) {
        sdir = fVec3(1, 0, 0);
        tdir = fVec3(0, 0, 1);
    }
    else {
        transform[1][1] = transform[2][2] = 0.0f;
        transform[1][2] = transform[2][1] = 1.0f;
        sdir = fVec3(1, 0, 0);
        tdir = fVec3(0, 1, 0);
    }
    sdir = ~cameraTransform.R()*sdir;
    tdir = ~cameraTransform.R()*tdir;
    cameraPosition = ~cameraTransform.R()*cameraPosition;
    glUniform3f(glGetUniformLocation(groundProgram, "sdirection"), (GLfloat) sdir[0], (GLfloat) sdir[1], (GLfloat) sdir[2]);
    glUniform3f(glGetUniformLocation(groundProgram, "tdirection"), (GLfloat) tdir[0], (GLfloat) tdir[1], (GLfloat) tdir[2]);
    glUniform3f(glGetUniformLocation(groundProgram, "cameraPosition"), (GLfloat) cameraPosition[0], (GLfloat) cameraPosition[1], (GLfloat) cameraPosition[2]);
    glUniform3f(glGetUniformLocation(groundProgram, "color1"), 0.3f, 0.2f, 0.0f);
    transform[groundAxis][3] = groundHeight;
    if (showShadows) {
        // Draw shadows on the ground.

        Mat<4, 4, GLfloat> transform2(1.0f);
        transform2[0][1] = 0.2f;
        transform2[1][1] = 0.0f;
        transform2[2][1] = 0.2f;
        transform2 = transform*transform2;
        glPushMatrix();
        glMultMatrixf(&transform2[0][0]);
        glDepthRange(0.0, 0.9999);
        // Solid and transparent shadows are the same color (sorry). Trying to
        // mix light and dark shadows is much harder and any simple attempts
        // (e.g. put light shadows on top of dark ones) look terrible.
        glUniform3f(glGetUniformLocation(groundProgram, "color2"), 0.5f, 0.4f, 0.35f);
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
            scene->solidMeshes[i].draw();
        for (int i = 0; i < (int) scene->transparentMeshes.size(); i++)
            scene->transparentMeshes[i].draw();
        glPopMatrix();
        glDepthRange(0.0, 1.0);
    }
    glUniform3f(glGetUniformLocation(groundProgram, "color2"), 1.0f, 0.8f, 0.7f);
    glDisable(GL_CULL_FACE);
    glPushMatrix();
    glMultMatrixf(&transform[0][0]);
    glBegin(GL_QUADS);
    glVertex3d(corner1[0], corner1[1], corner1[2]);
    glVertex3d(corner2[0], corner2[1], corner2[2]);
    glVertex3d(corner3[0], corner3[1], corner3[2]);
    glVertex3d(corner4[0], corner4[1], corner4[2]);
    glEnd();
    glEnable(GL_CULL_FACE);
    glPopMatrix();
    glUseProgram(0);
}

static vector<Menu> menus;


static void changeSize(int width, int height) {
    if (height == 0)
        height = 1;
    viewWidth = width;
    viewHeight = height;
}

static void renderScene() {
    static bool firstTime = true;
    static GLfloat prevNearClip; // initialize to near & farClip
    static GLfloat prevFarClip;  
    if (firstTime) {
        prevNearClip = nearClip;
        prevFarClip = farClip;
        firstTime = false;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glDisableClientState(GL_NORMAL_ARRAY);

    pthread_mutex_lock(&sceneLock);             //-------- LOCK SCENE --------

    if (scene != NULL) {
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
        computeSceneBounds(scene, sceneRadius, sceneCenter);
        float sceneScale = std::max(0.1f, sceneRadius);

        float centerDepth = ~(cameraTransform.T()-sceneCenter)*(cameraTransform.R().col(2));
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
        fVec3 cameraPos = cameraTransform.T();
        fVec3 centerPos = cameraTransform.T()+cameraTransform.R()*fVec3(0, 0, -1);
        fVec3 upDir = cameraTransform.R()*fVec3(0, 1, 0);
        gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2], centerPos[0], centerPos[1], centerPos[2], upDir[0], upDir[1], upDir[2]);

        // Render the objects in the scene.

        for (int i = 0; i < (int) scene->lines.size(); i++)
            scene->lines[i].draw();
        glLineWidth(2);
        for (int i = 0; i < (int) scene->strings.size(); i++)
            scene->strings[i].draw();
        glLineWidth(1);
        glEnableClientState(GL_NORMAL_ARRAY);
        for (int i = 0; i < (int) scene->drawnMeshes.size(); i++)
            scene->drawnMeshes[i].draw();
        glEnable(GL_LIGHTING);
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
            scene->solidMeshes[i].draw();
        if (showGround)
            drawGroundAndSky(farClipDistance);
        glEnable(GL_BLEND);
        glDepthMask(GL_FALSE);
        vector<pair<float, int> > order(scene->transparentMeshes.size());
        for (int i = 0; i < (int) order.size(); i++)
            order[i] = make_pair((float)(~cameraTransform.R()*scene->transparentMeshes[i].getTransform().T())[2], i);
        sort(order.begin(), order.end());
        for (int i = 0; i < (int) order.size(); i++)
            scene->transparentMeshes[order[i].second].draw();

        // Update the frame rate counter.

        int time = glutGet(GLUT_ELAPSED_TIME);
        if (time-fpsBaseTime > 1000) {
            fps = 1000.0f*fpsCounter/(time-fpsBaseTime);
            fpsBaseTime = time;
            fpsCounter = 0;
        }

        scene->sceneHasBeenDrawn = true;
    }

    // Signal in case the listener is waiting.
    pthread_cond_signal(&sceneHasBeenDrawn);    //----- SIGNAL CONDITION -----
    pthread_mutex_unlock(&sceneLock);           //----- UNLOCK SCENE ---------
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
    renderScene();

    // Draw menus.

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

    // Draw the frame rate counter.

    if (showFPS) {
        char fpstxt[64]; sprintf(fpstxt, "FPS: %.1f", fps); // 1 decimal place
        char cnttxt[64]; sprintf(cnttxt, "#%d", frameCounter);
        glColor3f(1.0f, 0.5f, 0.0f);
        glRasterPos2f(10, 25);
        for (const char* p = fpstxt; *p; ++p)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
        glRasterPos2f(10, 25+18);
        for (const char* p = cnttxt; *p; ++p)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
    }

    // Draw a message overlay.

    if (displayOverlayMessage) {
        glColor3f(1.0f, 0.5f, 0.0f);
        int width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char*) overlayMessage.c_str());
        glRasterPos2f(GLfloat(max(0, (viewWidth-width)/2)), GLfloat(viewHeight/2));
        for (int i = 0; i < (int) overlayMessage.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, overlayMessage[i]);
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glutSwapBuffers();
    ++fpsCounter;
    ++frameCounter;
}

// These are set when a mouse button is clicked and then referenced
// while the mouse is dragged with that button down.
static int clickModifiers;
static int clickButton;
static int clickX;
static int clickY;
static bool rotateCenterValid = false;
static fVec3 rotateCenter;
static float sceneScale; // for scaling translations

// Handle the initial press of a mouse button. Standard glut does not
// deal with the mouse wheel, but freeglut and some tweaked gluts treat
// it as a "button press" with button codes 3 (up) and 4 (down). This
// function handles those properly if they are present.
static void mouseButtonPressedOrReleased(int button, int state, int x, int y) {
    const int GlutWheelUp = 3, GlutWheelDown = 4;

    float radius;
    fVec3 sceneCenter;
    computeSceneBounds(scene, radius, sceneCenter);
    sceneScale = std::max(radius, 0.1f);

    // "state" (pressed/released) is irrelevant for mouse wheel. However, if 
    // we're being called by freeglut we'll get called twice, while (patched) 
    // glut calls just once, with state=GLUT_UP. 
    if ((button == GlutWheelUp || button == GlutWheelDown)
        && (state==GLUT_UP)) // TODO: does this work OK on Mac GLUT?
    {
        // Scroll wheel.
        const float ZoomFractionPerWheelClick = 0.15f;   // 15% scene radius

        const int   direction = button==GlutWheelUp ? -1 : 1;
        const float zoomBy    = direction * (ZoomFractionPerWheelClick * sceneScale);

        pthread_mutex_lock(&sceneLock);         //------ LOCK SCENE ----------
        cameraTransform.updT() += cameraTransform.R()*fVec3(0, 0, zoomBy);
        glutPostRedisplay();                    //------ POST REDISPLAY ------
        pthread_mutex_unlock(&sceneLock);       //------ UNLOCK SCENE --------
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

    // Left button is rotation; when it is first pressed we calcuate the center
    // of rotation; we'll do the actual rotating in mouseDragged().
    if (clickButton == GLUT_LEFT_BUTTON) {
        float distToCenter = (sceneCenter-cameraTransform.T()).norm();
        float defaultDistance = sceneScale;
        if (distToCenter > defaultDistance)
            rotateCenter = sceneCenter;
        else {
            fVec3 cameraDir = cameraTransform.R()*fVec3(0, 0, -1);
            fVec3 lookAt = cameraTransform.T()+defaultDistance*cameraDir;
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
        pthread_mutex_lock(&sceneLock);         //------ LOCK SCENE ----------
        cameraTransform.updT() += translatePerPixel*cameraTransform.R()*fVec3(dx, -dy, 0);
        pthread_mutex_unlock(&sceneLock);       //------ UNLOCK SCENE --------
    }

    // zoom: middle button or alt-left button (or mouse wheel; see above)
    else if (  clickButton == GLUT_MIDDLE_BUTTON 
           || (clickButton == GLUT_LEFT_BUTTON && clickModifiers & GLUT_ACTIVE_ALT))
    {
        pthread_mutex_lock(&sceneLock);         //------ LOCK SCENE ----------
        cameraTransform.updT() += translatePerPixel* cameraTransform.R()*fVec3(0, 0, dy);
        pthread_mutex_unlock(&sceneLock);       //------ UNLOCK SCENE --------
    }

    // rotate: left button alone: rotate scene left/right or up/down
    //         ctrl-left button:  roll scene about camera direction
    else if (clickButton == GLUT_LEFT_BUTTON && rotateCenterValid) {
        fVec3 cameraPos = cameraTransform.T();
        fVec3 cameraDir = cameraTransform.R()*fVec3(0, 0, -1);
        fVec3 upDir = cameraTransform.R()*fVec3(0, 1, 0);
        fRotation r;
        if (clickModifiers & GLUT_ACTIVE_CTRL)
            r.setRotationFromAngleAboutAxis(AnglePerPixel*(dy-dx), ZAxis);
        else
            r.setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence, 
                AnglePerPixel*dy, XAxis, AnglePerPixel*dx, YAxis);
        r = cameraTransform.R()*r*~cameraTransform.R();
        cameraPos = r*(cameraPos-rotateCenter)+rotateCenter;
        cameraDir = r*cameraDir;
        upDir = r*upDir;

        pthread_mutex_lock(&sceneLock);         //------ LOCK SCENE ----------
        cameraTransform.updT() = cameraPos;
        cameraTransform.updR().setRotationFromTwoAxes(fUnitVec3(-cameraDir), ZAxis, upDir, YAxis);
        pthread_mutex_unlock(&sceneLock);       //------ UNLOCK SCENE --------
    }
    else
        return;

    // Something changed.

    clickX = x; clickY = y;
    glutPostRedisplay();                    //-------- POST REDISPLAY --------
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
    char command = KEY_PRESSED;
    WRITE(outPipe, &command, 1);
    unsigned char buffer[2];
    buffer[0] = key;
    buffer[1] = 0;
    unsigned char modifiers = glutGetModifiers();
    if ((modifiers & GLUT_ACTIVE_SHIFT) != 0)
        buffer[1] += VisualizationEventListener::ShiftIsDown;
    if ((modifiers & GLUT_ACTIVE_CTRL) != 0)
        buffer[1] += VisualizationEventListener::ControlIsDown;
    if ((modifiers & GLUT_ACTIVE_ALT) != 0)
        buffer[1] += VisualizationEventListener::AltIsDown;
    WRITE(outPipe, buffer, 2);
}

static void specialKeyPressed(int key, int x, int y) {
    char command = KEY_PRESSED;
    WRITE(outPipe, &command, 1);
    unsigned char buffer[2];
    buffer[0] = (unsigned char)key; // this is the special key code
    buffer[1] = VisualizationEventListener::IsSpecialKey;
    unsigned char modifiers = glutGetModifiers();
    if ((modifiers & GLUT_ACTIVE_SHIFT) != 0)
        buffer[1] &= VisualizationEventListener::ShiftIsDown;
    if ((modifiers & GLUT_ACTIVE_CTRL) != 0)
        buffer[1] &= VisualizationEventListener::ControlIsDown;
    if ((modifiers & GLUT_ACTIVE_ALT) != 0)
        buffer[1] &= VisualizationEventListener::AltIsDown;
    WRITE(outPipe, buffer, 2);
}

static void disableOverlayTimer(int value) {
    displayOverlayMessage = false;
    glutPostRedisplay();                    //-------- POST REDISPLAY --------
}

static void setOverlayMessage(const string& message) {
    overlayMessage = message;
    displayOverlayMessage = true;
    glutTimerFunc(5000, disableOverlayTimer, 0);
}

static void saveImage() {
    // Create offscreen buffers for rendering the image.

    GLuint frameBuffer, colorBuffer, depthBuffer;
    glGenFramebuffersEXT(1, &frameBuffer);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBuffer);
    glGenRenderbuffersEXT(1, &colorBuffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, colorBuffer);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGB8, viewWidth, viewHeight);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, colorBuffer);
    glGenRenderbuffersEXT(1, &depthBuffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depthBuffer);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24, viewWidth, viewHeight);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depthBuffer);

    // Render the image and load it into memory.

    renderScene();
    vector<unsigned char> data(viewWidth*viewHeight*3);
    glReadPixels(0, 0, viewWidth, viewHeight, GL_RGB, GL_UNSIGNED_BYTE, &data[0]);
    glDeleteRenderbuffersEXT(1, &colorBuffer);
    glDeleteRenderbuffersEXT(1, &depthBuffer);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glDeleteFramebuffersEXT(1, &frameBuffer);

    // Flip the image vertically, since OpenGL and PNG use different row orders.

    const int rowLength = 3*viewWidth;
    for (int row = 0; row < viewHeight/2; ++row) {
        const int base1 = row*rowLength;
        const int base2 = (viewHeight-1-row)*rowLength;
        for (int i = 0; i < rowLength; i++) {
            unsigned char temp = data[base1+i];
            data[base1+i] = data[base2+i];
            data[base2+i] = temp;
        }
    }

    // Save it to disk.

    struct stat statInfo;
    int counter = 0;
    string filename;
    do {
        counter++;
        stringstream namestream;
        namestream << "Saved Image ";
        namestream << counter;
        namestream << ".png";
        filename = namestream.str();
    } while (stat(filename.c_str(), &statInfo) == 0);
    LodePNG::encode(filename, data, viewWidth, viewHeight, 2, 8);
    setOverlayMessage("Saved as: "+filename);
    redrawDisplay();
}

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

static void makeBox()  {
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
    meshes.push_back(new Mesh(vertices, normals, faces));
}

static void makeSphere() {
    const int numLatitude = 8;
    const int numLongitude = 12;
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
    meshes.push_back(new Mesh(vertices, normals, faces));
}

static void makeCylinder() {
    const int numSides = 12;
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
    meshes.push_back(new Mesh(vertices, normals, faces));
}

static void makeCircle() {
    const int numSides = 12;
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
    meshes.push_back(new Mesh(vertices, normals, faces));
}

// Read a particular number of bytes from the inPipe to the given data.
// This will hang until the expected number of bytes has been received.
static void readData(char* buffer, int bytes) {
    int totalRead = 0;
    while (totalRead < bytes)
        totalRead += READ(inPipe, buffer+totalRead, bytes-totalRead);
}

// We have just processed a START_OF_SCENE command. Read in all the scene
// elements until we see an END_OF_SCENE command. We allocate a new Scene
// object to hold the scene and return a pointer to it. Don't forget to
// delete that object when you are done with it.
static Scene* readNewScene() {
    char buffer[256];
    float*          floatBuffer = (float*)          buffer;
    int*            intBuffer   = (int*)            buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;

    Scene* newScene = new Scene;

    bool finished = false;
    while (!finished) {
        readData(buffer, 1);
        char command = buffer[0];

        switch (command) {

        case END_OF_SCENE:
            finished = true;
            break;

        // Add a scene element that uses an already-cached mesh.
        case ADD_POINT_MESH:
        case ADD_WIREFRAME_MESH:
        case ADD_SOLID_MESH: {
            readData(buffer, 13*sizeof(float)+sizeof(short));
            fTransform position;
            position.updR().setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
            position.updT() = fVec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 scale = fVec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            fVec4 color = fVec4(floatBuffer[9], floatBuffer[10], floatBuffer[11], floatBuffer[12]);
            short representation = (command == ADD_POINT_MESH ? DecorativeGeometry::DrawPoints : (command == ADD_WIREFRAME_MESH ? DecorativeGeometry::DrawWireframe : DecorativeGeometry::DrawSurface));
            RenderedMesh mesh(position, scale, color, representation, shortBuffer[13*sizeof(float)/sizeof(short)]);
            if (command != ADD_SOLID_MESH)
                newScene->drawnMeshes.push_back(mesh);
            else if (color[3] == 1)
                newScene->solidMeshes.push_back(mesh);
            else
                newScene->transparentMeshes.push_back(mesh);
            break;
        }

        case ADD_LINE: {
            readData(buffer, 10*sizeof(float));
            fVec3 color = fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            float thickness = floatBuffer[3];
            int index;
            int numLines = newScene->lines.size();
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

        case ADD_TEXT: {
            readData(buffer, 7*sizeof(float)+sizeof(short));
            fVec3 position = fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            float scale = floatBuffer[3];
            fVec3 color = fVec3(floatBuffer[4], floatBuffer[5], floatBuffer[6]);
            short length = shortBuffer[7*sizeof(float)/sizeof(short)];
            readData(buffer, length);
            newScene->strings.push_back(RenderedText(position, scale, color, string(buffer, length)));
            break;
        }

        case ADD_FRAME: {
            readData(buffer, 10*sizeof(float));
            fRotation rotation;
            rotation.setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
            fVec3 position(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            float axisLength = floatBuffer[6];
            float textScale = 0.2f*axisLength;
            float lineThickness = 1;
            fVec3 color = fVec3(floatBuffer[7], floatBuffer[8], floatBuffer[9]);
            int index;
            int numLines = newScene->lines.size();
            for (index = 0; index < numLines && (color != newScene->lines[index].getColor() || newScene->lines[index].getThickness() != lineThickness); index++)
                ;
            if (index == numLines)
                newScene->lines.push_back(RenderedLine(color, lineThickness));
            vector<GLfloat>& line = newScene->lines[index].getLines();
            fVec3 end = position+rotation*fVec3(axisLength, 0, 0);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->strings.push_back(RenderedText(end, textScale, color, "X"));
            end = position+rotation*fVec3(0, axisLength, 0);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->strings.push_back(RenderedText(end, textScale, color, "Y"));
            end = position+rotation*fVec3(0, 0, axisLength);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->strings.push_back(RenderedText(end, textScale, color, "Z"));
            break;
        }

        // Define a new mesh that will be assigned the next available mesh
        // index. It will be cached here and then can be referenced in this
        // scene and others by using it mesh index.
        case DEFINE_MESH: {
            readData(buffer, 2*sizeof(short));
            PendingMesh* mesh = new PendingMesh(); // assigns next mesh index
            int numVertices = shortBuffer[0];
            int numFaces = shortBuffer[1];
            mesh->vertices.resize(3*numVertices, 0);
            mesh->normals.resize(3*numVertices);
            mesh->faces.resize(3*numFaces);
            readData((char*) &mesh->vertices[0], mesh->vertices.size()*sizeof(float));
            readData((char*) &mesh->faces[0], mesh->faces.size()*sizeof(short));

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
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE --------
            pendingCommands.insert(pendingCommands.begin(), mesh);
            pthread_mutex_unlock(&sceneLock);   //------ UNLOCK SCENE --------
            break;
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
void* listenForInput(void* args) {
    char buffer[256];
    float* floatBuffer = (float*) buffer;
    int* intBuffer = (int*) buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;
    while (true) {
        // Read commands from the simulator.
        readData(buffer, 1);
        switch (buffer[0]) {
        case SET_CAMERA: {
            readData(buffer, 6*sizeof(float));
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            cameraTransform.updR().setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
            cameraTransform.updT() = fVec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case ZOOM_CAMERA: {
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            pendingCommands.push_back(new PendingCameraZoom());
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case LOOK_AT: {
            readData(buffer, 6*sizeof(float));
            fVec3 point(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            fVec3 updir(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            cameraTransform.updR().setRotationFromTwoAxes(fUnitVec3(cameraTransform.T()-point), ZAxis, updir, YAxis);
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case SET_FIELD_OF_VIEW: {
            readData(buffer, sizeof(float));
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            fieldOfView = floatBuffer[0];
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case SET_CLIP_PLANES: {
            readData(buffer, 2*sizeof(float));
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            nearClip = floatBuffer[0];
            farClip = floatBuffer[1];
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case SET_GROUND_POSITION: {
            readData(buffer, sizeof(float)+sizeof(short));
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            groundHeight = floatBuffer[0];
            groundAxis = shortBuffer[sizeof(float)/sizeof(short)];
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case DEFINE_MENU: {
            readData(buffer, sizeof(short));
            int titleLength = shortBuffer[0];
            vector<char> titleBuffer(titleLength);
            readData(&titleBuffer[0], titleLength);
            string title(&titleBuffer[0], titleLength);
            readData(buffer, sizeof(short));
            int numItems = shortBuffer[0];
            vector<pair<string, int> > items(numItems);
            for (int index = 0; index < numItems; index++) {
                readData(buffer, 2*sizeof(int));
                items[index].second = intBuffer[0];
                vector<char> textBuffer(intBuffer[1]);
                readData(&textBuffer[0], intBuffer[1]);
                items[index].first = string(&textBuffer[0], intBuffer[1]);
            }
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            menus.push_back(Menu(title, items, menuSelected));
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }
        case START_OF_SCENE: {
            Scene* newScene = readNewScene();
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            if (scene != NULL) {
                while (!scene->sceneHasBeenDrawn) {
                    // -------- WAIT FOR CONDITION --------
                    // Give up the lock and wait for notice.
                    pthread_cond_wait(&sceneHasBeenDrawn,
                                        &sceneLock);
                    // Now we hold the lock again, but this might
                    // be a spurious wakeup so check again.
                }
                // Previous scene has been drawn.
                delete scene; scene = 0;
            }
            // Swap in the new scene.
            scene = newScene;
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            break;
        }

        default:
            SimTK_ASSERT_ALWAYS(false, "Unexpected data sent to visualizer");
        }

        // Do this after every received command.
        glutPostRedisplay();                    //------- POST REDISPLAY -----
    }
    return 0;
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

static const int MENU_SAVE_IMAGE = 11;

static const int MENU_VSYNC_ENABLE = 12;
static const int MENU_VSYNC_DISABLE = 13;
static const int MENU_KEEPALIVE_ENABLE = 14;
static const int MENU_KEEPALIVE_DISABLE = 15;

// This is the handler for our built-in "View" pull down menu.
void viewMenuSelected(int option) {
    fRotation groundRotation;
    if (groundAxis == 0)
        groundRotation.setRotationFromTwoAxes(fUnitVec3(0, 1, 0), ZAxis, fVec3(1, 0, 0), YAxis);
    else if (groundAxis == 2)
        groundRotation.setRotationFromTwoAxes(fUnitVec3(1, 0, 0), ZAxis, fVec3(0, 0, 1), YAxis);
    
    switch (option) {
    case MENU_VIEW_FRONT:
        cameraTransform.updR().setRotationToIdentityMatrix();
        cameraTransform.updR() = groundRotation*cameraTransform.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_BACK:
        cameraTransform.updR().setRotationFromAngleAboutY((float)SimTK_PI);
        cameraTransform.updR() = groundRotation*cameraTransform.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_LEFT:
        cameraTransform.updR().setRotationFromAngleAboutY(-(float)(SimTK_PI/2));
        cameraTransform.updR() = groundRotation*cameraTransform.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_RIGHT:
        cameraTransform.updR().setRotationFromAngleAboutY((float)(SimTK_PI/2));
        cameraTransform.updR() = groundRotation*cameraTransform.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_TOP:
        cameraTransform.updR().setRotationFromAngleAboutX(-(float)(SimTK_PI/2));
        cameraTransform.updR() = groundRotation*cameraTransform.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_VIEW_BOTTOM:
        cameraTransform.updR().setRotationFromAngleAboutX((float)(SimTK_PI/2));
        cameraTransform.updR() = groundRotation*cameraTransform.R();
        zoomCameraToShowWholeScene();
        break;
    case MENU_BACKGROUND_BLACK:
        showGround = false;
        glClearColor(0, 0, 0, 1);
        break;
    case MENU_BACKGROUND_WHITE:
        showGround = false;
        glClearColor(1, 1, 1, 1);
        break;
    case MENU_BACKGROUND_SKY:
        showGround = true;
        break;
    case MENU_SHOW_SHADOWS:
        showShadows = !showShadows;
        break;
    case MENU_SHOW_FPS:
        showFPS = !showFPS;
        break;
    case MENU_SAVE_IMAGE:
        saveImage();
        break;
    case MENU_VSYNC_ENABLE: setVsync(true); break;
    case MENU_VSYNC_DISABLE: setVsync(false); break;
    case MENU_KEEPALIVE_ENABLE: setKeepAlive(true); break;
    case MENU_KEEPALIVE_DISABLE: setKeepAlive(false); break;
    }

    glutPostRedisplay();                    //-------- POST REDISPLAY --------
}

static const int DefaultWindowWidth  = 600;
static const int DefaultWindowHeight = 500;

// This seems to be necessary on Mac and Linux where the built-in glut
// idle hangs if there is no activity in the gl window. Not needed on Windows.
static void keepAliveIdleFunc() {
    struct timespec ts;
    // Wait 10ms and then check for something to do.
    nsToTimespec(secToNs(1./100), ts);
    nanosleep(&ts, 0);
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

int main(int argc, char** argv) {
    SimTK_ASSERT_ALWAYS(argc >= 3, "VisualizationGUI: must be at least two command line arguments (pipes)");

    stringstream(argv[1]) >> inPipe;
    stringstream(argv[2]) >> outPipe;

    string title("SimTK Visualizer");
    if (argc >= 4 && argv[3])
        title += ": " + string(argv[3]);


    // Initialize GLUT.

    glutInit(&argc, argv);

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

    printf( "OpenGL version information:\n"
            "---------------------------\n"
            "Vendor:   %s\n" 
            "Renderer: %s\n"
            "Version:  %s\n"
            "GLSL:     %s\n", 
        glGetString (GL_VENDOR), glGetString (GL_RENDERER),
        glGetString (GL_VERSION), glGetString (GL_SHADING_LANGUAGE_VERSION));

    // Set up callback funtions.
    glutDisplayFunc(redrawDisplay);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mouseButtonPressedOrReleased);
    glutMotionFunc(mouseDragged);
    glutPassiveMotionFunc(mouseMoved);
    glutKeyboardFunc(ordinaryKeyPressed);
    glutSpecialFunc(specialKeyPressed);

    // On some systems (Windows at least), some of the gl functions may 
    // need to be loaded dynamically.
    initGlextFuncPointersIfNeeded();

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
    glEnable(GL_RESCALE_NORMAL);
    makeBox();
    makeSphere();
    makeCylinder();
    makeCircle();
    nextMeshIndex = meshes.size();
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
    items.push_back(make_pair("Save Image", MENU_SAVE_IMAGE));
    items.push_back(make_pair("Vsync/Enable (default)", MENU_VSYNC_ENABLE));
    items.push_back(make_pair("Vsync/Disable", MENU_VSYNC_DISABLE));
    items.push_back(make_pair("Keep Alive/Enable", MENU_KEEPALIVE_ENABLE));
    items.push_back(make_pair("Keep Alive/Disable", MENU_KEEPALIVE_DISABLE));
    menus.push_back(Menu("View", items, viewMenuSelected));

    // Initialize pthread lock and condition variable.
    pthread_mutex_init(&sceneLock, NULL);
    pthread_cond_init(&sceneHasBeenDrawn, NULL);

    // Spawn the listener thread. After this it runs independently.
    pthread_t thread;
    pthread_create(&thread, NULL, listenForInput, NULL);

    // Avoid hangs on Mac & Linux.
#ifndef _WIN32
    setKeepAlive(true);
#endif

    // Enter the main loop. Note that there is no idle function. We expect
    // the main thread to go to sleep when there is nothing to do and that
    // the various events will call glutPostRedisplay() if anything happens
    // that requires a screen update.
    glutMainLoop();
    return 0;
}


// Initialize function pointers for Windows GL extensions.
static void initGlextFuncPointersIfNeeded() { 
#ifdef _WIN32
    wglSwapIntervalEXT = (PFNWGLSWAPINTERVALFARPROC)wglGetProcAddress( "wglSwapIntervalEXT" );

    glGenBuffers    = (PFNGLGENBUFFERSPROC) glutGetProcAddress("glGenBuffers");
    glBindBuffer    = (PFNGLBINDBUFFERPROC) glutGetProcAddress("glBindBuffer");
    glBufferData    = (PFNGLBUFFERDATAPROC) glutGetProcAddress("glBufferData");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) glutGetProcAddress("glCreateProgram");
    glCreateShader  = (PFNGLCREATESHADERPROC) glutGetProcAddress("glCreateShader");
    glShaderSource  = (PFNGLSHADERSOURCEPROC) glutGetProcAddress("glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) glutGetProcAddress("glCompileShader");
    glAttachShader  = (PFNGLATTACHSHADERPROC) glutGetProcAddress("glAttachShader");
    glLinkProgram   = (PFNGLLINKPROGRAMPROC) glutGetProcAddress("glLinkProgram");
    glUseProgram    = (PFNGLUSEPROGRAMPROC) glutGetProcAddress("glUseProgram");
    glUniform3f     = (PFNGLUNIFORM3FPROC) glutGetProcAddress("glUniform3f");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) glutGetProcAddress("glGetUniformLocation");
    // Using the "EXT" names here means we only require OpenGL 2.0.
    glGenFramebuffersEXT    = (PFNGLGENFRAMEBUFFERSEXTPROC) glutGetProcAddress("glGenFramebuffersEXT");
    glGenRenderbuffersEXT   = (PFNGLGENRENDERBUFFERSEXTPROC) glutGetProcAddress("glGenRenderbuffersEXT");
    glBindFramebufferEXT    = (PFNGLBINDFRAMEBUFFEREXTPROC) glutGetProcAddress("glBindFramebufferEXT");
    glBindRenderbufferEXT   = (PFNGLBINDRENDERBUFFEREXTPROC) glutGetProcAddress("glBindRenderbufferEXT");
    glRenderbufferStorageEXT = (PFNGLRENDERBUFFERSTORAGEEXTPROC) glutGetProcAddress("glRenderbufferStorageEXT");
    glFramebufferRenderbufferEXT = (PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC) glutGetProcAddress("glFramebufferRenderbufferEXT");
    glDeleteRenderbuffersEXT = (PFNGLDELETERENDERBUFFERSEXTPROC) glutGetProcAddress("glDeleteRenderbuffersEXT");
    glDeleteFramebuffersEXT  = (PFNGLDELETEFRAMEBUFFERSEXTPROC) glutGetProcAddress("glDeleteFramebuffersEXT");

    //TODO: check that all functions were present.
#endif
}

