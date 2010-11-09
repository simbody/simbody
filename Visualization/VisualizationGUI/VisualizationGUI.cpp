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
#include "simbody/internal/VisualizationProtocol.h"
#include "simbody/internal/VisualizationEventListener.h"
#include "lodepng.h"
#ifdef _WIN32
    // This file will not compile unless windows.h is
    // restricted by these defines due to name conflicts.
    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
    #include <windows.h>
    #include <io.h>
    #include <GL/gl.h>
    #include "glext.h"
    #include "glut.h"
    #define READ _read
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
    PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
    PFNGLGENRENDERBUFFERSPROC glGenRenderbuffers;
    PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
    PFNGLBINDRENDERBUFFERPROC glBindRenderbuffer;
    PFNGLRENDERBUFFERSTORAGEPROC glRenderbufferStorage;
    PFNGLFRAMEBUFFERRENDERBUFFERPROC glFramebufferRenderbuffer;
    PFNGLDELETERENDERBUFFERSPROC glDeleteRenderbuffers;
    PFNGLDELETEFRAMEBUFFERSPROC glDeleteFramebuffers;
#else
    #include <unistd.h>
    #ifdef __APPLE__
        #include <GLUT/glut.h>
    #else
        #define GL_GLEXT_PROTOTYPES
        #include <GL/gl.h>
        #include <GL/glext.h>
        #include <GL/glut.h>
    #endif
    #define READ read
#endif
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

using namespace SimTK;
using namespace std;

// gcc 4.4.3 complains bitterly if you don't check the return
// status from the write() system call. This avoids those 
// warnings and maybe, someday, will catch an error.
#define WRITE(pipeno, buf, len) \
   {int status=write((pipeno), (buf), (len)); \
    SimTK_ERRCHK4_ALWAYS(status!=-1, "VisualizationGUI",  \
    "An attempt to write() %d bytes to pipe %d failed with errno=%d (%s).", \
    (len),(pipeno),errno,strerror(errno));}

static Transform cameraTransform(Rotation(), Vec3(0, 0, 10));
static int clickModifiers;
static int clickButton;
static int clickX;
static int clickY;
static Vec3 rotateCenter;
static int inPipe, outPipe;
static bool needRedisplay;

static void computeBoundingSphereForVertices(const vector<float>& vertices, Real& radius, Vec3& center) {
    Vec3 lower(vertices[0], vertices[1], vertices[2]);
    Vec3 upper = lower;
    for (int i = 3; i < (int) vertices.size(); i += 3) {
        for (int j = 0; j < 3; j++) {
            lower[j] = min(lower[j], (Real) vertices[i+j]);
            upper[j] = max(upper[j], (Real) vertices[i+j]);
        }
    }
    center = 0.5*(lower+upper);
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
    void getBoundingSphere(Real& radius, Vec3& center) {
        radius = this->radius;
        center = this->center;
    }
private:
    int numVertices;
    GLuint vertBuffer, normBuffer;
    vector<GLushort> edges, faces;
    Vec3 center;
    Real radius;
};

static vector<Mesh*> meshes;

class RenderedMesh {
public:
    RenderedMesh(const Transform& transform, const Vec3& scale, const Vec4& color, short representation, unsigned short meshIndex) :
            transform(transform), scale(scale), representation(representation), meshIndex(meshIndex) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
        this->color[3] = color[3];
    }
    void draw() {
        glPushMatrix();
        glTranslated(transform.T()[0], transform.T()[1], transform.T()[2]);
        Vec4 rot = transform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        if (representation == DecorativeGeometry::DrawSurface)
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
        else
            glColor3fv(color);
        meshes[meshIndex]->draw(representation);
        glPopMatrix();
    }
    const Transform& getTransform() const {
        return transform;
    }
    void computeBoundingSphere(Real& radius, Vec3& center) {
        meshes[meshIndex]->getBoundingSphere(radius, center);
        center += transform.T();
        radius *= max(abs(scale[0]), max(abs(scale[1]), abs(scale[2])));
    }
private:
    Transform transform;
    Vec3 scale;
    GLfloat color[4];
    short representation;
    unsigned short meshIndex;
};

class RenderedLine {
public:
    RenderedLine(const Vec3& color, float thickness) : color(color), thickness(thickness) {
    }
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
    const Vec3& getColor() const {
        return color;
    }
    float getThickness() const {
        return thickness;
    }
    void computeBoundingSphere(Real& radius, Vec3& center) {
        computeBoundingSphereForVertices(lines, radius, center);
    }
private:
    Vec3 color;
    float thickness;
    vector<GLfloat> lines;
};

class RenderedText {
public:
    RenderedText(const Vec3& position, float scale, const Vec3& color, const string& text) :
            position(position), scale(scale/119), text(text) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
    }
    void draw() {
        glPushMatrix();
        glTranslated(position[0], position[1], position[2]);
        Vec4 rot = cameraTransform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale, scale, scale);
        glColor3fv(color);
        for (int i = 0; i < (int) text.size(); i++)
            glutStrokeCharacter(GLUT_STROKE_ROMAN, text[i]);
        glPopMatrix();
    }
    void computeBoundingSphere(Real& radius, Vec3& center) {
        center = position;
        radius = glutStrokeLength(GLUT_STROKE_ROMAN, (unsigned char*) text.c_str())*scale;
    }
private:
    Vec3 position;
    float scale;
    GLfloat color[3];
    string text;
};

class Scene {
public:
    vector<RenderedMesh> drawnMeshes;
    vector<RenderedMesh> solidMeshes;
    vector<RenderedMesh> transparentMeshes;
    vector<RenderedLine> lines;
    vector<RenderedText> strings;
};

class PendingCommand {
public:
    virtual void execute() = 0;
};

static int viewWidth, viewHeight;
static GLdouble fieldOfView = SimTK_PI/4;
static GLdouble nearClip = 1;
static GLdouble farClip = 1000;
static GLdouble groundHeight = 0;
static int groundAxis = 1;
static bool showGround = true, showShadows = true, showFPS = false;
static vector<PendingCommand*> pendingCommands;
static float fps = 0.0f;
static int fpsBaseTime = 0, fpsCounter = 0, nextMeshIndex;
static Scene* scene;
static string overlayMessage;
static int hideMessageTime = 0;

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

static void computeSceneBounds(Real& radius, Vec3& center) {
    // Record the bounding sphere of every object in the scene.

    vector<Vec3> centers;
    vector<Real> radii;
    for (int i = 0; i < (int) scene->drawnMeshes.size(); i++) {
        Vec3 center;
        Real radius;
        scene->drawnMeshes[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->solidMeshes.size(); i++) {
        Vec3 center;
        Real radius;
        scene->solidMeshes[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->transparentMeshes.size(); i++) {
        Vec3 center;
        Real radius;
        scene->transparentMeshes[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->lines.size(); i++) {
        Vec3 center;
        Real radius;
        scene->lines[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }
    for (int i = 0; i < (int) scene->strings.size(); i++) {
        Vec3 center;
        Real radius;
        scene->strings[i].computeBoundingSphere(radius, center);
        centers.push_back(center);
        radii.push_back(radius);
    }

    // Find the overall bounding sphere of the scene.

    if (centers.size() == 0) {
        radius = 0;
        center = Vec3(0);
    }
    else {
        Vec3 lower = centers[0]-radii[0];
        Vec3 upper = centers[0]+radii[0];
        for (int i = 1; i < (int) centers.size(); i++) {
            for (int j = 0; j < 3; j++) {
                lower[j] = min(lower[j], centers[i][j]-radii[i]);
                upper[j] = max(upper[j], centers[i][j]+radii[i]);
            }
        }
        center = 0.5*(lower+upper);
        radius = 0;
        for (int i = 0; i < (int) centers.size(); i++)
            radius = max(radius, (centers[i]-center).norm()+radii[i]);
    }
}

static void zoomCameraToShowWholeScene() {
    Real radius;
    Vec3 center;
    computeSceneBounds(radius, center);
    Real viewDistance = radius/tan(0.5*min(fieldOfView, fieldOfView*viewWidth/viewHeight));
    cameraTransform.updT() = center+cameraTransform.R()*Vec3(0, 0, viewDistance);
}

class PendingCameraZoom : public PendingCommand {
public:
    void execute() {
        zoomCameraToShowWholeScene();
    }
};

void menuSelected(int option) {
    char command = MENU_SELECTED;
    WRITE(outPipe, &command, 1);
    WRITE(outPipe, &option, sizeof(int));
}

class Menu {
public:
    Menu(string title, vector<pair<string, int> > items, void(*handler)(int)) : title(title), items(items), handler(handler), hasCreated(false) {
    }
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
        glRasterPos2f(x+2, y);
        for (int i = 0; i < (int) title.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, title[i]);
        glBegin(GL_TRIANGLES);
        glVertex2i(maxx-2, y-8);
        glVertex2i(maxx-10,y-8);
        glVertex2i(maxx-6, y-4);
        glEnd();
        return width+25;
    }
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

static void drawGroundAndSky(Real farClipDistance) {
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

    Real viewDistance = 0.9999*farClipDistance;
    Real xwidth = viewDistance*tan(0.5*fieldOfView*viewWidth/viewHeight);
    Real ywidth = viewDistance*tan(0.5*fieldOfView);
    Vec3 center = cameraTransform.T()-cameraTransform.R()*Vec3(0, 0, viewDistance);
    Vec3 corner1 = center+cameraTransform.R()*Vec3(-xwidth, -ywidth, 0);
    Vec3 corner2 = center+cameraTransform.R()*Vec3(xwidth, -ywidth, 0);
    Vec3 corner3 = center+cameraTransform.R()*Vec3(xwidth, ywidth, 0);
    Vec3 corner4 = center+cameraTransform.R()*Vec3(-xwidth, ywidth, 0);
    Vec3 cameraPosition, upDirection;
    if (groundAxis == 0) {
        cameraPosition = Vec3(groundHeight, cameraTransform.T()[1], cameraTransform.T()[2]);
        upDirection = Vec3(1, 0, 0);
    }
    else if (groundAxis == 1) {
        cameraPosition = Vec3(cameraTransform.T()[0], groundHeight, cameraTransform.T()[2]);
        upDirection = Vec3(0, 1, 0);
    }
    else {
        cameraPosition = Vec3(cameraTransform.T()[0], cameraTransform.T()[1], groundHeight);
        upDirection = Vec3(0, 0, 1);
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
    corner1 = center+Vec3(-farClipDistance, 0, -farClipDistance);
    corner2 = center+Vec3(farClipDistance, 0, -farClipDistance);
    corner3 = center+Vec3(farClipDistance, 0, farClipDistance);
    corner4 = center+Vec3(-farClipDistance, 0, farClipDistance);
    glUseProgram(groundProgram);
    Mat<4, 4, GLfloat> transform(1.0f);
    Vec3 sdir, tdir;
    if (groundAxis == 0) {
        transform[0][0] = transform[1][1] = 0.0f;
        transform[0][1] = transform[1][0] = 1.0f;
        sdir = Vec3(0, 1, 0);
        tdir = Vec3(0, 0, 1);
    }
    else if (groundAxis == 1) {
        sdir = Vec3(1, 0, 0);
        tdir = Vec3(0, 0, 1);
    }
    else {
        transform[1][1] = transform[2][2] = 0.0f;
        transform[1][2] = transform[2][1] = 1.0f;
        sdir = Vec3(1, 0, 0);
        tdir = Vec3(0, 1, 0);
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
        transform2[0][1] = 0.2;
        transform2[1][1] = 0.0;
        transform2[2][1] = 0.2;
        transform2 = transform*transform2;
        glPushMatrix();
        glMultMatrixf(&transform2[0][0]);
        glUniform3f(glGetUniformLocation(groundProgram, "color2"), 0.5f, 0.4f, 0.35f);
        glDepthRange(0.0, 0.9999);
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
            scene->solidMeshes[i].draw();
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
static pthread_mutex_t sceneLock;

static void changeSize(int width, int height) {
    if (height == 0)
        height = 1;
    viewWidth = width;
    viewHeight = height;
}

static void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glDisableClientState(GL_NORMAL_ARRAY);
    pthread_mutex_lock(&sceneLock);
    if (scene != NULL) {
        needRedisplay = false;

        // Execute any pending commands that need to be executed on the rendering thread.

        for (int i = 0; i < (int) pendingCommands.size(); i++) {
            pendingCommands[i]->execute();
            delete pendingCommands[i];
        }
        pendingCommands.clear();

        // Set up the viewpoint.

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glViewport(0, 0, viewWidth, viewHeight);
        Real sceneRadius;
        Vec3 sceneCenter;
        computeSceneBounds(sceneRadius, sceneCenter);
        Real centerDepth = ~(cameraTransform.T()-sceneCenter)*(cameraTransform.R().col(2));
        Real nearClipDistance, farClipDistance;
        if (showGround) {
            nearClipDistance = nearClip;
            farClipDistance = min(farClip, max(100.0, centerDepth+sceneRadius));
        }
        else {
            nearClipDistance = max(nearClip, centerDepth-sceneRadius);
            farClipDistance = min(farClip, centerDepth+sceneRadius);
        }
        gluPerspective(fieldOfView*SimTK_RADIAN_TO_DEGREE, (GLdouble) viewWidth/viewHeight, nearClipDistance, farClipDistance);
        glMatrixMode(GL_MODELVIEW);
        Vec3 cameraPos = cameraTransform.T();
        Vec3 centerPos = cameraTransform.T()+cameraTransform.R()*Vec3(0, 0, -1);
        Vec3 upDir = cameraTransform.R()*Vec3(0, 1, 0);
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
    }
    pthread_mutex_unlock(&sceneLock);
}

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
    glTranslatef(0, -viewHeight, 0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    int menux = 10;
    for (int i = 0; i < (int) menus.size(); i++)
        menux += menus[i].draw(menux, viewHeight-10);

    // Draw the frame rate counter.

    if (showFPS) {
        stringstream fpsstream;
        fpsstream << "Frames/Second: ";
        fpsstream << fps;
        string fps = fpsstream.str();
        glColor3f(1.0f, 0.5f, 0.0f);
        glRasterPos2f(10, 25);
        for (int i = 0; i < (int) fps.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, fps[i]);
    }

    // Draw a message overlay.

    if (hideMessageTime > 0) {
        glColor3f(1.0f, 0.5f, 0.0f);
        int width = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char*) overlayMessage.c_str());
        glRasterPos2f(std::max(0, (viewWidth-width)/2), viewHeight/2);
        for (int i = 0; i < (int) overlayMessage.size(); i++)
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, overlayMessage[i]);
    }
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glutSwapBuffers();
}

static void mousePressed(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        clickModifiers = glutGetModifiers();
        clickButton = button;
        clickX = x;
        clickY = y;
        if (clickButton == GLUT_LEFT_BUTTON) {
            Real radius;
            Vec3 sceneCenter;
            computeSceneBounds(radius, sceneCenter);
            Real distToCenter = (sceneCenter-cameraTransform.T()).norm();
            Real defaultDistance = radius;
            if (distToCenter > defaultDistance)
                rotateCenter = sceneCenter;
            else {
                Vec3 cameraDir = cameraTransform.R()*Vec3(0, 0, -1);
                Vec3 lookAt = cameraTransform.T()+defaultDistance*cameraDir;
                Real fract = (defaultDistance-distToCenter)/defaultDistance;
                rotateCenter = fract*lookAt+(1.0-fract)*sceneCenter;
            }
        }
        if (clickButton == 3 || clickButton == 4) {
            // Scroll wheel.

            cameraTransform.updT() += cameraTransform.R()*Vec3(0, 0, clickButton == 3 ? -0.5 : 0.5);
            pthread_mutex_lock(&sceneLock);
            needRedisplay = true;
            pthread_mutex_unlock(&sceneLock);
        }
    }
}

static void mouseDragged(int x, int y) {
    pthread_mutex_lock(&sceneLock);
    if (clickButton == GLUT_RIGHT_BUTTON || (clickButton == GLUT_LEFT_BUTTON && clickModifiers & GLUT_ACTIVE_SHIFT))
       cameraTransform.updT() += cameraTransform.R()*Vec3(0.01*(clickX-x), 0.01*(y-clickY), 0);
    else if (clickButton == GLUT_MIDDLE_BUTTON || (clickButton == GLUT_LEFT_BUTTON && clickModifiers & GLUT_ACTIVE_ALT))
       cameraTransform.updT() += cameraTransform.R()*Vec3(0, 0, 0.05*(clickY-y));
    else if (clickButton == GLUT_LEFT_BUTTON) {
        Vec3 cameraPos = cameraTransform.T();
        Vec3 cameraDir = cameraTransform.R()*Vec3(0, 0, -1);
        Vec3 upDir = cameraTransform.R()*Vec3(0, 1, 0);
        Rotation r;
        if (clickModifiers & GLUT_ACTIVE_CTRL)
            r.setRotationFromAngleAboutAxis(0.01*(clickY-y)-0.01*(clickX-x), ZAxis);
        else
            r.setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence, 0.01*(clickY-y), XAxis, 0.01*(clickX-x), YAxis);
        r = cameraTransform.R()*r*~cameraTransform.R();
        cameraPos = r*(cameraPos-rotateCenter)+rotateCenter;
        cameraDir = r*cameraDir;
        upDir = r*upDir;
        cameraTransform.updT() = cameraPos;
        cameraTransform.updR().setRotationFromTwoAxes(UnitVec3(-cameraDir), ZAxis, upDir, YAxis);
    }
    else
        return;
    clickX = x;
    clickY = y;
    needRedisplay = true;
    pthread_mutex_unlock(&sceneLock);
}

static void mouseMoved(int x, int y) {
    for (int i = 0; i < (int) menus.size(); i++)
        menus[i].mouseMoved(x, y);
}

static void keyPressed(unsigned char key, int x, int y) {
    char command = KEY_PRESSED;
    WRITE(outPipe, &command, 1);
    unsigned char buffer[2];
    buffer[0] = key;
    buffer[1] = 0;
    int modifiers = glutGetModifiers();
    if ((modifiers & GLUT_ACTIVE_SHIFT) != 0)
        buffer[1] += VisualizationEventListener::SHIFT_DOWN;
    if ((modifiers & GLUT_ACTIVE_CTRL) != 0)
        buffer[1] += VisualizationEventListener::CONTROL_DOWN;
    if ((modifiers & GLUT_ACTIVE_ALT) != 0)
        buffer[1] += VisualizationEventListener::ALT_DOWN;
    WRITE(outPipe, buffer, 2);
}

static void animateDisplay(int value) {
    if (hideMessageTime > 0 && hideMessageTime < glutGet(GLUT_ELAPSED_TIME)) {
        hideMessageTime = 0;
        needRedisplay = true;
    }
    if (needRedisplay)
        redrawDisplay();
    glutTimerFunc(33, animateDisplay, 0);
}

static void setOverlayMessage(const string& message) {
    overlayMessage = message;
    hideMessageTime = glutGet(GLUT_ELAPSED_TIME)+5000;
}

static void saveImage() {
    // Create offscreen buffers for rendering the image.

    GLuint frameBuffer, colorBuffer, depthBuffer;
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
    glGenRenderbuffers(1, &colorBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, colorBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGB8, viewWidth, viewHeight);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, colorBuffer);
    glGenRenderbuffers(1, &depthBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, viewWidth, viewHeight);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);

    // Render the image and load it into memory.

    renderScene();
    vector<unsigned char> data(viewWidth*viewHeight*3);
    glReadPixels(0, 0, viewWidth, viewHeight, GL_RGB, GL_UNSIGNED_BYTE, &data[0]);
    glDeleteRenderbuffers(1, &colorBuffer);
    glDeleteRenderbuffers(1, &depthBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDeleteFramebuffers(1, &frameBuffer);

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

void readData(char* buffer, int bytes) {
    int totalRead = 0;
    while (totalRead < bytes)
        totalRead += READ(inPipe, buffer+totalRead, bytes-totalRead);
}

void* listenForInput(void* args) {
    char buffer[256];
    float* floatBuffer = (float*) buffer;
    int* intBuffer = (int*) buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;
    while (true) {
        // Read a new scene.
        readData(buffer, 1);
        switch (buffer[0]) {
            case SET_CAMERA: {
                readData(buffer, 6*sizeof(float));
                pthread_mutex_lock(&sceneLock);
                cameraTransform.updR().setRotationToBodyFixedXYZ(Vec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
                cameraTransform.updT() = Vec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
                pthread_mutex_unlock(&sceneLock);
                break;
            }
            case ZOOM_CAMERA: {
                pthread_mutex_lock(&sceneLock);
                pendingCommands.push_back(new PendingCameraZoom());
                pthread_mutex_unlock(&sceneLock);
                break;
            }
            case LOOK_AT: {
                readData(buffer, 6*sizeof(float));
                Vec3 point(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
                Vec3 updir(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
                pthread_mutex_lock(&sceneLock);
                cameraTransform.updR().setRotationFromTwoAxes(UnitVec3(cameraTransform.T()-point), ZAxis, updir, YAxis);
                pthread_mutex_unlock(&sceneLock);
                break;
            }
            case SET_FIELD_OF_VIEW: {
                readData(buffer, sizeof(float));
                pthread_mutex_lock(&sceneLock);
                fieldOfView = floatBuffer[0];
                pthread_mutex_unlock(&sceneLock);
                break;
            }
            case SET_CLIP_PLANES: {
                readData(buffer, 2*sizeof(float));
                pthread_mutex_lock(&sceneLock);
                nearClip = floatBuffer[0];
                farClip = floatBuffer[1];
                pthread_mutex_unlock(&sceneLock);
                break;
            }
            case SET_GROUND_POSITION: {
                readData(buffer, sizeof(float)+sizeof(short));
                pthread_mutex_lock(&sceneLock);
                groundHeight = floatBuffer[0];
                groundAxis = shortBuffer[sizeof(float)/sizeof(short)];
                pthread_mutex_unlock(&sceneLock);
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
                pthread_mutex_lock(&sceneLock);
                menus.push_back(Menu(title, items, menuSelected));
                pthread_mutex_unlock(&sceneLock);
                break;
            }
            case START_OF_SCENE: {
                Scene* newScene = new Scene();
                bool finished = false;
                while (!finished) {
                    readData(buffer, 1);
                    char command = buffer[0];
                    switch (command) {
                        case END_OF_SCENE:
                            pthread_mutex_lock(&sceneLock);
                            if (scene != NULL)
                                delete scene;
                            scene = newScene;
                            needRedisplay = true;
                            fpsCounter++;
                            pthread_mutex_unlock(&sceneLock);
                            finished = true;
                            break;
                        case ADD_POINT_MESH:
                        case ADD_WIREFRAME_MESH:
                        case ADD_SOLID_MESH: {
                            readData(buffer, 13*sizeof(float)+sizeof(short));
                            Transform position;
                            position.updR().setRotationToBodyFixedXYZ(Vec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
                            position.updT() = Vec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
                            Vec3 scale = Vec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
                            Vec4 color = Vec4(floatBuffer[9], floatBuffer[10], floatBuffer[11], floatBuffer[12]);
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
                            Vec3 color = Vec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
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
                            Vec3 position = Vec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
                            float scale = floatBuffer[3];
                            Vec3 color = Vec3(floatBuffer[4], floatBuffer[5], floatBuffer[6]);
                            short length = shortBuffer[7*sizeof(float)/sizeof(short)];
                            readData(buffer, length);
                            newScene->strings.push_back(RenderedText(position, scale, color, string(buffer, length)));
                            break;
                        }
                        case ADD_FRAME: {
                            readData(buffer, 10*sizeof(float));
                            Rotation rotation;
                            rotation.setRotationToBodyFixedXYZ(Vec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
                            Vec3 position(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
                            float axisLength = floatBuffer[6];
                            float textScale = 0.2f*axisLength;
                            float lineThickness = 1;
                            Vec3 color = Vec3(floatBuffer[7], floatBuffer[8], floatBuffer[9]);
                            int index;
                            int numLines = newScene->lines.size();
                            for (index = 0; index < numLines && (color != newScene->lines[index].getColor() || newScene->lines[index].getThickness() != lineThickness); index++)
                                ;
                            if (index == numLines)
                                newScene->lines.push_back(RenderedLine(color, lineThickness));
                            vector<GLfloat>& line = newScene->lines[index].getLines();
                            Vec3 end = position+rotation*Vec3(axisLength, 0, 0);
                            line.push_back(position[0]);
                            line.push_back(position[1]);
                            line.push_back(position[2]);
                            line.push_back(end[0]);
                            line.push_back(end[1]);
                            line.push_back(end[2]);
                            newScene->strings.push_back(RenderedText(end, textScale, color, "X"));
                            end = position+rotation*Vec3(0, axisLength, 0);
                            line.push_back(position[0]);
                            line.push_back(position[1]);
                            line.push_back(position[2]);
                            line.push_back(end[0]);
                            line.push_back(end[1]);
                            line.push_back(end[2]);
                            newScene->strings.push_back(RenderedText(end, textScale, color, "Y"));
                            end = position+rotation*Vec3(0, 0, axisLength);
                            line.push_back(position[0]);
                            line.push_back(position[1]);
                            line.push_back(position[2]);
                            line.push_back(end[0]);
                            line.push_back(end[1]);
                            line.push_back(end[2]);
                            newScene->strings.push_back(RenderedText(end, textScale, color, "Z"));
                            break;
                        }
                        case DEFINE_MESH: {
                            readData(buffer, 2*sizeof(short));
                            PendingMesh* mesh = new PendingMesh();
                            int numVertices = shortBuffer[0];
                            int numFaces = shortBuffer[1];
                            mesh->vertices.resize(3*numVertices, 0);
                            mesh->normals.resize(3*numVertices);
                            mesh->faces.resize(3*numFaces);
                            readData((char*) &mesh->vertices[0], mesh->vertices.size()*sizeof(float));
                            readData((char*) &mesh->faces[0], mesh->faces.size()*sizeof(short));

                            // Compute normal vectors for the mesh.

                            vector<Vec3> normals(numVertices, Vec3(0));
                            for (int i = 0; i < numFaces; i++) {
                                int v1 = mesh->faces[3*i];
                                int v2 = mesh->faces[3*i+1];
                                int v3 = mesh->faces[3*i+2];
                                Vec3 vert1(mesh->vertices[3*v1], mesh->vertices[3*v1+1], mesh->vertices[3*v1+2]);
                                Vec3 vert2(mesh->vertices[3*v2], mesh->vertices[3*v2+1], mesh->vertices[3*v2+2]);
                                Vec3 vert3(mesh->vertices[3*v3], mesh->vertices[3*v3+1], mesh->vertices[3*v3+2]);
                                Vec3 norm = (vert2-vert1)%(vert3-vert1);
                                Real length = norm.norm();
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

                            // A real mesh will be generated from this the next time the screen is redrawn.

                            pthread_mutex_lock(&sceneLock);
                            pendingCommands.insert(pendingCommands.begin(), mesh);
                            pthread_mutex_unlock(&sceneLock);
                            break;
                        }
                        default:
                            SimTK_ASSERT_ALWAYS(false, "Unexpected data sent to visualizer");
                    }
                }
                break;
            }
            default:
                SimTK_ASSERT_ALWAYS(false, "Unexpected data sent to visualizer");
        }
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

void viewMenuSelected(int option) {
    Rotation groundRotation;
    if (groundAxis == 0)
        groundRotation.setRotationFromTwoAxes(UnitVec3(0, 1, 0), ZAxis, Vec3(1, 0, 0), YAxis);
    else if (groundAxis == 2)
        groundRotation.setRotationFromTwoAxes(UnitVec3(1, 0, 0), ZAxis, Vec3(0, 0, 1), YAxis);
    switch (option) {
        case MENU_VIEW_FRONT:
            cameraTransform.updR().setRotationToIdentityMatrix();
            cameraTransform.updR() = groundRotation*cameraTransform.R();
            zoomCameraToShowWholeScene();
            break;
        case MENU_VIEW_BACK:
            cameraTransform.updR().setRotationFromAngleAboutY(SimTK_PI);
            cameraTransform.updR() = groundRotation*cameraTransform.R();
            zoomCameraToShowWholeScene();
            break;
        case MENU_VIEW_LEFT:
            cameraTransform.updR().setRotationFromAngleAboutY(-0.5*SimTK_PI);
            cameraTransform.updR() = groundRotation*cameraTransform.R();
            zoomCameraToShowWholeScene();
            break;
        case MENU_VIEW_RIGHT:
            cameraTransform.updR().setRotationFromAngleAboutY(0.5*SimTK_PI);
            cameraTransform.updR() = groundRotation*cameraTransform.R();
            zoomCameraToShowWholeScene();
            break;
        case MENU_VIEW_TOP:
            cameraTransform.updR().setRotationFromAngleAboutX(-0.5*SimTK_PI);
            cameraTransform.updR() = groundRotation*cameraTransform.R();
            zoomCameraToShowWholeScene();
            break;
        case MENU_VIEW_BOTTOM:
            cameraTransform.updR().setRotationFromAngleAboutX(0.5*SimTK_PI);
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
    }
    needRedisplay = true;
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
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(500, 500);
    glutCreateWindow(title.c_str());
    glutDisplayFunc(redrawDisplay);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mousePressed);
    glutMotionFunc(mouseDragged);
    glutPassiveMotionFunc(mouseMoved);
    glutKeyboardFunc(keyPressed);
    glutTimerFunc(33, animateDisplay, 0);

    // On Windows, initialize function pointers for GL "extensions".
#ifdef _WIN32
    glGenBuffers    = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");
    glBindBuffer    = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");
    glBufferData    = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
    glCreateShader  = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
    glShaderSource  = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
    glAttachShader  = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
    glLinkProgram   = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
    glUseProgram    = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
    glUniform3f     = (PFNGLUNIFORM3FPROC) wglGetProcAddress("glUniform3f");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
    glGenFramebuffers    = (PFNGLGENFRAMEBUFFERSPROC) wglGetProcAddress("glGenFramebuffers");
    glGenRenderbuffers   = (PFNGLGENRENDERBUFFERSPROC) wglGetProcAddress("glGenRenderbuffers");
    glBindFramebuffer    = (PFNGLBINDFRAMEBUFFERPROC) wglGetProcAddress("glBindFramebuffer");
    glBindRenderbuffer   = (PFNGLBINDRENDERBUFFERPROC) wglGetProcAddress("glBindRenderbuffer");
    glRenderbufferStorage = (PFNGLRENDERBUFFERSTORAGEPROC) wglGetProcAddress("glRenderbufferStorage");
    glFramebufferRenderbuffer = (PFNGLFRAMEBUFFERRENDERBUFFERPROC) wglGetProcAddress("glFramebufferRenderbuffer");
    glDeleteRenderbuffers = (PFNGLDELETERENDERBUFFERSPROC) wglGetProcAddress("glDeleteRenderbuffers");
    glDeleteFramebuffers  = (PFNGLDELETEFRAMEBUFFERSPROC) wglGetProcAddress("glDeleteFramebuffers");
#endif

    // Set up lighting.

    GLfloat light_diffuse[] = {0.8, 0.8, 0.8, 1};
    GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
    GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1};
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
    menus.push_back(Menu("View", items, viewMenuSelected));

    // Spawn the listener thread.

    pthread_mutex_init(&sceneLock, NULL);
    pthread_t thread;
    pthread_create(&thread, NULL, listenForInput, NULL);

    // Enter the main loop.

    glutMainLoop();
    return 0;
}
