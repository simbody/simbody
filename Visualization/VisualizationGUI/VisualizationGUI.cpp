#include "SimTKcommon.h"
#include "simbody/internal/VisualizationProtocol.h"
#include "simbody/internal/VisualizationEventListener.h"
#ifdef _WIN32
    #include <windows.h>
    #include <io.h>
    #include <GL/gl.h>
    #include "glext.h"
    #include "glut.h"
    #define READ _read
    PFNGLGENBUFFERSPROC glGenBuffers;
    PFNGLBINDBUFFERPROC glBindBuffer;
    PFNGLBUFFERDATAPROC glBufferData;
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
#include <algorithm>
#include <pthread.h>
#include <set>
#include <vector>
#include <stdio.h>
#include <utility>
#include <string>

using namespace SimTK;
using namespace std;

static Transform cameraTransform(Rotation(), Vec3(0, 0, 10));
static int clickModifiers;
static int clickButton;
static int clickX;
static int clickY;
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
static GLdouble farClip = 100;
static vector<PendingCommand*> pendingCommands;
static Scene* scene;

class PendingMesh : public PendingCommand {
public:
    void execute() {
        meshes.push_back(new Mesh(vertices, normals, faces));
    }
    vector<float> vertices;
    vector<float> normals;
    vector<GLushort> faces;
};

class PendingCameraZoom : public PendingCommand {
public:
    void execute() {
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
        if (centers.size() > 0) {
            // Find the overall bounding sphere of the scene.

            Vec3 lower = centers[0]-radii[0];
            Vec3 upper = centers[0]+radii[0];
            for (int i = 1; i < (int) centers.size(); i++) {
                for (int j = 0; j < 3; j++) {
                    lower[j] = min(lower[j], centers[i][j]-radii[i]);
                    upper[j] = max(upper[j], centers[i][j]+radii[i]);
                }
            }
            Vec3 center = 0.5*(lower+upper);
            Real radius = 0;
            for (int i = 0; i < (int) centers.size(); i++)
                radius = max(radius, (centers[i]-center).norm()+radii[i]);
            Real viewDistance = radius/tan(0.5*min(fieldOfView, fieldOfView*viewWidth/viewHeight));
            cameraTransform.updT() = center+cameraTransform.R()*Vec3(0, 0, viewDistance);
        }
    }
};

void menuSelected(int option) {

}

class Menu {
public:
    Menu(string title, vector<pair<string, int> > items) : title(title), items(items), hasCreated(false) {
    }
    int draw(int x, int y) {
        if (!hasCreated) {
            id = glutCreateMenu(menuSelected);
            for (int i = 0; i < (int) items.size(); i++)
                glutAddMenuEntry(items[i].first.c_str(), items[i].second);
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
    int id, minx, miny, maxx, maxy;
    bool hasCreated;
    static int currentMenu;
};

int Menu::currentMenu = -1;

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
        gluPerspective(fieldOfView*SimTK_RADIAN_TO_DEGREE, (GLdouble) viewWidth/viewHeight, nearClip, farClip);
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
        glEnable(GL_BLEND);
        glDepthMask(GL_FALSE);
        vector<pair<float, int> > order(scene->transparentMeshes.size());
        for (int i = 0; i < (int) order.size(); i++)
            order[i] = make_pair((float)(~cameraTransform.R()*scene->transparentMeshes[i].getTransform().T())[2], i);
        sort(order.begin(), order.end());
        for (int i = 0; i < (int) order.size(); i++)
            scene->transparentMeshes[order[i].second].draw();
    }
    pthread_mutex_unlock(&sceneLock);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);

    // Draw menus.

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
        Vec3 centerPos = cameraTransform.T()+10*cameraDir;
        Vec3 upDir = cameraTransform.R()*Vec3(0, 1, 0);
        Rotation r;
        if (clickModifiers & GLUT_ACTIVE_CTRL)
            r.setRotationFromAngleAboutAxis(0.01*(clickY-y)-0.01*(clickX-x), ZAxis);
        else
            r.setRotationFromTwoAnglesTwoAxes(SpaceRotationSequence, 0.01*(clickY-y), XAxis, 0.01*(clickX-x), YAxis);
        r = cameraTransform.R()*r*~cameraTransform.R();
        cameraPos = r*(cameraPos-centerPos)+centerPos;
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
    write(outPipe, &command, 1);
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
    write(outPipe, buffer, 2);
}

static void animateDisplay(int value) {
    if (needRedisplay)
        renderScene();
    glutTimerFunc(33, animateDisplay, 0);
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

int main(int argc, char** argv) {
    stringstream(argv[1]) >> inPipe;
    stringstream(argv[2]) >> outPipe;

    // Initialize GLUT.

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(500, 500);
    glutCreateWindow("SimTK Visualizer");
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mousePressed);
    glutMotionFunc(mouseDragged);
    glutPassiveMotionFunc(mouseMoved);
    glutKeyboardFunc(keyPressed);
    glutTimerFunc(33, animateDisplay, 0);

    // On Windows, initialize function pointers for GL "extensions".
#ifdef _WIN32
    glGenBuffers = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");
    glBindBuffer = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");
    glBufferData = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
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
    makeBox();
    makeSphere();
    makeCylinder();
    makeCircle();
    scene = NULL;

    // Spawn the listener thread.

    pthread_mutex_init(&sceneLock, NULL);
    pthread_t thread;
    pthread_create(&thread, NULL, listenForInput, NULL);

    // Enter the main loop.

    glutMainLoop();
    return 0;
}
