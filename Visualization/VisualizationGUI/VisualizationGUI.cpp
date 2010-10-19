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

using namespace SimTK;
using namespace std;

static Transform cameraTransform(Rotation(), Vec3(0, 0, 10));
static int clickModifiers;
static int clickButton;
static int clickX;
static int clickY;
static int inPipe, outPipe;
static bool needRedisplay;

class Mesh {
public:
    Mesh(vector<float>& vertices, vector<float>& normals, vector<GLushort>& faces) : numVertices(vertices.size()), faces(faces) {
        // Build OpenGL buffers.

        GLuint buffers[2];
        glGenBuffers(3, buffers);
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
    }
    void draw(short representation) const {
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glNormalPointer(GL_FLOAT, 0, 0);
        if (representation == DecorativeGeometry::DrawSurface)
            glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_SHORT, &faces[0]);
        else if (representation == DecorativeGeometry::DrawPoints)
            glDrawArrays(GL_POINTS, 0, numVertices);
        else if (representation == DecorativeGeometry::DrawWireframe)
            glDrawElements(GL_LINES, edges.size(), GL_UNSIGNED_SHORT, &edges[0]);
    }
private:
    int numVertices;
    GLuint vertBuffer, normBuffer;
    vector<GLushort> edges, faces;
};

class RenderedMesh {
public:
    RenderedMesh(const Transform& transform, const Vec3& scale, const Vec4& color, short representation, const Mesh& mesh) :
            transform(transform), scale(scale), representation(representation), mesh(&mesh) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
        this->color[3] = color[4];
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
        mesh->draw(representation);
        glPopMatrix();
    }
    const Transform& getTransform() const {
        return transform;
    }
private:
    Transform transform;
    Vec3 scale;
    GLfloat color[4];
    short representation;
    const Mesh* mesh;
};

class Scene {
public:
    vector<RenderedMesh> drawnMeshes;
    vector<RenderedMesh> solidMeshes;
    vector<RenderedMesh> transparentMeshes;
};

static vector<Mesh*> meshes;
static Scene* scene;
static pthread_mutex_t sceneLock;

static void changeSize(int width, int height) {
    if (height == 0)
        height = 1;
    float ratio = 1.0*width/height;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, width, height);
    gluPerspective(45,ratio,1,100);
}

static void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vec3 cameraPos = cameraTransform.T();
    Vec3 centerPos = cameraTransform.T()+cameraTransform.R()*Vec3(0, 0, -1);
    Vec3 upDir = cameraTransform.R()*Vec3(0, 1, 0);
    gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2], centerPos[0], centerPos[1], centerPos[2], upDir[0], upDir[1], upDir[2]);
    glClear(GL_COLOR_BUFFER_BIT);
    pthread_mutex_lock(&sceneLock);
    needRedisplay = false;
    glDisable(GL_LIGHTING);
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
    pthread_mutex_unlock(&sceneLock);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glutSwapBuffers();
}

static void mousePressed(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        clickModifiers = glutGetModifiers();
        clickButton = button;
        clickX = x;
        clickY = y;
    }
}

static void mouseDragged(int x, int y) {
    if (clickButton == GLUT_LEFT_BUTTON) {
        if (clickModifiers & GLUT_ACTIVE_SHIFT)
           cameraTransform.updT() += cameraTransform.R()*Vec3(0.01*(clickX-x), 0.01*(y-clickY), 0);
        else {
            Vec3 cameraPos = cameraTransform.T();
            Vec3 cameraDir = cameraTransform.R()*Vec3(0, 0, -1);
            Vec3 centerPos = cameraTransform.T()+10*cameraDir;
            Vec3 upDir = cameraTransform.R()*Vec3(0, 1, 0);
            Rotation r(SpaceRotationSequence, 0.01*(clickY-y), XAxis, 0.01*(clickX-x), YAxis);
            r = cameraTransform.R()*r*~cameraTransform.R();
            cameraPos = r*(cameraPos-centerPos)+centerPos;
            cameraDir = r*cameraDir;
            upDir = r*upDir;
            cameraTransform.updT() = cameraPos;
            cameraTransform.updR().setRotationFromTwoAxes(UnitVec3(-cameraDir), ZAxis, upDir, YAxis);
        }
        clickX = x;
        clickY = y;
        pthread_mutex_lock(&sceneLock);
        needRedisplay = true;
        pthread_mutex_unlock(&sceneLock);
    }
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
    const int numLatitude = 10;
    const int numLongitude = 15;
    const float radius = 1.0f;
    vector<GLfloat> vertices;
    vector<GLfloat> normals;
    vector<GLushort> faces;
    addVec(vertices, 0, radius, 0);
    addVec(normals, 0, 1, 0);
    for (int i = 0; i < numLatitude; i++)
    {
      float phi = (float) (((i+1)*SimTK_PI)/(numLatitude+1));
      float sphi = sin(phi);
      float cphi = cos(phi);
      float y = radius*cphi;
      float r = radius*sphi;
      for (int j = 0; j < numLongitude; j++)
      {
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
    for (int i = 1; i < numLatitude; i++)
    {
      int base = (i-1)*numLongitude+1;
      for (int j = 0; j < numLongitude; j++)
      {
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
        SimTK_ASSERT_ALWAYS(buffer[0] == START_OF_SCENE, "Unexpected data sent to visualizer");
        Scene* newScene = new Scene();
        bool finished = false;
        while (!finished) {
            readData(buffer, 1);
            char command = buffer[0];
            switch (command) {
                case END_OF_SCENE:
                    pthread_mutex_lock(&sceneLock);
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
                    RenderedMesh mesh(position, scale, color, representation, *meshes[shortBuffer[13*sizeof(float)/sizeof(short)]]);
                    if (command != ADD_SOLID_MESH)
                        newScene->drawnMeshes.push_back(mesh);
                    else if (color[3] == 1)
                        newScene->solidMeshes.push_back(mesh);
                    else
                        newScene->transparentMeshes.push_back(mesh);
                    break;
                }
                default:
                    SimTK_ASSERT_ALWAYS(false, "Unexpected data sent to visualizer");
            }
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
    glutInitWindowSize(320,320);
    glutCreateWindow("SimTK Visualizer");
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mousePressed);
    glutMotionFunc(mouseDragged);
    glutKeyboardFunc(keyPressed);
    glutTimerFunc(33, animateDisplay, 0);

    // On Windows, initialize function pointers for GL "extensions".
#ifdef _WIN32
    glGenBuffers = (PFNGLGENBUFFERSPROC) wglGetProcAddress("glGenBuffers");
    glBindBuffer = (PFNGLBINDBUFFERPROC) wglGetProcAddress("glBindBuffer");
    glBufferData = (PFNGLBUFFERDATAPROC) wglGetProcAddress("glBufferData");
#endif

    // Set up lighting.

    GLfloat light_diffuse[] = {0.5, 0.5, 0.5, 1};
    GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
    GLfloat light_ambient[] = {0.3, 0.3, 0.3, 1};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambient);
    glEnable(GL_LIGHT0);

    // Initialize miscellaneous OpenGL state.

    glEnable(GL_DEPTH_TEST);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_CULL_FACE);
    makeBox();
    makeSphere();
    scene = new Scene();

    // Spawn the listener thread.

    pthread_mutex_init(&sceneLock, NULL);
    pthread_t thread;
    pthread_create(&thread, NULL, listenForInput, NULL);

    // Enter the main loop.

    glutMainLoop();
    return 0;
}
