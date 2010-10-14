#include "SimTKcommon.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include <stdio.h>
using namespace SimTK;
using namespace std;

class Mesh {
public:
    Mesh(vector<float>& vertices, vector<float>& normals, vector<GLushort>& faces) : faces(faces) {
        GLuint buffers[2];
        glGenBuffers(2, buffers);
        vertBuffer = buffers[0];
        normBuffer = buffers[1];
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), &vertices[0], GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glBufferData(GL_ARRAY_BUFFER, normals.size()*sizeof(float), &normals[0], GL_STATIC_DRAW);
    }
    void draw() const {
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glNormalPointer(GL_FLOAT, 0, 0);
        glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_SHORT, &faces[0]);
    }
private:
    GLuint vertBuffer, normBuffer;
    vector<GLushort> faces;
};

class RenderedMesh {
public:
    RenderedMesh(const Transform& transform, const Vec3& scale, const Vec3& color, const Mesh& mesh) :
            transform(transform), scale(scale), mesh(&mesh) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
        this->color[3] = 1;
    }
    void draw() {
        glPushMatrix();
        glTranslated(transform.T()[0], transform.T()[1], transform.T()[2]);
        Vec4 rot = transform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, color);
        mesh->draw();
        glPopMatrix();
    }
private:
    Transform transform;
    Vec3 scale;
    GLfloat color[4];
    const Mesh* mesh;
};

static vector<Mesh*> meshes;
static vector<RenderedMesh> scene;
static pthread_mutex_t sceneLock;

static Vec3 cameraPos(0, 0, 30);
static Vec3 centerPos(0, 0, 0);
static Vec3 upDir(0, 1, 0);
static int clickModifiers;
static int clickButton;
static int clickX;
static int clickY;
static int pipes[2];
static bool needRedisplay;

static const char START_OF_SCENE = 0;
static const char END_OF_SCENE = 1;
static const char ADD_MESH = 2;

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
    gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2], centerPos[0], centerPos[1], centerPos[2], upDir[0], upDir[1], upDir[2]);
    glClear(GL_COLOR_BUFFER_BIT);
    pthread_mutex_lock(&sceneLock);
    needRedisplay = false;
    for (int i = 0; i < (int) scene.size(); i++)
        scene[i].draw();
    pthread_mutex_unlock(&sceneLock);
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
        Rotation r(SpaceRotationSequence, 0.01*(clickY-y), XAxis, 0.01*(clickX-x), YAxis);
        cameraPos = r*cameraPos;
        upDir = r*upDir;
        clickX = x;
        clickY = y;
        pthread_mutex_lock(&sceneLock);
        needRedisplay = true;
        pthread_mutex_unlock(&sceneLock);
    }
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
        totalRead += fread(buffer+totalRead, 1, bytes-totalRead, stdin);
}

void* listenForInput(void* args) {
    char buffer[256];
    float* floatBuffer = (float*) buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;
    while (true) {
        // Read a new scene.
        readData(buffer, 1);
        SimTK_ASSERT_ALWAYS(buffer[0] == START_OF_SCENE, "Unexpected data sent to visualizer");
        vector<RenderedMesh> newScene;
        bool finished = false;
        while (!finished) {
            readData(buffer, 1);
            Transform position;
            Vec3 scale, color;
            switch (buffer[0]) {
                case END_OF_SCENE:
                    pthread_mutex_lock(&sceneLock);
                    scene = newScene;
                    needRedisplay = true;
                    pthread_mutex_unlock(&sceneLock);
                    finished = true;
                    break;
                case ADD_MESH:
                    readData(buffer, 12*sizeof(float)+sizeof(unsigned short));
                    position.updR().setRotationToBodyFixedXYZ(Vec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
                    position.updT() = Vec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
                    scale = Vec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
                    color = Vec3(floatBuffer[9], floatBuffer[10], floatBuffer[11]);
                    newScene.push_back(RenderedMesh(position, scale, color, *meshes[shortBuffer[12*sizeof(float)/sizeof(unsigned short)]]));
                    break;
                default:
                    SimTK_ASSERT_ALWAYS(false, "Unexpected data sent to visualizer");
            }
        }
    }
    return 0;
}

int main(int argc, char* argv) {
    // Initialize GLUT.

    glutInit(&argc, &argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(320,320);
    glutCreateWindow("SimTK Visualizer");
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mousePressed);
    glutMotionFunc(mouseDragged);
    glutTimerFunc(33, animateDisplay, 0);

    // Set up lighting.

    GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Initialize miscellaneous OpenGL state.

    glEnable(GL_DEPTH_TEST);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    makeBox();
    makeSphere();

    // Spawn the listener thread.

    pthread_mutex_init(&sceneLock, NULL);
    pthread_t thread;
    pthread_create(&thread, NULL, listenForInput, NULL);

    // Enter the main loop.

    glutMainLoop();
    return 0;
}
