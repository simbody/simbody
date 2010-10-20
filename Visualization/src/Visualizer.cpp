#include "SimTKsimbody.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/VisualizationProtocol.h"
#include "simbody/internal/VisualizationEventListener.h"
#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <string>

using namespace SimTK;
using namespace std;

#ifdef _WIN32
    #include <fcntl.h>
    #include <io.h>
    #include <process.h>
    #define READ _read
#else
    #include <spawn.h>
    #include <unistd.h>
    #ifdef __APPLE__
        #include <crt_externs.h>
        static char** environ = (*_NSGetEnviron());
    #else
        extern char** environ;
    #endif
    #define READ read
#endif

static int inPipe;

// Create a pipe, using the right call for this platform.
static int createPipe(int pipeHandles[2]) {
    const int status =
#ifdef _WIN32
        _pipe(pipeHandles, 16384, _O_BINARY);
#else
        pipe(pipeHandles);
#endif
    return status;
}

// Spawn the visualizer GUI executable, using the right method for
// this platform. We take two executables to try in order,
// and return as soon as one of them fails. If neither works, we return
// status -1, otherwise we return the status code from the successful
// spawn.
static int spawnViz(const char* localPath, const char* installPath,
                    const char* appName, int toSimPipe, int fromSimPipe)
{
    int status;
    char vizPipeToSim[32], vizPipeFromSim[32];
    sprintf(vizPipeToSim, "%d", toSimPipe);
    sprintf(vizPipeFromSim, "%d", fromSimPipe);
#ifdef _WIN32
    status = _spawnl(P_NOWAIT, localPath, appName, vizPipeToSim, vizPipeFromSim, NULL);
    if (status == -1)
        status = _spawnl(P_NOWAIT, installPath, appName, vizPipeToSim, vizPipeFromSim, NULL);
    SimTK_ERRCHK2_ALWAYS(status != -1, "Visualizer::ctor()",
        "Unable to spawn the Visualization GUI; tried '%s' and '%s'.", localPath, installPath);
#else
    pid_t pid;
    const char* const argv[] = {appName, vizPipeToSim, vizPipeFromSim, NULL};
    posix_spawn_file_actions_t fileActions;
    status = posix_spawn(&pid, localPath, NULL, NULL, 
                         (char* const*)argv, environ);
    if (status != 0)
        status = posix_spawn(&pid, installPath, NULL, NULL, 
                             (char* const*)argv, environ);
    SimTK_ERRCHK2_ALWAYS(status == 0, "Visualizer::ctor()",
        "Unable to spawn the Visualization GUI; tried '%s' and '%s'.", localPath, installPath);
#endif
    return status;
}

static void readData(char* buffer, int bytes) {
    int totalRead = 0;
    while (totalRead < bytes)
        totalRead += READ(inPipe, buffer+totalRead, bytes-totalRead);
}

static void* listenForVisualizationEvents(void* arg) {
    Visualizer& visualizer = *reinterpret_cast<Visualizer*>(arg);
    const vector<VisualizationEventListener*>& listeners = visualizer.getEventListeners();
    char buffer[256];
    while (true) {
        // Receive an event.

        readData(buffer, 1);
        switch (buffer[0]) {
            case KEY_PRESSED:
                readData(buffer, 2);
                for (int i = 0; i < (int) listeners.size(); i++)
                    listeners[i]->keyPressed(buffer[0], buffer[1]);
                break;
            default:
                SimTK_ASSERT_ALWAYS(false, "Unexpected data received from visualizer");
        }
    }
    return 0;
}

namespace SimTK {

Visualizer::Visualizer() {
    // Launch the GUI application. We'll first look for one in the same directory
    // as the running executable; then if that doesn't work we'll look in the 
    // bin subdirectory of the SimTK installation.

    const char* GuiAppName = "VisualizationGUI";
    const String localPath = Pathname::getThisExecutableDirectory() + GuiAppName;
    const String installPath = 
        Pathname::addDirectoryOffset(Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK"), 
                                     "bin") + GuiAppName;

    int sim2vizPipe[2], viz2simPipe[2], status;

    // Create pipe pair for communication from simulator to visualizer.
    status = createPipe(sim2vizPipe);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    outPipe = sim2vizPipe[1];

    // Create pipe pair for communication from visualizer to simulator.
    status = createPipe(viz2simPipe);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    inPipe = viz2simPipe[0];

    // Spawn the visualizer gui, trying local first then installed version.
    spawnViz(localPath.c_str(), installPath.c_str(),
                      GuiAppName, sim2vizPipe[0], viz2simPipe[1]);

    // Spawn the thread to listen for events.

    pthread_t thread;
    pthread_create(&thread, NULL, listenForVisualizationEvents, this);
}

void Visualizer::addEventListener(VisualizationEventListener* listener) {
    listeners.push_back(listener);
}

const vector<VisualizationEventListener*>& Visualizer::getEventListeners() const {
    return listeners;
}

void Visualizer::beginScene() const {
    char command = START_OF_SCENE;
    write(outPipe, &command, 1);
}

void Visualizer::finishScene() const {
    char command = END_OF_SCENE;
    write(outPipe, &command, 1);
}

void Visualizer::drawBox(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const {
    drawMesh(transform, scale, color, (short) representation, 0);
}

void Visualizer::drawEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const {
    drawMesh(transform, scale, color, (short) representation, 1);
}

void Visualizer::drawCylinder(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const {
    drawMesh(transform, scale, color, (short) representation, 2);
}

void Visualizer::drawCircle(const Transform& transform, const Vec3& scale, const Vec4& color, int representation) const {
    drawMesh(transform, scale, color, (short) representation, 3);
}

void Visualizer::drawPolygonalMesh(const PolygonalMesh& mesh, const Transform& transform, Real scale, const Vec4& color, int representation) const {
    const void* impl = &mesh.getImpl();
    int index;
    map<const void*, int>::const_iterator iter = meshes.find(impl);
    if (iter == meshes.end()) {
        // This is a new mesh, so we need to send it to the visualizer.  Build lists of vertices and faces,
        // triangulating as necessary.

        vector<float> vertices;
        vector<unsigned short> faces;
        for (int i = 0; i < mesh.getNumVertices(); i++) {
            Vec3 pos = mesh.getVertexPosition(i);
            vertices.push_back((float) pos[0]);
            vertices.push_back((float) pos[1]);
            vertices.push_back((float) pos[2]);
        }
        for (int i = 0; i < mesh.getNumFaces(); i++) {
            int numVert = mesh.getNumVerticesForFace(i);
            if (numVert < 3)
                continue; // Ignore it.
            if (numVert == 3) {
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 1));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 2));
            }
            else if (numVert == 4) {
                // Split it into two triangles.

                faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 1));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 2));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 2));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 3));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
            }
            else {
                // Add a vertex at the center, then split it into triangles.

                Vec3 center(0);
                for (int j = 0; j < numVert; j++)
                    center += vertices[mesh.getFaceVertex(i, j)];
                center /= numVert;
                vertices.push_back((float) center[0]);
                vertices.push_back((float) center[1]);
                vertices.push_back((float) center[2]);
                int newIndex = vertices.size()-1;
                for (int j = 0; j < numVert-1; j++) {
                    faces.push_back((unsigned short) mesh.getFaceVertex(i, j));
                    faces.push_back((unsigned short) mesh.getFaceVertex(i, j+1));
                    faces.push_back((unsigned short) newIndex);
                }
            }
        }
        SimTK_ASSERT_ALWAYS(vertices.size() <= 65536*3, "DecorativeMesh cannot have more than 65536 vertices");
        SimTK_ASSERT_ALWAYS(faces.size() <= 65536*3, "DecorativeMesh cannot have more than 65536 faces");
        index = meshes.size()+4;
        meshes[impl] = index;
        char command = DEFINE_MESH;
        write(outPipe, &command, 1);
        unsigned short numVertices = vertices.size()/3;
        unsigned short numFaces = faces.size()/3;
        write(outPipe, &numVertices, sizeof(short));
        write(outPipe, &numFaces, sizeof(short));
        write(outPipe, &vertices[0], vertices.size()*sizeof(float));
        write(outPipe, &faces[0], faces.size()*sizeof(short));
    }
    else
        index = iter->second;
    drawMesh(transform, Vec3(scale), color, (short) representation, index);
}

void Visualizer::drawMesh(const Transform& transform, const Vec3& scale, const Vec4& color, short representation, short meshIndex) const {
    char command = (representation == DecorativeGeometry::DrawPoints ? ADD_POINT_MESH : (representation == DecorativeGeometry::DrawWireframe ? ADD_WIREFRAME_MESH : ADD_SOLID_MESH));
    write(outPipe, &command, 1);
    float buffer[13];
    Vec3 rot = transform.R().convertRotationToBodyFixedXYZ();
    buffer[0] = (float) rot[0];
    buffer[1] = (float) rot[1];
    buffer[2] = (float) rot[2];
    buffer[3] = (float) transform.T()[0];
    buffer[4] = (float) transform.T()[1];
    buffer[5] = (float) transform.T()[2];
    buffer[6] = (float) scale[0];
    buffer[7] = (float) scale[1];
    buffer[8] = (float) scale[2];
    buffer[9] = (float) color[0];
    buffer[10] = (float) color[1];
    buffer[11] = (float) color[2];
    buffer[12] = (float) color[3];
    write(outPipe, buffer, 13*sizeof(float));
    write(outPipe, &meshIndex, sizeof(short));
}

void Visualizer::drawLine(const Vec3& end1, const Vec3& end2, const Vec4& color, Real thickness) const {
    char command = ADD_LINE;
    write(outPipe, &command, 1);
    float buffer[10];
    buffer[0] = (float) color[0];
    buffer[1] = (float) color[1];
    buffer[2] = (float) color[2];
    buffer[3] = (float) thickness;
    buffer[4] = (float) end1[0];
    buffer[5] = (float) end1[1];
    buffer[6] = (float) end1[2];
    buffer[7] = (float) end2[0];
    buffer[8] = (float) end2[1];
    buffer[9] = (float) end2[2];
    write(outPipe, buffer, 10*sizeof(float));
}

void Visualizer::drawText(const Vec3& position, Real scale, const Vec4& color, const string& string) const {
    SimTK_ASSERT_ALWAYS(string.size() <= 256, "DecorativeText cannot be longer than 256 characters");
    char command = ADD_TEXT;
    write(outPipe, &command, 1);
    float buffer[7];
    buffer[0] = (float) position[0];
    buffer[1] = (float) position[1];
    buffer[2] = (float) position[2];
    buffer[3] = (float) scale;
    buffer[4] = (float) color[0];
    buffer[5] = (float) color[1];
    buffer[6] = (float) color[2];
    write(outPipe, buffer, 7*sizeof(float));
    short length = string.size();
    write(outPipe, &length, sizeof(short));
    write(outPipe, &string[0], length);
}

void Visualizer::drawFrame(const Transform& transform, Real axisLength, const Vec4& color) const {
    char command = ADD_FRAME;
    write(outPipe, &command, 1);
    float buffer[10];
    Vec3 rot = transform.R().convertRotationToBodyFixedXYZ();
    buffer[0] = (float) rot[0];
    buffer[1] = (float) rot[1];
    buffer[2] = (float) rot[2];
    buffer[3] = (float) transform.T()[0];
    buffer[4] = (float) transform.T()[1];
    buffer[5] = (float) transform.T()[2];
    buffer[6] = (float) axisLength;
    buffer[7] = (float) color[0];
    buffer[8] = (float) color[1];
    buffer[9] = (float) color[2];
    write(outPipe, buffer, 10*sizeof(float));
}

}
