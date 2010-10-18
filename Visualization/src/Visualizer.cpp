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
#else
    pid_t pid;
    const char* const argv[] = {appName, vizPipeToSim, vizPipeFromSim, NULL};
    posix_spawn_file_actions_t fileActions;
    status = posix_spawn(&pid, localPath), NULL, NULL, argv, environ);
    if (status == -1)
        status = posix_spawn(&pid, installPath, NULL, NULL, argv, environ);
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
    float* floatBuffer = (float*) buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;
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
    status = spawnViz(localPath.c_str(), installPath.c_str(),
                      GuiAppName, sim2vizPipe[0], viz2simPipe[1]);

    // status==-1 means we failed to spawn either executable.
    SimTK_ERRCHK2_ALWAYS(status != -1, "Visualizer::ctor()", 
        "Unable to spawn the Visualization GUI; tried '%s' and '%s'.",
        localPath.c_str(), installPath.c_str());

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

void Visualizer::drawBox(const Transform& transform, const Vec3& scale, const Vec3& color) const {
    drawMesh(transform, scale, color, 0);
}

void Visualizer::drawEllipsoid(const Transform& transform, const Vec3& scale, const Vec3& color) const {
    drawMesh(transform, scale, color, 1);
}

void Visualizer::drawMesh(const Transform& transform, const Vec3& scale, const Vec3& color, int meshIndex) const {
    char command = ADD_MESH;
    write(outPipe, &command, 1);
    float buffer[12];
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
    write(outPipe, buffer, 12*sizeof(float));
    write(outPipe, &meshIndex, sizeof(unsigned short));
}

}
