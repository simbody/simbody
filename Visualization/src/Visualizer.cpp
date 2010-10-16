#include "SimTKsimbody.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/VisualizationProtocol.h"
#include "simbody/internal/VisualizationEventListener.h"
#include <cstdlib>
#include <pthread.h>
#include <sstream>
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
    // Launch the GUI application.

    const char* GUI_APP_NAME = "VisualizationGUI";
    const String path = Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK")+"bin/"+GUI_APP_NAME;

    int pipes[2];
#ifdef _WIN32
    int status = _pipe(pipes, 16384, _O_BINARY);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    outPipe = pipes[1];
    String vizPipeToSim(pipes[0]); // convert pipe number to string
    status = _pipe(pipes, 16384, _O_BINARY);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    inPipe = pipes[0];
    String vizPipeFromSim(pipes[1]); // convert pipe number to string
    status = _spawnl(P_NOWAIT, path.c_str(), GUI_APP_NAME, vizPipeToSim.c_str(), vizPipeFromSim.c_str(), NULL);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to launch GUI");
#else
    pid_t pid;
    int status = pipe(pipes);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    outPipe = pipes[1];
    stringstream outPipeString, inPipeString;
    outPipeString << pipes[0];
    char outPipeArg[100], inPipeArg[100];
    outPipeString >> outPipeArg;
    status = pipe(pipes);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    inPipe = pipes[0];
    inPipeString << pipes[1];
    inPipeString >> inPipeArg;
    char* const argv[] = {GUI_APP_NAME, outPipeArg, inPipeArg, NULL};
    posix_spawn_file_actions_t fileActions;
    status = posix_spawn(&pid, path.c_str(), NULL, NULL, argv, environ);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to launch GUI");
#endif

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
