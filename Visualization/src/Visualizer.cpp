#include "SimTKsimbody.h"
#include "simbody/internal/Visualizer.h"
#include <cstdlib>
#include <sstream>
#include <string>

using namespace SimTK;
using namespace std;

static const char START_OF_SCENE = 0;
static const char END_OF_SCENE = 1;
static const char ADD_MESH = 2;

#ifdef _WIN32
    #include <fcntl.h>
    #include <io.h>
    #include <process.h>
#else
    #include <spawn.h>
    #include <unistd.h>
    #ifdef __APPLE__
        #include <crt_externs.h>
        static char** environ = (*_NSGetEnviron());
    #else
        extern char** environ;
    #endif
#endif

namespace SimTK {

Visualizer::Visualizer() {
    char* GUI_APP_NAME = "VisualizationGUI";
    int pipes[2];
#ifdef _WIN32
    int status = _pipe(pipes, 16384, _O_BINARY);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    outPipe = pipes[1];
    string path = Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK")+"bin\\"+GUI_APP_NAME;
    stringstream pipeString;
    pipeString << pipes[0];
    status = _spawnl(P_NOWAIT, path.c_str(), GUI_APP_NAME, pipeString.str().c_str(), NULL);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to launch GUI");
#else
    pid_t pid;
    char* const argv[] = {GUI_APP_NAME, NULL};
    string path = Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK")+"bin/"+GUI_APP_NAME;
    int status = pipe(pipes);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to open pipe");
    outPipe = pipes[1];
    posix_spawn_file_actions_t fileActions;
    status = posix_spawn_file_actions_init(&fileActions);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to launch GUI");
    status = posix_spawn_file_actions_adddup2(&fileActions, pipes[0], 0);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to launch GUI");
    status = posix_spawn(&pid, path.c_str(), &fileActions, NULL, argv, environ);
    SimTK_ASSERT_ALWAYS(status != -1, "Visualizer: Failed to launch GUI");
#endif
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
