#include "SimTKsimbody.h"
#include "simbody/internal/Visualizer.h"
#include <spawn.h>
#include <cstdlib>
#include <string>
#include <unistd.h>

using namespace SimTK;
using namespace std;

static const char START_OF_SCENE = 0;
static const char END_OF_SCENE = 1;
static const char ADD_MESH = 2;

int outPipe;
#ifdef __APPLE__
#include <crt_externs.h>
char** environ = (*_NSGetEnviron());
#else
extern char** environ;
#endif

namespace SimTK {

Visualizer::Visualizer() {
    pid_t pid;
    char* GUI_APP_NAME = "VisualizationGUI";
    char* const argv[] = {GUI_APP_NAME, NULL};
    string path = Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK")+"bin/"+GUI_APP_NAME;
    int pipes[2];
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
    buffer[0] = rot[0];
    buffer[1] = rot[1];
    buffer[2] = rot[2];
    buffer[3] = transform.T()[0];
    buffer[4] = transform.T()[1];
    buffer[5] = transform.T()[2];
    buffer[6] = scale[0];
    buffer[7] = scale[1];
    buffer[8] = scale[2];
    buffer[9] = color[0];
    buffer[10] = color[1];
    buffer[11] = color[2];
    write(outPipe, buffer, 12*sizeof(float));
    write(outPipe, &meshIndex, sizeof(unsigned short));
}

}
