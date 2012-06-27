/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "simbody/internal/common.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/Visualizer_InputListener.h"
#include "VisualizerProtocol.h"

#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <string>

using namespace SimTK;
using namespace std;

#ifdef _WIN32
    #include <fcntl.h>
    #include <io.h>
    #include <process.h>
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
    SimTK_ERRCHK4_ALWAYS(status!=-1, "VisualizerProtocol",  \
    "An attempt to write() %d bytes to pipe %d failed with errno=%d (%s).", \
    (len),(pipeno),errno,strerror(errno));}

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
// and return after the first one succeeds. If neither works, we throw
// an error that is hopefully helful.
static void spawnViz(const Array_<String>& searchPath,
                     const String& appName, int toSimPipe, int fromSimPipe)
{
    int status;
    char vizPipeToSim[32], vizPipeFromSim[32];
    sprintf(vizPipeToSim, "%d", toSimPipe);
    sprintf(vizPipeFromSim, "%d", fromSimPipe);

    String exePath; // search path + appName

#ifdef _WIN32
    intptr_t handle;
    for (unsigned i=0; i < searchPath.size(); ++i) {
        exePath = searchPath[i] + appName;
        handle = _spawnl(P_NOWAIT, exePath.c_str(), appName.c_str(), 
                         vizPipeToSim, vizPipeFromSim, (const char*)0);
        if (handle != -1)
            break; // success!
    }
    status = (handle==-1) ? -1 : 0;
#else
    const pid_t pid = fork();
    if (pid == 0) {
        // child process
        for (unsigned i=0; i < searchPath.size(); ++i) {
            exePath = searchPath[i] + appName;
            status = execl(exePath.c_str(), appName.c_str(), 
                           vizPipeToSim, vizPipeFromSim, (const char*)0); 
            // if we get here the execl() failed
        }
        // fall through -- we failed on every try
    } else {
        // parent process
        status = (pid==-1) ? -1 : 0;
    }
#endif

    if (status != 0) {
        // Create a chronicle of our failures above.
        String failedPath;
        for (unsigned i=0; i < searchPath.size(); ++i)
            failedPath += ("  " + searchPath[i] + "\n");

        SimTK_ERRCHK4_ALWAYS(status == 0, "VisualizerProtocol::ctor()",
            "Unable to spawn executable '%s' from directories:\n%s"
            "Final system error was errno=%d (%s).", 
            appName.c_str(), failedPath.c_str(), errno, strerror(errno));
    }
}

static void readDataFromPipe(int srcPipe, unsigned char* buffer, int bytes) {
    int totalRead = 0;
    while (totalRead < bytes)
        totalRead += READ(srcPipe, buffer+totalRead, bytes-totalRead);
}

static void readData(unsigned char* buffer, int bytes) 
{
    readDataFromPipe(inPipe, buffer, bytes);
}

static void* listenForVisualizerEvents(void* arg) {
    Visualizer& visualizer = *reinterpret_cast<Visualizer*>(arg);
    unsigned char buffer[256];

    try
  { while (true) {
        // Receive a user input event.

        readData(buffer, 1);
        switch (buffer[0]) {
        case KeyPressed: {
            readData(buffer, 2);
            const Array_<Visualizer::InputListener*>& listeners = visualizer.getInputListeners();
            unsigned keyCode = buffer[0];
            if (buffer[1] & Visualizer::InputListener::IsSpecialKey)
                keyCode += Visualizer::InputListener::SpecialKeyOffset;
            for (int i = 0; i < (int) listeners.size(); i++)
                if (listeners[i]->keyPressed(keyCode, (unsigned)(buffer[1])))
                    break; // key press has been handled
            break;
        }
        case MenuSelected: {
            int menu, item;
            readData((unsigned char*) &menu, sizeof(int));
            readData((unsigned char*) &item, sizeof(int));
            const Array_<Visualizer::InputListener*>& listeners = visualizer.getInputListeners();
            for (int i = 0; i < (int) listeners.size(); i++)
                if (listeners[i]->menuSelected(menu, item))
                    break; // menu event has been handled
            break;
        }
        case SliderMoved: {
            int slider;
            readData((unsigned char*) &slider, sizeof(int));
            float value;
            readData((unsigned char*) &value, sizeof(float));
            const Array_<Visualizer::InputListener*>& listeners = visualizer.getInputListeners();
            for (int i = 0; i < (int) listeners.size(); i++)
                if (listeners[i]->sliderMoved(slider, value))
                    break; // slider event has been handled
            break;
        }
        default:
            SimTK_ERRCHK1_ALWAYS(false, "listenForVisualizerEvents()",
                "Unexpected command %u received from VisualizerGUI. Can't continue.",
                (unsigned)buffer[0]);
        }
    }
  } catch (const std::exception& e) {
        std::cout << "Visualizer listenerThread: unrecoverable error:\n";
        std::cout << e.what() << std::endl;
        return (void*)1;
    }
    return (void*)0;
}

VisualizerProtocol::VisualizerProtocol
   (Visualizer& visualizer, const Array_<String>& userSearchPath) 
{
    // Launch the GUI application. We'll first look for one in the same directory
    // as the running executable; then if that doesn't work we'll look in the
    // bin subdirectory of the SimTK installation.

    const char* GuiAppName = "VisualizerGUI";

    Array_<String> actualSearchPath;
    // Always start with the current executable's directory.
    actualSearchPath.push_back(Pathname::getThisExecutableDirectory());
    // User's stuff comes next, if any directories were provided. We're going
    // to turn these into absolute pathnames, interpreting them as defined
    // by Pathname, which includes executable-relative names. The "bin"
    // subdirectory if any must already be present in the directory names.
    for (unsigned i=0; i < userSearchPath.size(); ++i)
        actualSearchPath.push_back
           (Pathname::getAbsoluteDirectoryPathname(userSearchPath[i]));

    if (Pathname::environmentVariableExists("SIMBODY_HOME")) {
        const std::string e = Pathname::getAbsoluteDirectoryPathname(
                Pathname::getEnvironmentVariable("SIMBODY_HOME"));
        actualSearchPath.push_back(Pathname::addDirectoryOffset(e,"bin"));
    } else if (Pathname::environmentVariableExists("SimTK_INSTALL_DIR")) {
        const std::string e = Pathname::getAbsoluteDirectoryPathname(
            Pathname::getEnvironmentVariable("SimTK_INSTALL_DIR"));
        actualSearchPath.push_back(Pathname::addDirectoryOffset(e,"bin"));
    } else {
        // No environment variables set. Our last desperate attempts will
        // be  <platformDefaultInstallDir>/Simbody/bin
        // and <platformDefaultInstallDir>/SimTK/bin
        const std::string def = Pathname::getDefaultInstallDir();

        actualSearchPath.push_back(
            Pathname::addDirectoryOffset(def,
                Pathname::addDirectoryOffset("Simbody", "bin")));
        actualSearchPath.push_back(
            Pathname::addDirectoryOffset(def,
                Pathname::addDirectoryOffset("SimTK", "bin")));
    }

    int sim2vizPipe[2], viz2simPipe[2], status;

    // Create pipe pair for communication from simulator to visualizer.
    status = createPipe(sim2vizPipe);
    SimTK_ASSERT_ALWAYS(status != -1, "VisualizerProtocol: Failed to open pipe");
    outPipe = sim2vizPipe[1];

    // Create pipe pair for communication from visualizer to simulator.
    status = createPipe(viz2simPipe);
    SimTK_ASSERT_ALWAYS(status != -1, "VisualizerProtocol: Failed to open pipe");
    inPipe = viz2simPipe[0];

    // Spawn the visualizer gui, trying local first then installed version.
    spawnViz(actualSearchPath,
             GuiAppName, sim2vizPipe[0], viz2simPipe[1]);

    // Before we do anything else, attempt to exchange handshake messages with
    // the VisualizerGUI. This will throw an exception if anything goes wrong.
    // Note that this is done on the main thread.
    shakeHandsWithGUI(outPipe, inPipe);

    // Spawn the thread to listen for events.

    pthread_mutex_init(&sceneLock, NULL);
    pthread_t thread;
    pthread_create(&thread, NULL, listenForVisualizerEvents, &visualizer);
}

// This is executed on the main thread at GUI startup and thus does not
// require locking.
void VisualizerProtocol::shakeHandsWithGUI(int toGUIPipe, int fromGUIPipe) {

        // First send handshake message to GUI.

    // The first two items must never change: the handshake command value, and
    // an unsigned int containing the version number. Anything else might vary so both
    // sides must stop reading if the version numbers are not compatible.
    WRITE(outPipe, &StartupHandshake, 1);
    WRITE(outPipe, &ProtocolVersion, sizeof(unsigned int));

    // Send the current Simbody version number.
    int major, minor, patch;
    SimTK_version_simbody(&major, &minor, &patch);
    WRITE(outPipe, &major, sizeof(int));
    WRITE(outPipe, &minor, sizeof(int));
    WRITE(outPipe, &patch, sizeof(int));

    // Send the name of the current executable for use as a default window
    // title and in "about" info.
    bool isAbsolutePath;
    std::string directory, fileName, extension;
    Pathname::deconstructPathname(Pathname::getThisExecutablePath(),
        isAbsolutePath, directory, fileName, extension);
    // We're just sending the file name, not a full path. Keep it short.
    unsigned nameLength = std::min((unsigned)fileName.size(), (unsigned)255);
    WRITE(outPipe, &nameLength, sizeof(unsigned));
    WRITE(outPipe, fileName.c_str(), nameLength);

        // Now wait for handshake response from GUI.

    unsigned char handshakeCommand;
    readDataFromPipe(fromGUIPipe, &handshakeCommand, 1);
    SimTK_ERRCHK2_ALWAYS(handshakeCommand == ReturnHandshake,
        "VisualizerProtocol::shakeHandsWithGUI()",
        "Expected initial handshake command %u but received %u. Can't continue.",
        (unsigned)ReturnHandshake, (unsigned)handshakeCommand);

    unsigned int GUIversion;
    readDataFromPipe(fromGUIPipe, (unsigned char*)&GUIversion, sizeof(unsigned int));
    SimTK_ERRCHK2_ALWAYS(GUIversion == ProtocolVersion,
        "VisualizerProtocol::shakeHandsWithGUI()",
        "VisualizerGUI protocol version %u is not compatible with the Simbody"
        " Visualizer class protocol %u; this may be an installation problem."
        " Can't continue.",
        GUIversion, ProtocolVersion);

    // Handshake was successful.
}


void VisualizerProtocol::beginScene(Real time) {
    pthread_mutex_lock(&sceneLock);
    char command = StartOfScene;
    WRITE(outPipe, &command, 1);
    float fTime = (float)time;
    WRITE(outPipe, &fTime, sizeof(float));
}

void VisualizerProtocol::finishScene() {
    char command = EndOfScene;
    WRITE(outPipe, &command, 1);
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::drawBox(const Transform& X_GB, const Vec3& scale, const Vec4& color, int representation) {
    drawMesh(X_GB, scale, color, (short) representation, MeshBox, 0);
}

void VisualizerProtocol::drawEllipsoid(const Transform& X_GB, const Vec3& scale, const Vec4& color, int representation, unsigned short resolution) {
    drawMesh(X_GB, scale, color, (short) representation, MeshEllipsoid, resolution);
}

void VisualizerProtocol::drawCylinder(const Transform& X_GB, const Vec3& scale, const Vec4& color, int representation, unsigned short resolution) {
    drawMesh(X_GB, scale, color, (short) representation, MeshCylinder, resolution);
}

void VisualizerProtocol::drawCircle(const Transform& X_GB, const Vec3& scale, const Vec4& color, int representation, unsigned short resolution) {
    drawMesh(X_GB, scale, color, (short) representation, MeshCircle, resolution);
}

void VisualizerProtocol::drawPolygonalMesh(const PolygonalMesh& mesh, const Transform& X_GM, const Vec3& scale, const Vec4& color, int representation) {
    const void* impl = &mesh.getImpl();
    map<const void*, unsigned short>::const_iterator iter = meshes.find(impl);

    if (iter != meshes.end()) {
        // This mesh was already cached; just reference it by index number.
        drawMesh(X_GM, scale, color, (short)representation, iter->second, 0);
        return;
    }

    // This is a new mesh, so we need to send it to the visualizer. Build lists
    // of vertices and faces, triangulating as necessary.

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
            for (int j = 0; j < numVert; j++) {
                Vec3 pos = mesh.getVertexPosition(mesh.getFaceVertex(i,j));
                center += pos;
            }
            center /= numVert;
            vertices.push_back((float) center[0]);
            vertices.push_back((float) center[1]);
            vertices.push_back((float) center[2]);
            const unsigned newIndex = vertices.size()/3-1;
            for (int j = 0; j < numVert-1; j++) {
                faces.push_back((unsigned short) mesh.getFaceVertex(i, j));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, j+1));
                faces.push_back((unsigned short) newIndex);
            }
        }
    }
    SimTK_ERRCHK1_ALWAYS(vertices.size() <= 65535*3, 
        "VisualizerProtocol::drawPolygonalMesh()",
        "Can't display a DecorativeMesh with more than 65535 vertices;"
        " received one with %llu.", (unsigned long long)vertices.size());
    SimTK_ERRCHK1_ALWAYS(faces.size() <= 65535*3, 
        "VisualizerProtocol::drawPolygonalMesh()",
        "Can't display a DecorativeMesh with more than 65535 vertices;"
        " received one with %llu.", (unsigned long long)faces.size());

    const int index = NumPredefinedMeshes + (int)meshes.size();
    SimTK_ERRCHK_ALWAYS(index <= 65535,
        "VisualizerProtocol::drawPolygonalMesh()",
        "Too many unique DecorativeMesh objects; max is 65535.");
    
    meshes[impl] = (unsigned short)index;    // insert new mesh
    WRITE(outPipe, &DefineMesh, 1);
    unsigned short numVertices = (unsigned)vertices.size()/3;
    unsigned short numFaces = (unsigned)faces.size()/3;
    WRITE(outPipe, &numVertices, sizeof(short));
    WRITE(outPipe, &numFaces, sizeof(short));
    WRITE(outPipe, &vertices[0], (unsigned)(vertices.size()*sizeof(float)));
    WRITE(outPipe, &faces[0], (unsigned)(faces.size()*sizeof(short)));

    drawMesh(X_GM, scale, color, (short) representation, index, 0);
}

void VisualizerProtocol::
drawMesh(const Transform& X_GM, const Vec3& scale, const Vec4& color, 
         short representation, unsigned short meshIndex, unsigned short resolution)
{
    char command = (representation == DecorativeGeometry::DrawPoints 
                    ? AddPointMesh 
                    : (representation == DecorativeGeometry::DrawWireframe 
                        ? AddWireframeMesh : AddSolidMesh));
    WRITE(outPipe, &command, 1);
    float buffer[13];
    Vec3 rot = X_GM.R().convertRotationToBodyFixedXYZ();
    buffer[0] = (float) rot[0];
    buffer[1] = (float) rot[1];
    buffer[2] = (float) rot[2];
    buffer[3] = (float) X_GM.p()[0];
    buffer[4] = (float) X_GM.p()[1];
    buffer[5] = (float) X_GM.p()[2];
    buffer[6] = (float) scale[0];
    buffer[7] = (float) scale[1];
    buffer[8] = (float) scale[2];
    buffer[9] = (float) color[0];
    buffer[10] = (float) color[1];
    buffer[11] = (float) color[2];
    buffer[12] = (float) color[3];
    WRITE(outPipe, buffer, 13*sizeof(float));
    unsigned short buffer2[2];
    buffer2[0] = meshIndex;
    buffer2[1] = resolution;
    WRITE(outPipe, buffer2, 2*sizeof(unsigned short));
}

void VisualizerProtocol::
drawLine(const Vec3& end1, const Vec3& end2, const Vec4& color, Real thickness)
{
    WRITE(outPipe, &AddLine, 1);
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
    WRITE(outPipe, buffer, 10*sizeof(float));
}

void VisualizerProtocol::
drawText(const Vec3& position, const Vec3& scale, const Vec4& color, 
         const string& string, bool faceCamera) {
    SimTK_ERRCHK1_ALWAYS(string.size() <= 256,
        "VisualizerProtocol::drawText()",
        "Can't display DecorativeText longer than 256 characters;"
        " received text of length %u.", (unsigned)string.size());
    WRITE(outPipe, &AddText, 1);
    float buffer[9];
    buffer[0] = (float) position[0];
    buffer[1] = (float) position[1];
    buffer[2] = (float) position[2];
    buffer[3] = (float) scale[0];
    buffer[4] = (float) scale[1];
    buffer[5] = (float) scale[2];
    buffer[6] = (float) color[0];
    buffer[7] = (float) color[1];
    buffer[8] = (float) color[2];
    WRITE(outPipe, buffer, 9*sizeof(float));
    short face = (short)faceCamera;
    WRITE(outPipe, &face, sizeof(short));
    short length = (short)string.size();
    WRITE(outPipe, &length, sizeof(short));
    WRITE(outPipe, &string[0], length);
}

void VisualizerProtocol::
drawCoords(const Transform& X_GF, const Vec3& axisLengths, const Vec4& color) {
    WRITE(outPipe, &AddCoords, 1);
    float buffer[12];
    Vec3 rot = X_GF.R().convertRotationToBodyFixedXYZ();
    buffer[0] = (float) rot[0];
    buffer[1] = (float) rot[1];
    buffer[2] = (float) rot[2];
    buffer[3] = (float) X_GF.p()[0];
    buffer[4] = (float) X_GF.p()[1];
    buffer[5] = (float) X_GF.p()[2];
    buffer[6] = (float) axisLengths[0];
    buffer[7] = (float) axisLengths[1];
    buffer[8] = (float) axisLengths[2];
    buffer[9] = (float) color[0];
    buffer[10]= (float) color[1];
    buffer[11]= (float) color[2];
    WRITE(outPipe, buffer, 12*sizeof(float));
}

void VisualizerProtocol::
addMenu(const String& title, int id, const Array_<pair<String, int> >& items) {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &DefineMenu, 1);
    short titleLength = title.size();
    WRITE(outPipe, &titleLength, sizeof(short));
    WRITE(outPipe, title.c_str(), titleLength);
    WRITE(outPipe, &id, sizeof(int));
    short numItems = items.size();
    WRITE(outPipe, &numItems, sizeof(short));
    for (int i = 0; i < numItems; i++) {
        int buffer[] = {items[i].second, items[i].first.size()};
        WRITE(outPipe, buffer, 2*sizeof(int));
        WRITE(outPipe, items[i].first.c_str(), items[i].first.size());
    }
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::
addSlider(const String& title, int id, Real minVal, Real maxVal, Real value) {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &DefineSlider, 1);
    short titleLength = title.size();
    WRITE(outPipe, &titleLength, sizeof(short));
    WRITE(outPipe, title.c_str(), titleLength);
    WRITE(outPipe, &id, sizeof(int));
    float buffer[3];
    buffer[0] = (float) minVal;
    buffer[1] = (float) maxVal;
    buffer[2] = (float) value;
    WRITE(outPipe, buffer, 3*sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}


void VisualizerProtocol::setSliderValue(int id, Real newValue) const {
    const float value = (float)newValue;
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetSliderValue, 1);
    WRITE(outPipe, &id, sizeof(int));
    WRITE(outPipe, &value, sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setSliderRange(int id, Real newMin, Real newMax) const {
    float buffer[2];
    buffer[0] = (float)newMin; buffer[1] = (float)newMax;
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetSliderRange, 1);
    WRITE(outPipe, &id, sizeof(int));
    WRITE(outPipe, buffer, 2*sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setWindowTitle(const String& title) const {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetWindowTitle, 1);
    short titleLength = title.size();
    WRITE(outPipe, &titleLength, sizeof(short));
    WRITE(outPipe, title.c_str(), titleLength);
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setMaxFrameRate(Real rate) const {
    const float frameRate = (float)rate;
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetMaxFrameRate, 1);
    WRITE(outPipe, &frameRate, sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}


void VisualizerProtocol::setBackgroundColor(const Vec3& color) const {
    float buffer[3];
    buffer[0] = (float)color[0]; 
    buffer[1] = (float)color[1]; 
    buffer[2] = (float)color[2];
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetBackgroundColor, 1);
    WRITE(outPipe, buffer, 3*sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setShowShadows(bool shouldShow) const {
    const short show = (short)shouldShow; // 0 or 1
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetShowShadows, 1);
    WRITE(outPipe, &show, sizeof(short));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setBackgroundType(Visualizer::BackgroundType type) const {
    const short backgroundType = (short)type;
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetBackgroundType, 1);
    WRITE(outPipe, &backgroundType, sizeof(short));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setCameraTransform(const Transform& X_GC) const {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetCamera, 1);
    float buffer[6];
    Vec3 rot = X_GC.R().convertRotationToBodyFixedXYZ();
    buffer[0] = (float) rot[0];
    buffer[1] = (float) rot[1];
    buffer[2] = (float) rot[2];
    buffer[3] = (float) X_GC.p()[0];
    buffer[4] = (float) X_GC.p()[1];
    buffer[5] = (float) X_GC.p()[2];
    WRITE(outPipe, buffer, 6*sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::zoomCamera() const {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &ZoomCamera, 1);
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::lookAt(const Vec3& point, const Vec3& upDirection) const {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &LookAt, 1);
    float buffer[6];
    buffer[0] = (float) point[0];
    buffer[1] = (float) point[1];
    buffer[2] = (float) point[2];
    buffer[3] = (float) upDirection[0];
    buffer[4] = (float) upDirection[1];
    buffer[5] = (float) upDirection[2];
    WRITE(outPipe, buffer, 6*sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setFieldOfView(Real fov) const {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetFieldOfView, 1);
    float buffer[1];
    buffer[0] = (float)fov;
    WRITE(outPipe, buffer, sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setClippingPlanes(Real near, Real far) const {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetClipPlanes, 1);
    float buffer[2];
    buffer[0] = (float)near;
    buffer[1] = (float)far;
    WRITE(outPipe, buffer, 2*sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::
setSystemUpDirection(const CoordinateDirection& upDir) {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetSystemUpDirection, 1);
    const unsigned char axis = (unsigned char)upDir.getAxis();
    const signed char   sign = (signed char)upDir.getDirection();
    WRITE(outPipe, &axis, 1);
    WRITE(outPipe, &sign, 1);
    pthread_mutex_unlock(&sceneLock);
}

void VisualizerProtocol::setGroundHeight(Real height) {
    pthread_mutex_lock(&sceneLock);
    WRITE(outPipe, &SetGroundHeight, 1);
    float heightBuffer = (float) height;
    WRITE(outPipe, &heightBuffer, sizeof(float));
    pthread_mutex_unlock(&sceneLock);
}


