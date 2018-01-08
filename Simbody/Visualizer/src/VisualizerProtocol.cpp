/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-15 Stanford University and the Authors.        *
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
    #define WRITEFUNC _write
    #define CLOSE _close
#else
    #include <unistd.h>
    #define READ read
    #define WRITEFUNC write
    #define CLOSE close
#endif

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strcat, sprintf, etc.
#endif

// gcc 4.4.3 complains bitterly if you don't check the return
// status from the write() system call. This avoids those 
// warnings and maybe, someday, will catch an error.
#define WRITE(pipeno, buf, len) \
   {int status=WRITEFUNC((pipeno), (buf), (len)); \
    SimTK_ERRCHK4_ALWAYS(status!=-1, "VisualizerProtocol",  \
    "An attempt to write() %d bytes to pipe %d failed with errno=%d (%s).", \
    (len),(pipeno),errno,strerror(errno));}

static int inPipe;

// Create the pipe going *from* simulator *to* visualizer, so only the
// read end should be inherited by the visualizer.
static int createPipeSim2Viz(int sim2viz[2]) {
    const int status =
#ifdef _WIN32
        _pipe(sim2viz, 16384, _O_BINARY|_O_NOINHERIT);
        if (status==0) {
            const int readfd = _dup(sim2viz[0]); // copy is inheritable
            _close(sim2viz[0]); // close the non-inheritable version
            sim2viz[0] = readfd; // replace with inheritable copy
        }
#else
        pipe(sim2viz);
#endif
    return status;
}


// Create the pipe going *to* simulator *from* visualizer, so only the
// write end should be inherited by the visualizer.
static int createPipeViz2Sim(int viz2sim[2]) {
    const int status =
#ifdef _WIN32
        _pipe(viz2sim, 16384, _O_BINARY|_O_NOINHERIT);
        if (status==0) {
            const int writefd = _dup(viz2sim[1]); // copy is inheritable
            _close(viz2sim[1]); // close the non-inheritable version
            viz2sim[1] = writefd; // replace with inheritable copy
        }
#else
        pipe(viz2sim);
#endif
    return status;
}

// Spawn the visualizer GUI executable, using the right method for
// this platform. We take two executables to try in order,
// and return after the first one succeeds. If neither works, we throw
// an error that is hopefully helful.
static void spawnViz(const Array_<String>& searchPath, const String& appName, 
                     int sim2vizPipe[2], int viz2simPipe[2])
{
    int status;

    // Pass pipe numbers as command line arguments to the visualizer.
    char vizReadFromSim[32], vizWriteToSim[32];
    sprintf(vizReadFromSim, "%d", sim2vizPipe[0]);
    sprintf(vizWriteToSim, "%d", viz2simPipe[1]);

    String exePath; // search path + appName

#ifdef _WIN32
    intptr_t handle;
    for (unsigned i=0; i < searchPath.size(); ++i) {
        exePath = searchPath[i] + appName;
        handle = _spawnl(P_NOWAIT, exePath.c_str(), appName.c_str(), 
                         vizReadFromSim, vizWriteToSim, (const char*)0);
        if (handle != -1) {
            // success: visualizer is running
            _close(sim2vizPipe[0]); // read end (belongs to visualizer)
            _close(viz2simPipe[1]); // write end (belongs to visualizer)
            break; // success!
        }
    }
    status = (handle==-1) ? -1 : 0;
#else
    const pid_t pid = fork();
    if (pid == 0) {
        // child process (the visualizer)
        close(sim2vizPipe[1]); // write end (belongs to simulator)
        close(viz2simPipe[0]); // read end(belongs to simulator)
        for (unsigned i=0; i < searchPath.size(); ++i) {
            exePath = searchPath[i] + appName;
            status = execl(exePath.c_str(), appName.c_str(), 
                           vizReadFromSim, vizWriteToSim, (const char*)0); 
            // if we get here the execl() failed
        }
        // fall through -- we failed on every try
    } else {
        // parent process (the simulator)
        status = (pid==-1) ? -1 : 0;
        if (status == 0) {
            close(sim2vizPipe[0]); // read end (belongs to visualizer)
            close(viz2simPipe[1]); // write end (belongs to visualizer)
        }
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

class ReadingInterrupted : public std::exception {};

// This will hang until the expected number of bytes has been received.
// Throws ReadingInterrupted if the srcPipe is closed.
static void readDataFromPipe(int srcPipe, unsigned char* buffer, int bytes) {
    int totalRead = 0;
    while (totalRead < bytes) {
        auto retval = READ(srcPipe, buffer + totalRead, bytes - totalRead);
        SimTK_ERRCHK4_ALWAYS(retval!=-1, "VisualizerProtocol",
            "An attempt to read() %d bytes from pipe %d failed with errno=%d (%s).", 
            bytes - totalRead, srcPipe, errno, strerror(errno));
        // The pipe was closed, perhaps because simbody-visualizer was closed,
        // or shutdownGUI() or ~VisualizerProtocol() was called.
        // Without this check, we end up in an infinite loop if the user closes
        // simbody-visualizer.
        if (retval == 0) throw ReadingInterrupted();

        totalRead += retval;
    }
}

// Throws ReadingInterrupted if inPipe is closed.
static void readData(unsigned char* buffer, int bytes)
{
    readDataFromPipe(inPipe, buffer, bytes);
}

static void listenForVisualizerEvents(Visualizer& visualizer) {
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
                "Unexpected command %u received from simbody-visualizer. Can't continue.",
                (unsigned)buffer[0]);
        }
    }
  } catch (const ReadingInterrupted&) {
        // We were told (by the main thread) to stop listening, or
        // simbody-visualizer closed.
  } catch (const std::exception& e) {
        std::cout << "Visualizer listenerThread: unrecoverable error:\n";
        std::cout << e.what() << std::endl;
    }
}

VisualizerProtocol::VisualizerProtocol
   (Visualizer& visualizer, const Array_<String>& userSearchPath) 
{
    // Launch the GUI application. We'll first look for one in the same
    // directory as the running executable; then if that doesn't work we'll
    // look in the bin subdirectory of the SimTK installation.

    String vizExecutableName;
    if (Pathname::environmentVariableExists("SIMBODY_VISUALIZER_NAME")) {
        vizExecutableName =
            Pathname::getEnvironmentVariable("SIMBODY_VISUALIZER_NAME");
    } else {
        vizExecutableName = "simbody-visualizer";
        #ifndef NDEBUG
            vizExecutableName += "_d";
        #endif
    }

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
        actualSearchPath.push_back(Pathname::addDirectoryOffset(e,
                    SIMBODY_VISUALIZER_REL_INSTALL_DIR));
    } else if (Pathname::environmentVariableExists("SimTK_INSTALL_DIR")) {
        const std::string e = Pathname::getAbsoluteDirectoryPathname(
            Pathname::getEnvironmentVariable("SimTK_INSTALL_DIR"));
        actualSearchPath.push_back(Pathname::addDirectoryOffset(e,
                    SIMBODY_VISUALIZER_REL_INSTALL_DIR));
    }

    // Try using the location of SimTKsimbody combined with the path from the
    // SimTKsimbody library to the simbody-visualizer install location (only
    // available on non-Windows platforms). We must provide a function that
    // resides in SimTKsimbody.
    std::string SimTKsimbodyDir;
    if (Pathname::getFunctionLibraryDirectory((void*)createPipeSim2Viz,
                                              SimTKsimbodyDir)) {
        // We have the path to SimTKsimbody; now we combine this with the path
        // from SimTKsimbody to simbody-visualizer (this assumes the
        // installation has not been reorganized from what CMake originally
        // specified).
        std::string absPathToVizDir =
          Pathname::getAbsoluteDirectoryPathnameUsingSpecifiedWorkingDirectory(
                  SimTKsimbodyDir, SIMBODY_PATH_FROM_LIBDIR_TO_VIZ_DIR);
        actualSearchPath.push_back(absPathToVizDir);
    }

    // Try the build-time install location:
    actualSearchPath.push_back(SIMBODY_VISUALIZER_INSTALL_DIR);

    // Our last desperate attempts will
    // be  <platformDefaultInstallDir>/Simbody/bin
    // and <platformDefaultInstallDir>/SimTK/bin
    const std::string def = Pathname::getDefaultInstallDir();

    actualSearchPath.push_back(
        Pathname::addDirectoryOffset(def,
            Pathname::addDirectoryOffset("Simbody", SIMBODY_VISUALIZER_REL_INSTALL_DIR)));
    actualSearchPath.push_back(
        Pathname::addDirectoryOffset(def,
            Pathname::addDirectoryOffset("SimTK", SIMBODY_VISUALIZER_REL_INSTALL_DIR)));

    // Pipe[0] is the read end, Pipe[1] is the write end.
    int sim2vizPipe[2], viz2simPipe[2], status;

    // Create pipe pair for communication from simulator to visualizer.
    status = createPipeSim2Viz(sim2vizPipe);
    SimTK_ASSERT_ALWAYS(status != -1, "VisualizerProtocol: Failed to open pipe");
    outPipe = sim2vizPipe[1]; // write here to talk to visualizer

    // Create pipe pair for communication from visualizer to simulator.
    status = createPipeViz2Sim(viz2simPipe);
    SimTK_ASSERT_ALWAYS(status != -1, "VisualizerProtocol: Failed to open pipe");
    inPipe = viz2simPipe[0]; // read from here to receive from visualizer

    // Spawn the visualizer gui, trying local first then installed version.
    spawnViz(actualSearchPath, vizExecutableName, sim2vizPipe, viz2simPipe);

    // Before we do anything else, attempt to exchange handshake messages with
    // the visualizer. This will throw an exception if anything goes wrong.
    // Note that this is done on the main thread.
    shakeHandsWithGUI(outPipe, inPipe);

    // Spawn the thread to listen for events.
    eventListenerThread = std::thread(listenForVisualizerEvents,
            std::ref(visualizer));
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
        "simbody-visualizer protocol version %u is not compatible with the Simbody"
        " Visualizer class protocol %u; this may be an installation problem."
        " Can't continue.",
        GUIversion, ProtocolVersion);

    // Handshake was successful.
}

void VisualizerProtocol::shutdownGUI() {
    // Don't wait for scene completion; kill GUI now.
    
    // We no longer need to listen for events from the GUI. Stop the listener
    // thread before writing to the out pipe, because attempting to write to
    // the pipe may throw an exception. Shutting down the listener thread was
    // added to solve an issue with OpenSim MATLAB bindings, wherein MATLAB
    // would use more and more CPU each time a Visualizer was created.
    stopListeningIfNecessary();
    
    char command = Shutdown;
    WRITE(outPipe, &command, 1);
}

VisualizerProtocol::~VisualizerProtocol() {
    // If shutdownGUI() was not called, then the listener thread is still
    // running and we should kill it.
    stopListeningIfNecessary();
    int retval = CLOSE(outPipe); // TODO(chrisdembia) is this necessary?
    if (retval == -1) {
        std::cout << "Warning in Simbody VisualizerProtocol: "
            << "An attempt to close() pipe " << outPipe
            << " failed with errno=" << errno << " (" << strerror(errno) << ")."
            << std::endl;
    }
}

void VisualizerProtocol::stopListeningIfNecessary() {
    if (eventListenerThread.joinable()) {
        // Shut down the listener thread cleanly. Tell the GUI to tell the
        // simulator's listener thread to stop listening, which will allow the
        // the (simulator's) listener thread to die.
        WRITE(outPipe, &StopCommunication, 1);
        eventListenerThread.join();
    }
}

void VisualizerProtocol::beginScene(Real time) {
    sceneLockBeginFinishScene.lock();
    char command = StartOfScene;
    WRITE(outPipe, &command, 1);
    float fTime = (float)time;
    WRITE(outPipe, &fTime, sizeof(float));
    // The sceneMutex is NOT unlocked at the end of this scope
    // (sceneLockBeginFinishScene is a member variable); see finishScene().
}

void VisualizerProtocol::finishScene() {
    char command = EndOfScene;
    WRITE(outPipe, &command, 1);
    sceneLockBeginFinishScene.unlock();
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
            const unsigned newIndex = (unsigned)(vertices.size()/3-1);
            for (int j = 0; j < numVert-1; j++) {
                faces.push_back((unsigned short) mesh.getFaceVertex(i, j));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, j+1));
                faces.push_back((unsigned short) newIndex);
            }
            // Close the face (thanks, Alexandra Zobova).
            faces.push_back((unsigned short) mesh.getFaceVertex(i, numVert-1));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
            faces.push_back((unsigned short) newIndex);
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
    unsigned short numVertices = (unsigned short)(vertices.size()/3);
    unsigned short numFaces = (unsigned short)(faces.size()/3);
    WRITE(outPipe, &numVertices, sizeof(short));
    WRITE(outPipe, &numFaces, sizeof(short));
    WRITE(outPipe, &vertices[0], (unsigned)(vertices.size()*sizeof(float)));
    WRITE(outPipe, &faces[0], (unsigned)(faces.size()*sizeof(short)));

    drawMesh(X_GM, scale, color, (short) representation, (unsigned short)index, 0);
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
drawText(const Transform& X_GT, const Vec3& scale, const Vec4& color, 
         const string& string, bool faceCamera, bool isScreenText) {
    SimTK_ERRCHK1_ALWAYS(string.size() <= 256,
        "VisualizerProtocol::drawText()",
        "Can't display DecorativeText longer than 256 characters;"
        " received text of length %u.", (unsigned)string.size());
    WRITE(outPipe, &AddText, 1);
    float buffer[12];
    const Vec3 rot = X_GT.R().convertRotationToBodyFixedXYZ();
    buffer[0] = (float) rot[0];
    buffer[1] = (float) rot[1];
    buffer[2] = (float) rot[2];
    buffer[3] = (float) X_GT.p()[0];
    buffer[4] = (float) X_GT.p()[1];
    buffer[5] = (float) X_GT.p()[2];
    buffer[6] = (float) scale[0];
    buffer[7] = (float) scale[1];
    buffer[8] = (float) scale[2];
    buffer[9] = (float) color[0];
    buffer[10]= (float) color[1];
    buffer[11]= (float) color[2];
    WRITE(outPipe, buffer, 12*sizeof(float));
    short face = (short)faceCamera;
    WRITE(outPipe, &face, sizeof(short));
    short screen = (short)isScreenText;
    WRITE(outPipe, &screen, sizeof(short));
    short length = (short)string.size();
    WRITE(outPipe, &length, sizeof(short));
    WRITE(outPipe, &string[0], length);
}

void VisualizerProtocol::
drawCoords(const Transform& X_GF, const Vec3& axisLengths, const Vec4& color) {
    WRITE(outPipe, &AddCoords, 1);
    float buffer[12];
    const Vec3 rot = X_GF.R().convertRotationToBodyFixedXYZ();
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
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &DefineMenu, 1);
    short titleLength = (short)title.size();
    WRITE(outPipe, &titleLength, sizeof(short));
    WRITE(outPipe, title.c_str(), titleLength);
    WRITE(outPipe, &id, sizeof(int));
    short numItems = (short)items.size();
    WRITE(outPipe, &numItems, sizeof(short));
    for (int i = 0; i < numItems; i++) {
        int buffer[] = {items[i].second, items[i].first.size()};
        WRITE(outPipe, buffer, 2*sizeof(int));
        WRITE(outPipe, items[i].first.c_str(), items[i].first.size());
    }
}

void VisualizerProtocol::
addSlider(const String& title, int id, Real minVal, Real maxVal, Real value) {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &DefineSlider, 1);
    short titleLength = (short)title.size();
    WRITE(outPipe, &titleLength, sizeof(short));
    WRITE(outPipe, title.c_str(), titleLength);
    WRITE(outPipe, &id, sizeof(int));
    float buffer[3];
    buffer[0] = (float) minVal;
    buffer[1] = (float) maxVal;
    buffer[2] = (float) value;
    WRITE(outPipe, buffer, 3*sizeof(float));
}


void VisualizerProtocol::setSliderValue(int id, Real newValue) const {
    const float value = (float)newValue;
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetSliderValue, 1);
    WRITE(outPipe, &id, sizeof(int));
    WRITE(outPipe, &value, sizeof(float));
}

void VisualizerProtocol::setSliderRange(int id, Real newMin, Real newMax) const {
    float buffer[2];
    buffer[0] = (float)newMin; buffer[1] = (float)newMax;
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetSliderRange, 1);
    WRITE(outPipe, &id, sizeof(int));
    WRITE(outPipe, buffer, 2*sizeof(float));
}

void VisualizerProtocol::setWindowTitle(const String& title) const {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetWindowTitle, 1);
    short titleLength = (short)title.size();
    WRITE(outPipe, &titleLength, sizeof(short));
    WRITE(outPipe, title.c_str(), titleLength);
}

void VisualizerProtocol::setMaxFrameRate(Real rate) const {
    const float frameRate = (float)rate;
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetMaxFrameRate, 1);
    WRITE(outPipe, &frameRate, sizeof(float));
}


void VisualizerProtocol::setBackgroundColor(const Vec3& color) const {
    float buffer[3];
    buffer[0] = (float)color[0]; 
    buffer[1] = (float)color[1]; 
    buffer[2] = (float)color[2];
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetBackgroundColor, 1);
    WRITE(outPipe, buffer, 3*sizeof(float));
}

void VisualizerProtocol::setShowShadows(bool shouldShow) const {
    const short show = (short)shouldShow; // 0 or 1
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetShowShadows, 1);
    WRITE(outPipe, &show, sizeof(short));
}

void VisualizerProtocol::setShowFrameRate(bool shouldShow) const {
    const short show = (short)shouldShow; // 0 or 1
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetShowFrameRate, 1);
    WRITE(outPipe, &show, sizeof(short));
}

void VisualizerProtocol::setShowSimTime(bool shouldShow) const {
    const short show = (short)shouldShow; // 0 or 1
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetShowSimTime, 1);
    WRITE(outPipe, &show, sizeof(short));
}

void VisualizerProtocol::setShowFrameNumber(bool shouldShow) const {
    const short show = (short)shouldShow; // 0 or 1
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetShowFrameNumber, 1);
    WRITE(outPipe, &show, sizeof(short));
}

void VisualizerProtocol::setBackgroundType(Visualizer::BackgroundType type) const {
    const short backgroundType = (short)type;
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetBackgroundType, 1);
    WRITE(outPipe, &backgroundType, sizeof(short));
}

void VisualizerProtocol::setCameraTransform(const Transform& X_GC) const {
    std::lock_guard<std::mutex> lock(sceneMutex);
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
}

void VisualizerProtocol::zoomCamera() const {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &ZoomCamera, 1);
}

void VisualizerProtocol::lookAt(const Vec3& point, const Vec3& upDirection) const {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &LookAt, 1);
    float buffer[6];
    buffer[0] = (float) point[0];
    buffer[1] = (float) point[1];
    buffer[2] = (float) point[2];
    buffer[3] = (float) upDirection[0];
    buffer[4] = (float) upDirection[1];
    buffer[5] = (float) upDirection[2];
    WRITE(outPipe, buffer, 6*sizeof(float));
}

void VisualizerProtocol::setFieldOfView(Real fov) const {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetFieldOfView, 1);
    float buffer[1];
    buffer[0] = (float)fov;
    WRITE(outPipe, buffer, sizeof(float));
}

void VisualizerProtocol::setClippingPlanes(Real near, Real far) const {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetClipPlanes, 1);
    float buffer[2];
    buffer[0] = (float)near;
    buffer[1] = (float)far;
    WRITE(outPipe, buffer, 2*sizeof(float));
}

void VisualizerProtocol::
setSystemUpDirection(const CoordinateDirection& upDir) {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetSystemUpDirection, 1);
    const unsigned char axis = (unsigned char)upDir.getAxis();
    const signed char   sign = (signed char)upDir.getDirection();
    WRITE(outPipe, &axis, 1);
    WRITE(outPipe, &sign, 1);
}

void VisualizerProtocol::setGroundHeight(Real height) {
    std::lock_guard<std::mutex> lock(sceneMutex);
    WRITE(outPipe, &SetGroundHeight, 1);
    float heightBuffer = (float) height;
    WRITE(outPipe, &heightBuffer, sizeof(float));
}


