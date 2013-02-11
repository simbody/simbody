#include "VisualizerBase.h"
#include <vector>
#include <utility>

using namespace std;
using namespace SimTK;

VisualizerBase::VisualizerBase()
{
	dummyDisplay = NULL;
	scene = NULL;
}

VisualizerBase::~VisualizerBase()
{
}

void VisualizerBase::shakeHandsWithSimulator(int inPipe, int outPipe)
{
	unsigned char handshakeCommand;
	int simbodyVersion[3];

    readDataFromPipe(inPipe, &handshakeCommand, 1);
    SimTK_ERRCHK2_ALWAYS(handshakeCommand == StartupHandshake,
        "VisualizerGUI::shakeHandsWithSimulator()",
        "Expected initial handshake command %u but received %u. Can't continue.",
        (unsigned)StartupHandshake, (unsigned)handshakeCommand);

    unsigned SimVersion;
    readDataFromPipe(inPipe, (unsigned char*)&SimVersion, sizeof(unsigned int));
    SimTK_ERRCHK2_ALWAYS(SimVersion == ProtocolVersion,
        "VisualizerGUI::shakeHandsWithSimulator()",
        "The Simbody Visualizer class protocol version %u is not compatible with "
        " VisualizerGUI protocol %u; this may be an installation problem."
        " Can't continue.",
        SimVersion, ProtocolVersion);

    // Get Simbody version number as major,minor,patch
    readDataFromPipe(inPipe, (unsigned char*)simbodyVersion, 3*sizeof(int));
    simbodyVersionStr = String(simbodyVersion[0]) + "." + String(simbodyVersion[1]);
    if (simbodyVersion[2]) simbodyVersionStr += "." + String(simbodyVersion[2]);

    unsigned exeNameLength;
    char exeNameBuf[256]; // just a file name, not a path name
    readDataFromPipe(inPipe, (unsigned char*)&exeNameLength, sizeof(unsigned));
    SimTK_ASSERT_ALWAYS(exeNameLength <= 255,
        "VisualizerGUI: executable name length violates protocol.");
    readDataFromPipe(inPipe, (unsigned char*)exeNameBuf, exeNameLength);
    exeNameBuf[exeNameLength] = (char)0;

    simulatorExecutableName = std::string(exeNameBuf, exeNameLength);

    WRITE(outPipe, &ReturnHandshake, 1);
    WRITE(outPipe, &ProtocolVersion, sizeof(unsigned));
}

void VisualizerBase::readDataFromPipe(int srcPipe, unsigned char* buffer, int bytes)
{
    int totalRead = 0;
    while (totalRead < bytes)
        totalRead += READ(srcPipe, buffer+totalRead, bytes-totalRead);

}

void VisualizerBase::setName(std::string appName)
{
	this->name = appName;
}

std::string VisualizerBase::getName()
{
	return this->name;
}

std::string VisualizerBase::getVersion()
{
	return simbodyVersionStr;
}

std::string VisualizerBase::getExecutableName()
{
	return simulatorExecutableName;
}

// We have just processed a StartOfScene command. Read in all the scene
// elements until we see an EndOfScene command. We allocate a new Scene
// object to hold the scene and return a pointer to it. Don't forget to
// delete that object when you are done with it.
Scene* VisualizerBase::readNewScene(int inPipe) {
    unsigned char buffer[256];
    float*          floatBuffer = (float*)          buffer;
    int*            intBuffer   = (int*)            buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;

    Scene* newScene = new Scene;

    // Simulated time for this frame comes first.
    readDataFromPipe(inPipe, buffer, sizeof(float));
    newScene->simTime = floatBuffer[0];

    bool finished = false;
    while (!finished) {
        readDataFromPipe(inPipe, buffer, 1);
        char command = buffer[0];

        switch (command) {

        case EndOfScene:
            finished = true;
            break;

        // Add a scene element that uses an already-cached mesh.
        case AddPointMesh:
        case AddWireframeMesh:
        case AddSolidMesh: {
            readDataFromPipe(inPipe, buffer, 13*sizeof(float)+2*sizeof(short));
            fTransform position;
            position.updR().setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]));
            position.updP() = fVec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 scale = fVec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            fVec4 color = fVec4(floatBuffer[9], floatBuffer[10], floatBuffer[11], floatBuffer[12]);
            short representation = (command == AddPointMesh ? DecorativeGeometry::DrawPoints : (command == AddWireframeMesh ? DecorativeGeometry::DrawWireframe : DecorativeGeometry::DrawSurface));
            unsigned short meshIndex = shortBuffer[13*sizeof(float)/sizeof(short)];
            unsigned short resolution = shortBuffer[13*sizeof(float)/sizeof(short)+1];
            RenderedMesh mesh(position, scale, color, representation, meshIndex, resolution);
            if (command != AddSolidMesh)
                newScene->drawnMeshes.push_back(mesh);
            else if (color[3] == 1)
                newScene->solidMeshes.push_back(mesh);
            else
                newScene->transparentMeshes.push_back(mesh);
//            if (meshIndex < NumPredefinedMeshes && (meshes[meshIndex].size() <= resolution || meshes[meshIndex][resolution] == NULL)) {
                // A real mesh will be generated from this the next
                // time the scene is redrawn.
//                pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE --------
//                pendingCommands.insert(pendingCommands.begin(), new PendingStandardMesh(meshIndex, resolution));
//                pthread_mutex_unlock(&sceneLock);   //------ UNLOCK SCENE --------
//            }
            break;
        }

        case AddLine: {
            readDataFromPipe(inPipe, buffer, 10*sizeof(float));
            fVec3 color = fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            float thickness = floatBuffer[3];
            int index;
            int numLines = (int)newScene->lines.size();
            for (index = 0; index < numLines && (color != newScene->lines[index].getColor() || thickness != newScene->lines[index].getThickness()); index++)
                ;
            if (index == numLines)
                newScene->lines.push_back(RenderedLine(color, thickness));
            vector<float>& line = newScene->lines[index].getLines();
            line.push_back(floatBuffer[4]);
            line.push_back(floatBuffer[5]);
            line.push_back(floatBuffer[6]);
            line.push_back(floatBuffer[7]);
            line.push_back(floatBuffer[8]);
            line.push_back(floatBuffer[9]);
            break;
        }

        case AddText: {
            readDataFromPipe(inPipe, buffer, 9*sizeof(float)+3*sizeof(short));
            fVec3 position = fVec3(floatBuffer[0], floatBuffer[1], floatBuffer[2]);
            fVec3 scale = fVec3(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 color = fVec3(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            unsigned short* shortp = &shortBuffer[9*sizeof(float)/sizeof(short)];
            bool faceCamera = (shortp[0] != 0);
            bool isScreenText = (shortp[1] != 0);
            short length = shortp[2];
            readDataFromPipe(inPipe, buffer, length);

            if (isScreenText)
                newScene->screenText.push_back(
                    ScreenText(string((char*)buffer, length)));
            else
                newScene->sceneText.push_back(
                    RenderedText(position, scale, color, string((char*)buffer, length), faceCamera));
            break;
        }

        case AddCoords: {
            readDataFromPipe(inPipe, buffer, 12*sizeof(float));
            fRotation rotation;
            rotation.setRotationToBodyFixedXYZ(fVec3(floatBuffer[0], 
                                                     floatBuffer[1], 
                                                     floatBuffer[2]));
            fVec3 position(floatBuffer[3], floatBuffer[4], floatBuffer[5]);
            fVec3 axisLengths(floatBuffer[6], floatBuffer[7], floatBuffer[8]);
            fVec3 textScale = fVec3(0.2f*min(axisLengths));
            float lineThickness = 1;
            fVec3 color = fVec3(floatBuffer[9], floatBuffer[10], floatBuffer[11]);
            int index;
            int numLines = (int)newScene->lines.size();
            for (index = 0; 
                 index < numLines && (color != newScene->lines[index].getColor() 
                                      || newScene->lines[index].getThickness() 
                                                              != lineThickness); 
                 index++)
                ;
            if (index == numLines)
                newScene->lines.push_back(RenderedLine(color, lineThickness));
            vector<float>& line = newScene->lines[index].getLines();
            fVec3 end = position+rotation*fVec3(axisLengths[0], 0, 0);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->sceneText.push_back(RenderedText(end, textScale, color, "X"));
            end = position+rotation*fVec3(0, axisLengths[1], 0);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->sceneText.push_back(RenderedText(end, textScale, color, "Y"));
            end = position+rotation*fVec3(0, 0, axisLengths[2]);
            line.push_back(position[0]);
            line.push_back(position[1]);
            line.push_back(position[2]);
            line.push_back(end[0]);
            line.push_back(end[1]);
            line.push_back(end[2]);
            newScene->sceneText.push_back(RenderedText(end, textScale, color, "Z"));
            break;
        }

        // Define a new mesh that will be assigned the next available mesh
        // index. It will be cached here and then can be referenced in this
        // scene and others by using it mesh index.
        case DefineMesh: {
			printf("=-=-=-=-= DefineMesh: =-=-=-=");
            readDataFromPipe(inPipe, buffer, 2*sizeof(short));
							 /*
            PendingMesh* mesh = new PendingMesh(); // assigns next mesh index
            int numVertices = shortBuffer[0];
            int numFaces = shortBuffer[1];
            mesh->vertices.resize(3*numVertices, 0);
            mesh->normals.resize(3*numVertices);
            mesh->faces.resize(3*numFaces);
            readData((unsigned char*)&mesh->vertices[0], (int)(mesh->vertices.size()*sizeof(float)));
            readData((unsigned char*)&mesh->faces[0], (int)(mesh->faces.size()*sizeof(short)));

            // Compute normal vectors for the mesh.

            vector<fVec3> normals(numVertices, fVec3(0));
            for (int i = 0; i < numFaces; i++) {
                int v1 = mesh->faces[3*i];
                int v2 = mesh->faces[3*i+1];
                int v3 = mesh->faces[3*i+2];
                fVec3 vert1(mesh->vertices[3*v1], mesh->vertices[3*v1+1], mesh->vertices[3*v1+2]);
                fVec3 vert2(mesh->vertices[3*v2], mesh->vertices[3*v2+1], mesh->vertices[3*v2+2]);
                fVec3 vert3(mesh->vertices[3*v3], mesh->vertices[3*v3+1], mesh->vertices[3*v3+2]);
                fVec3 norm = (vert2-vert1)%(vert3-vert1);
                float length = norm.norm();
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

            // A real mesh will be generated from this the next
            // time the scene is redrawn.
//            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE --------
//            pendingCommands.insert(pendingCommands.begin(), mesh);
//            pthread_mutex_unlock(&sceneLock);   //------ UNLOCK SCENE --------
            */
            break;
        }

        default: {
				 }
//            SimTK_ASSERT_ALWAYS(false, "Unexpected scene data sent to visualizer");
        }
    }

    return newScene;
}

//	Input listener
void VisualizerBase::listenForInput(int inPipe, int outPipe)
{
    unsigned char buffer[256];
    float* floatBuffer = (float*) buffer;
    int* intBuffer = (int*) buffer;
    unsigned short* shortBuffer = (unsigned short*) buffer;

    try
    { 
		while (true) 
		{
//        bool issuedActiveRedisplay = false;
        // Read commands from the simulator.
        readDataFromPipe(inPipe, buffer, 1);

        switch (buffer[0]) 
		{
        case DefineMenu: {
			readDataFromPipe(inPipe, buffer, sizeof(short));
            int titleLength = shortBuffer[0];
			std::vector<char> titleBuffer(titleLength);
            readDataFromPipe(inPipe, (unsigned char*)&titleBuffer[0], titleLength);
			std::string title(&titleBuffer[0], titleLength);
            readDataFromPipe(inPipe, buffer, sizeof(int));
            const int menuId = intBuffer[0];
            readDataFromPipe(inPipe, buffer, sizeof(short));
            int numItems = shortBuffer[0];
			std::vector<std::pair<std::string, int> > items(numItems);
            for (int index = 0; index < numItems; index++) {
                readDataFromPipe(inPipe, buffer, 2*sizeof(int));
                items[index].second = intBuffer[0];
				std::vector<char> textBuffer(intBuffer[1]);
                readDataFromPipe(inPipe, (unsigned char*)&textBuffer[0], intBuffer[1]);
                items[index].first = std::string(&textBuffer[0], intBuffer[1]);
            }
			std::cout << "===========DefineMenu was called===========\n";
			break;
        }
        case DefineSlider: {
			readDataFromPipe(inPipe, buffer, sizeof(short));
            int titleLength = shortBuffer[0];
			std::vector<char> titleBuffer(titleLength);
            readDataFromPipe(inPipe, (unsigned char*)&titleBuffer[0], titleLength);
			std::string title(&titleBuffer[0], titleLength);
            readDataFromPipe(inPipe, buffer, sizeof(int)+3*sizeof(float));
			std::cout << "===========DefineSlider was called==========\n";
			break;
        }
        case SetSliderValue: {
			std::cout << "============SetSliderValue was called=============\n";
			break;
        }
        case SetSliderRange: {
			std::cout << "===========SetSliderRange was called==============\n";
			break;
        }
        case SetCamera: {
			std::cout << "===========SetCamera was called================\n";
			break;
        }
        case ZoomCamera: {
			std::cout << "===========ZoomCamera was called===============\n";
			break;
        }
        case LookAt: {
			std::cout << "===========LookAt was called==============\n";
			break;
        }
        case SetFieldOfView: {
			std::cout << "===========SetFieldOfView was called============\n";
			break;
        }
        case SetClipPlanes: {
			std::cout << "===========SetClipPlanes was called=============\n";
			break;
        }
        case SetSystemUpDirection: {
			readDataFromPipe(inPipe, buffer, 2);
			std::cout << "===========SetSystemUpDirection was called============\n";
			break;
        }
        case SetGroundHeight: {
			std::cout << "===========SetGroundHeight was called===============\n";
			break;
        }
        case SetWindowTitle: {
			std::cout << "===========SetWindowTitle was called==============\n";
			break;
        }
        case SetMaxFrameRate: {
			readDataFromPipe(inPipe, buffer, sizeof(float));
			std::cout << "===========SetMaxFrameRate was called=============\n";
			break;
        }
        case SetBackgroundColor: {
			readDataFromPipe(inPipe, buffer, 3*sizeof(float));
			std::cout << "===========SetBackgroundColor was called============\n";
			break;
        }
        case SetShowShadows: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowShadows was called============\n";
			break;
        }
        case SetBackgroundType: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetBackgroundType was called===========\n";
			break;
        }
        case SetShowFrameRate: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowFrameRate was called===========\n";
			break;
        }
        case SetShowSimTime: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowSimTime was called============\n";
			break;
        }
        case SetShowFrameNumber: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowFrameNumber was called==========\n";
			break;
        }
        case StartOfScene: {
			printf( "===========StartOfScene was called===========\n");

            Scene* newScene = readNewScene(inPipe);
//            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            if (scene != NULL) {
                while (!scene->sceneHasBeenDrawn) {
                    // -------- WAIT FOR CONDITION --------
                    // Give up the lock and wait for notice.
//                    pthread_cond_wait(&sceneHasBeenDrawn,
//                                        &sceneLock);
                    // Now we hold the lock again, but this might
                    // be a spurious wakeup so check again.
                }
                // Previous scene has been drawn.
                delete scene; scene = 0;
            }
            // Swap in the new scene.
            scene = newScene;
//            saveNextFrameToMovie = savingMovie;
//            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
//            forceActiveRedisplay();             //------- ACTIVE REDISPLAY ---
//            issuedActiveRedisplay = true;
            break;
        }

        default:
            SimTK_ERRCHK1_ALWAYS(!"unrecognized command", "listenForInput()",
                "Unexpected command %u received from visualizer. Can't continue.",
                (unsigned)buffer[0]);
        }

        // Do this after every received command.
        // if (!issuedActiveRedisplay)
        //    requestPassiveRedisplay();         //------- PASSIVE REDISPLAY --
    	}
  } catch (const std::exception& e) {
        std::cout << "VisualizerBase listenerThread: unrecoverable error:\n";
        std::cout << e.what() << std::endl;
    }
}

void VisualizerBase::dumpMessage() 
{
    printf("\n\n=================== ABOUT SIMBODY VISUALIZER ===================\n");
    printf("Simbody(tm) %s VisualizerGUI (protocol rev. %u)\n",
        simbodyVersionStr.c_str(), ProtocolVersion);
    printf("\nName of invoking executable: %s\n", simulatorExecutableName.c_str());
    printf("================================================================\n\n");
}

