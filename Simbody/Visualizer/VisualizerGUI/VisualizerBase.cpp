#include "VisualizerBase.h"
#include <vector>
#include <utility>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/glx.h>

VisualizerBase::VisualizerBase()
{
	dummyDisplay = NULL;
}

VisualizerBase::~VisualizerBase()
{
}

bool VisualizerBase::CreateContext()
{
  bool result = true;

  try
  {
    this->dummyDisplay = XOpenDisplay(0);
    if (!this->dummyDisplay)
    {
		std::cerr << "Can't open display: " << XDisplayName(0) << "\n";
      return false;
    }

    int screen = DefaultScreen(this->dummyDisplay);

    int attribList[] = {GLX_RGBA, GLX_DOUBLEBUFFER, GLX_DEPTH_SIZE, 16,
      GLX_STENCIL_SIZE, 8, None };

    XVisualInfo *dummyVisual = glXChooseVisual(
        static_cast<Display*>(this->dummyDisplay),
        screen, static_cast<int *>(attribList));

    if (!dummyVisual)
    {
		std::cerr << "Unable to create glx visual\n";
      return false;
    }

    this->dummyWindowId = XCreateSimpleWindow(
        static_cast<Display*>(this->dummyDisplay),
        RootWindow(static_cast<Display*>(this->dummyDisplay), screen),
        0, 0, 1, 1, 0, 0, 0);

    this->dummyContext = glXCreateContext(
        static_cast<Display*>(this->dummyDisplay),
        dummyVisual, NULL, 1);

    if (!this->dummyContext)
    {
		std::cerr << "Unable to create glx context\n";
      return false;
    }

    glXMakeCurrent(static_cast<Display*>(this->dummyDisplay),
        this->dummyWindowId, static_cast<GLXContext>(this->dummyContext));
  }
  catch(...)
  {
    result = false;
  }


  return result;
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
        }
        case DefineSlider: {
			readDataFromPipe(inPipe, buffer, sizeof(short));
            int titleLength = shortBuffer[0];
			std::vector<char> titleBuffer(titleLength);
            readDataFromPipe(inPipe, (unsigned char*)&titleBuffer[0], titleLength);
			std::string title(&titleBuffer[0], titleLength);
            readDataFromPipe(inPipe, buffer, sizeof(int)+3*sizeof(float));
			std::cout << "===========DefineSlider was called==========\n";
        }
        case SetSliderValue: {
			std::cout << "============SetSliderValue was called=============\n";
        }
        case SetSliderRange: {
			std::cout << "===========SetSliderRange was called==============\n";
        }
        case SetCamera: {
			std::cout << "===========SetCamera was called================\n";
        }
        case ZoomCamera: {
			std::cout << "===========ZoomCamera was called===============\n";
        }
        case LookAt: {
			std::cout << "===========LookAt was called==============\n";
        }
        case SetFieldOfView: {
			std::cout << "===========SetFieldOfView was called============\n";
        }
        case SetClipPlanes: {
			std::cout << "===========SetClipPlanes was called=============\n";
        }
        case SetSystemUpDirection: {
			std::cout << "===========SetSystemUpDirection was called============\n";
        }
        case SetGroundHeight: {
			std::cout << "===========SetGroundHeight was called===============\n";
        }
        case SetWindowTitle: {
			std::cout << "===========SetWindowTitle was called==============\n";
        }
        case SetMaxFrameRate: {
			readDataFromPipe(inPipe, buffer, sizeof(float));
			std::cout << "===========SetMaxFrameRate was called=============\n";
        }
        case SetBackgroundColor: {
			readDataFromPipe(inPipe, buffer, 3*sizeof(float));
			std::cout << "===========SetBackgroundColor was called============\n";
        }
        case SetShowShadows: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowShadows was called============\n";
        }
        case SetBackgroundType: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetBackgroundType was called===========\n";
        }
        case SetShowFrameRate: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowFrameRate was called===========\n";
        }
        case SetShowSimTime: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowSimTime was called============\n";
        }
        case SetShowFrameNumber: {
            readDataFromPipe(inPipe, buffer, sizeof(short));
			std::cout << "===========SetShowFrameNumber was called==========\n";
        }
        case StartOfScene: {
			printf( "===========StartOfScene was called===========\n");
			
			this->drawGroundAndSky();
			this->startRendering();
			/*
            Scene* newScene = readNewScene();
            pthread_mutex_lock(&sceneLock);     //------- LOCK SCENE ---------
            if (scene != NULL) {
                while (!scene->sceneHasBeenDrawn) {
                    // -------- WAIT FOR CONDITION --------
                    // Give up the lock and wait for notice.
                    pthread_cond_wait(&sceneHasBeenDrawn,
                                        &sceneLock);
                    // Now we hold the lock again, but this might
                    // be a spurious wakeup so check again.
                }
                // Previous scene has been drawn.
                delete scene; scene = 0;
            }
            // Swap in the new scene.
            scene = newScene;
            saveNextFrameToMovie = savingMovie;
            pthread_mutex_unlock(&sceneLock);   //------- UNLOCK SCENE -------
            forceActiveRedisplay();             //------- ACTIVE REDISPLAY ---
            issuedActiveRedisplay = true;
            break;
			*/
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

