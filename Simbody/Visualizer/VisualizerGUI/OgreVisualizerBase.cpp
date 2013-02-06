/*
 * OgreVisualizerBase.cpp
 * Taken from Ogre wiki
*/
#include "OgreVisualizerBase.h"

//-------------------------------------------------------------------------------------
OgreVisualizerBase::OgreVisualizerBase(void)
    : mRoot(0),
    mCamera(0),
    mSceneMgr(0),
    mWindow(0),
    mTrayMgr(0),
    mCameraMan(0),
    mDetailsPanel(0),
    mCursorWasVisible(false),
    mShutDown(false),
    mInputManager(0),
    mMouse(0),
    mKeyboard(0)
{
}

//-------------------------------------------------------------------------------------
OgreVisualizerBase::~OgreVisualizerBase(void)
{
    if (mTrayMgr) delete mTrayMgr;
    if (mCameraMan) delete mCameraMan;

    //Remove ourself as a Window listener
    Ogre::WindowEventUtilities::removeWindowEventListener(mWindow, this);
    windowClosed(mWindow);
    delete mRoot;
}

//-------------------------------------------------------------------------------------
bool OgreVisualizerBase::configure(void)
{
	Ogre::RenderSystem* rs = NULL;
  const Ogre::RenderSystemList *rsList;

  	// Set parameters of render system (window size, etc.)
#if OGRE_VERSION_MAJOR == 1 && OGRE_VERSION_MINOR == 6
    rsList = mRoot->getAvailableRenderers();
#else
    rsList = &(mRoot->getAvailableRenderers());
#endif

  	int c = 0;
  	rs = NULL;

	do
	{
		if (c == static_cast<int>(rsList->size()))
		break;

		rs = rsList->at(c);
		c++;
	}
	while (rs &&
			 rs->getName().compare("OpenGL Rendering Subsystem") != 0);

    if (rs == NULL)
    	throw("unable to find rendering system");
		
	  mRoot->getRenderSystemByName("OpenGL Rendering Subsystem");
      // adjust its options
      rs->setConfigOption("Full Screen", "No");
      // define it has the system renderer
      mRoot->setRenderSystem(rs);

//	  mRoot->showConfigDialog();
      // init root
      mRoot->initialise(false);
      // create the main window
      mWindow = mRoot->createRenderWindow("Ogre Visualizer",640,480,false);
	
	  return true;
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::chooseSceneManager(void)
{
    // Get the SceneManager, in this case a generic one
    mSceneMgr = mRoot->createSceneManager(Ogre::ST_GENERIC);
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::createCamera(void)
{
    // Create the camera
    mCamera = mSceneMgr->createCamera("PlayerCam");

    // Position it at 500 in Z direction
    mCamera->setPosition(Ogre::Vector3(0,0,500));
    // Look back along -Z
    mCamera->setNearClipDistance(5);

    mCameraMan = new OgreBites::SdkCameraMan(mCamera);   // create a default camera controller

	mCameraMan->setStyle(OgreBites::CS_ORBIT);
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::createFrameListener(void)
{
    mLoggerMan->logMessage("*** Initializing OIS ***");
    OIS::ParamList pl;
    size_t windowHnd = 0;
    std::ostringstream windowHndStr;

    mWindow->getCustomAttribute("WINDOW", &windowHnd);
    windowHndStr << windowHnd;
    pl.insert(std::make_pair(std::string("WINDOW"), windowHndStr.str()));

	pl.insert(std::make_pair(std::string("x11_mouse_grab"), std::string("false")));
	pl.insert(std::make_pair(std::string("x11_mouse_hide"), std::string("false")));
	pl.insert(std::make_pair(std::string("x11_keyboard_grab"), std::string("false")));
	pl.insert(std::make_pair(std::string("XAutoRepeatOn"), std::string("true")));

    mInputManager = OIS::InputManager::createInputSystem( pl );

    mKeyboard = static_cast<OIS::Keyboard*>(mInputManager->createInputObject( OIS::OISKeyboard, true ));
    mMouse = static_cast<OIS::Mouse*>(mInputManager->createInputObject( OIS::OISMouse, true ));

    mMouse->setEventCallback(this);
    mKeyboard->setEventCallback(this);

    //Set initial mouse clipping size
    windowResized(mWindow);

    //Register as a Window listener
    Ogre::WindowEventUtilities::addWindowEventListener(mWindow, this);

    mTrayMgr = new OgreBites::SdkTrayManager("InterfaceName", mWindow, mMouse, this);
    mTrayMgr->showFrameStats(OgreBites::TL_BOTTOMLEFT);
    mTrayMgr->showLogo(OgreBites::TL_BOTTOMRIGHT);
    mTrayMgr->hideCursor();

    // create a params panel for displaying sample details
    Ogre::StringVector items;
    items.push_back("cam.pX");
    items.push_back("cam.pY");
    items.push_back("cam.pZ");
    items.push_back("");
    items.push_back("cam.oW");
    items.push_back("cam.oX");
    items.push_back("cam.oY");
    items.push_back("cam.oZ");
    items.push_back("");
    items.push_back("Filtering");
    items.push_back("Poly Mode");

    mDetailsPanel = mTrayMgr->createParamsPanel(OgreBites::TL_NONE, "DetailsPanel", 200, items);
    mDetailsPanel->setParamValue(9, "Bilinear");
    mDetailsPanel->setParamValue(10, "Solid");
    mDetailsPanel->hide();

    mRoot->addFrameListener(this);
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::destroyScene(void)
{
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::createViewports(void)
{
    // Create one viewport, entire window
    Ogre::Viewport* vp = mWindow->addViewport(mCamera);
    vp->setBackgroundColour(Ogre::ColourValue(0,0,0));

    // Alter the camera aspect ratio to match the viewport
    mCamera->setAspectRatio(
        Ogre::Real(vp->getActualWidth()) / Ogre::Real(vp->getActualHeight()));
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::setupResources(void)
{
	      // get the resource manager
      Ogre::ResourceGroupManager &resgroup = Ogre::ResourceGroupManager::getSingleton();
	  resgroup.addResourceLocation("/home/kevinux/projects/Visualizer/tools/Ogre/ogre_src_v1-8-1/Samples/Media", "FileSystem");
	  resgroup.addResourceLocation("/home/kevinux/projects/Visualizer/tools/Ogre/ogre_src_v1-8-1/Samples/Media/models", "FileSystem");
	  resgroup.addResourceLocation("/home/kevinux/projects/Visualizer/tools/Ogre/ogre_src_v1-8-1/Samples/Media/materials/scripts", "FileSystem");
	  resgroup.addResourceLocation("/home/kevinux/projects/Visualizer/tools/Ogre/ogre_src_v1-8-1/Samples/Media/materials/textures", "FileSystem");
	  resgroup.addResourceLocation("/home/kevinux/projects/Visualizer/tools/Ogre/ogre_src_v1-8-1/Samples/Media/materials/programs", "FileSystem");
	  resgroup.addResourceLocation("/home/kevinux/projects/Visualizer/tools/Ogre/ogre_src_v1-8-1/Samples/Media/packs/SdkTrays.zip", "Zip", "General" );
}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::createResourceListener(void)
{

}
//-------------------------------------------------------------------------------------
void OgreVisualizerBase::loadResources(void)
{
    Ogre::ResourceGroupManager::getSingleton().initialiseAllResourceGroups();
}

void OgreVisualizerBase::loadPlugins()
{
std::list<std::string>::iterator iter;

std::list<std::string> ogrePaths;
ogrePaths.push_back("/usr/local/lib/OGRE");
ogrePaths.push_back("/usr/lib/OGRE");

  for (iter = ogrePaths.begin();
       iter!= ogrePaths.end(); ++iter)
  {
    std::string path(*iter);
	/*
    DIR *dir = opendir(path.c_str());

    if (dir == NULL)
    {
      continue;
    }
    closedir(dir);
*/
    std::vector<std::string> plugins;
    std::vector<std::string>::iterator piter;

	std::cout << "path: " << path << std::endl;

    plugins.push_back(path+"/RenderSystem_GL.so");
    plugins.push_back(path+"/Plugin_ParticleFX.so");
    plugins.push_back(path+"/Plugin_BSPSceneManager.so");
    plugins.push_back(path+"/Plugin_OctreeSceneManager.so");
    plugins.push_back(path+"/Plugin_CgProgramManager.so");

    for (piter = plugins.begin(); piter!= plugins.end(); ++piter)
    {
      try
      {
        // Load the plugin into OGRE
        mRoot->loadPlugin(*piter);
      }
      catch(Ogre::Exception &e)
      {
        if ((*piter).find("RenderSystem") != std::string::npos)
        {
          std::string description("Unable to load Ogre Plugin[");
          description.append(*piter);
          description.append("]...This won't end well.");
//          gzerr << description << "\n";
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------
void OgreVisualizerBase::go(void)
{
    if (!setup())
        return;

    mRoot->startRendering();
//	mRoot->renderOneFrame();

    // clean up
    destroyScene();
}
//-------------------------------------------------------------------------------------
bool OgreVisualizerBase::setup(void)
{

	mLoggerMan = new Ogre::LogManager();
	mLoggerMan->createLog("ogre.log", true, false, false);

    mRoot = new Ogre::Root();

	loadPlugins();
    setupResources();

    bool carryOn = configure();
    if (!carryOn) return false;

    chooseSceneManager();
    createCamera();
    createViewports();

    // Set default mipmap level (NB some APIs ignore this)
    Ogre::TextureManager::getSingleton().setDefaultNumMipmaps(5);

    // Create any resource listeners (for loading screens)
    createResourceListener();
    // Load resources
    loadResources();

	mWindow->setFullscreen(false, 800, 600);
    // Create the scene
//    createScene();

    createFrameListener();

    return true;
};
//-------------------------------------------------------------------------------------
bool OgreVisualizerBase::frameRenderingQueued(const Ogre::FrameEvent& evt)
{
    if(mWindow->isClosed())
        return false;

    if(mShutDown)
        return false;

    //Need to capture/update each device
    mKeyboard->capture();
    mMouse->capture();

    mTrayMgr->frameRenderingQueued(evt);

    if (!mTrayMgr->isDialogVisible())
    {
        mCameraMan->frameRenderingQueued(evt);   // if dialog isn't up, then update the camera
        if (mDetailsPanel->isVisible())   // if details panel is visible, then update its contents
        {
            mDetailsPanel->setParamValue(0, Ogre::StringConverter::toString(mCamera->getDerivedPosition().x));
            mDetailsPanel->setParamValue(1, Ogre::StringConverter::toString(mCamera->getDerivedPosition().y));
            mDetailsPanel->setParamValue(2, Ogre::StringConverter::toString(mCamera->getDerivedPosition().z));
            mDetailsPanel->setParamValue(4, Ogre::StringConverter::toString(mCamera->getDerivedOrientation().w));
            mDetailsPanel->setParamValue(5, Ogre::StringConverter::toString(mCamera->getDerivedOrientation().x));
            mDetailsPanel->setParamValue(6, Ogre::StringConverter::toString(mCamera->getDerivedOrientation().y));
            mDetailsPanel->setParamValue(7, Ogre::StringConverter::toString(mCamera->getDerivedOrientation().z));
        }
    }

    return true;
}
//-------------------------------------------------------------------------------------
bool OgreVisualizerBase::keyPressed( const OIS::KeyEvent &arg )
{
    if (mTrayMgr->isDialogVisible()) return true;   // don't process any more keys if dialog is up

    if (arg.key == OIS::KC_F)   // toggle visibility of advanced frame stats
    {
        mTrayMgr->toggleAdvancedFrameStats();
    }
    else if (arg.key == OIS::KC_G)   // toggle visibility of even rarer debugging details
    {
        if (mDetailsPanel->getTrayLocation() == OgreBites::TL_NONE)
        {
            mTrayMgr->moveWidgetToTray(mDetailsPanel, OgreBites::TL_TOPRIGHT, 0);
            mDetailsPanel->show();
        }
        else
        {
            mTrayMgr->removeWidgetFromTray(mDetailsPanel);
            mDetailsPanel->hide();
        }
    }
    else if (arg.key == OIS::KC_T)   // cycle polygon rendering mode
    {
        Ogre::String newVal;
        Ogre::TextureFilterOptions tfo;
        unsigned int aniso;

        switch (mDetailsPanel->getParamValue(9).asUTF8()[0])
        {
        case 'B':
            newVal = "Trilinear";
            tfo = Ogre::TFO_TRILINEAR;
            aniso = 1;
            break;
        case 'T':
            newVal = "Anisotropic";
            tfo = Ogre::TFO_ANISOTROPIC;
            aniso = 8;
            break;
        case 'A':
            newVal = "None";
            tfo = Ogre::TFO_NONE;
            aniso = 1;
            break;
        default:
            newVal = "Bilinear";
            tfo = Ogre::TFO_BILINEAR;
            aniso = 1;
        }

        Ogre::MaterialManager::getSingleton().setDefaultTextureFiltering(tfo);
        Ogre::MaterialManager::getSingleton().setDefaultAnisotropy(aniso);
        mDetailsPanel->setParamValue(9, newVal);
    }
    else if (arg.key == OIS::KC_R)   // cycle polygon rendering mode
    {
        Ogre::String newVal;
        Ogre::PolygonMode pm;

        switch (mCamera->getPolygonMode())
        {
        case Ogre::PM_SOLID:
            newVal = "Wireframe";
            pm = Ogre::PM_WIREFRAME;
            break;
        case Ogre::PM_WIREFRAME:
            newVal = "Points";
            pm = Ogre::PM_POINTS;
            break;
        default:
            newVal = "Solid";
            pm = Ogre::PM_SOLID;
        }

        mCamera->setPolygonMode(pm);
        mDetailsPanel->setParamValue(10, newVal);
    }
    else if(arg.key == OIS::KC_F5)   // refresh all textures
    {
        Ogre::TextureManager::getSingleton().reloadAll();
    }
    else if (arg.key == OIS::KC_SYSRQ)   // take a screenshot
    {
        mWindow->writeContentsToTimestampedFile("screenshot", ".jpg");
    }
    else if (arg.key == OIS::KC_ESCAPE)
    {
        mShutDown = true;
    }

    mCameraMan->injectKeyDown(arg);
    return true;
}

bool OgreVisualizerBase::keyReleased( const OIS::KeyEvent &arg )
{
    mCameraMan->injectKeyUp(arg);
    return true;
}

bool OgreVisualizerBase::mouseMoved( const OIS::MouseEvent &arg )
{
    if (mTrayMgr->injectMouseMove(arg)) return true;
    mCameraMan->injectMouseMove(arg);
   return true;
}

bool OgreVisualizerBase::mousePressed( const OIS::MouseEvent &arg, OIS::MouseButtonID id )
{
    if (mTrayMgr->injectMouseDown(arg, id)) return true;
    mCameraMan->injectMouseDown(arg, id);
    return true;
}

bool OgreVisualizerBase::mouseReleased( const OIS::MouseEvent &arg, OIS::MouseButtonID id )
{
    if (mTrayMgr->injectMouseUp(arg, id)) return true;
    mCameraMan->injectMouseUp(arg, id);
    return true;
}

//Adjust mouse clipping area
void OgreVisualizerBase::windowResized(Ogre::RenderWindow* rw)
{
    unsigned int width, height, depth;
    int left, top;
    rw->getMetrics(width, height, depth, left, top);

    const OIS::MouseState &ms = mMouse->getMouseState();
    ms.width = width;
    ms.height = height;
}

//Unattach OIS before window shutdown (very important under Linux)
void OgreVisualizerBase::windowClosed(Ogre::RenderWindow* rw)
{
    //Only close for window that created OIS (the main window in these demos)
    if( rw == mWindow )
    {
        if( mInputManager )
        {
            mInputManager->destroyInputObject( mMouse );
            mInputManager->destroyInputObject( mKeyboard );

            OIS::InputManager::destroyInputSystem(mInputManager);
            mInputManager = 0;
        }
    }
}
