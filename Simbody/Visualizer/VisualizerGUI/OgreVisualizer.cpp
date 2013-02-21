/*
 * OgreVisualizer.cpp
 * Taken from Ogre Wiki
*/
#include "OgreVisualizer.h"
#include "DynamicLines.h"
#include <OgreSceneNode.h>
#include <pthread.h>
#include <sstream>

#include "MovableText.h"

//-------------------------------------------------------------------------------------
OgreVisualizer::OgreVisualizer() : lineCount(0)
{
}
//-------------------------------------------------------------------------------------
OgreVisualizer::~OgreVisualizer(void)
{
}

//-------------------------------------------------------------------------------------
void OgreVisualizer::createScene(void)
{
	//
//	mSceneMgr->setSkyDome(true, "Examples/CloudySky", 5, 8);	
	/*
	Ogre::Plane skyPlane;
	skyPlane.d = 1000;
	skyPlane.normal = Ogre::Vector3::NEGATIVE_UNIT_Y;

	mSceneMgr->setSkyPlane(true, skyPlane, "Examples/CloudySky", 1500, 40, true, 1.5f, 150, 150);
	*/
//   mSceneMgr->setAmbientLight(Ogre::ColourValue(0, 0, 0));
    mSceneMgr->setShadowTechnique(Ogre::SHADOWTYPE_STENCIL_ADDITIVE);

//    Ogre::Entity* entNinja = mSceneMgr->createEntity("Ninja", "ninja.mesh");
//    entNinja->setCastShadows(true);
//    mSceneMgr->getRootSceneNode()->createChildSceneNode()->attachObject(entNinja);

    Ogre::Plane plane(Ogre::Vector3::UNIT_Y, 0);
 
    Ogre::MeshManager::getSingleton().createPlane("ground", Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME,
        plane, 90000, 90000, 20, 20, true, 1, 2, 2, Ogre::Vector3::UNIT_Z);
 
    Ogre::Entity* entGround = mSceneMgr->createEntity("GroundEntity", "ground");

    mSceneMgr->getRootSceneNode()->createChildSceneNode()->attachObject(entGround);

	Ogre::ResourcePtr groundMat = Ogre::MaterialManager::getSingleton().create("ground", "General");
	Ogre::MaterialPtr mat = Ogre::MaterialPtr(groundMat);
	mat.get()->setAmbient(0.95, 0.64, 0.23);
	entGround->setMaterialName("ground");

    entGround->setCastShadows(false);
 
    Ogre::Light* pointLight = mSceneMgr->createLight("pointLight");
    pointLight->setType(Ogre::Light::LT_POINT);
    pointLight->setPosition(Ogre::Vector3(0, 150, 250));
 
    pointLight->setDiffuseColour(1.0, 0.0, 0.0);
    pointLight->setSpecularColour(1.0, 0.0, 0.0);

    Ogre::Light* directionalLight = mSceneMgr->createLight("directionalLight");
    directionalLight->setType(Ogre::Light::LT_DIRECTIONAL);
    directionalLight->setDiffuseColour(Ogre::ColourValue(.25, .25, 0));
    directionalLight->setSpecularColour(Ogre::ColourValue(.25, .25, 0));
 
    directionalLight->setDirection(Ogre::Vector3( 0, -1, 1 )); 
 
	/*
    Ogre::Light* spotLight = mSceneMgr->createLight("spotLight");
    spotLight->setType(Ogre::Light::LT_SPOTLIGHT);
    spotLight->setDiffuseColour(0, 0, 1.0);
    spotLight->setSpecularColour(0, 0, 1.0);
 
    spotLight->setDirection(-1, -1, 0);
    spotLight->setPosition(Ogre::Vector3(300, 300, 0));
 
    spotLight->setSpotlightRange(Ogre::Degree(35), Ogre::Degree(50));
	*/
	linesNode = mSceneMgr->getRootSceneNode()->createChildSceneNode();
}

void OgreVisualizer::drawLine(RenderedLine& line, std::string name)
{
	Ogre::SceneNode *sceneNode = NULL;
  	Ogre::ManualObject *obj = NULL;
  	bool attached = false;

	float r,g,b;
	r = line.getColor()[0];
	g = line.getColor()[1];
	b = line.getColor()[2];

	Ogre::MaterialPtr matPtr = Ogre::MaterialManager::getSingleton().create(name, "General");
	matPtr->setReceiveShadows(false);
	matPtr->getTechnique(0)->setLightingEnabled(true);
//	matPtr->getTechnique(0)->getPass(0)->setDiffuse(r,g,b,0);
	matPtr->getTechnique(0)->getPass(0)->setAmbient(r,g,b);
//	matPtr->getTechnique(0)->getPass(0)->setSelfIllumination(0,0,1);

  	if (this->mSceneMgr->hasManualObject(name))
  	{
    	sceneNode = this->mSceneMgr->getSceneNode(name);
   		obj = this->mSceneMgr->getManualObject(name);
    	attached = true;
  	}
  	else
  	{
    	sceneNode = this->mSceneMgr->getRootSceneNode()->createChildSceneNode(name);
    	obj = this->mSceneMgr->createManualObject(name);
  	}

  	sceneNode->setVisible(true);
  	obj->setVisible(true);

  	obj->begin(name, Ogre::RenderOperation::OT_LINE_LIST);
	
	for (int j = 0; j < line.getLines().size(); j+=3 )
	{
		obj->position(line.getLines().at(j),
						line.getLines().at(j+1),
						line.getLines().at(j+2));
	}

  	obj->end();

  	if (!attached)
	{
    	sceneNode->attachObject(obj);
	}
}

void OgreVisualizer::drawText(RenderedText& text, const std::string& name)
{
	Ogre::MovableText* msg = new Ogre::MovableText(name, text.getText());
//	msg->setGlobalTranslation(Ogre::Vector3(text.getPosition()[0], text.getPosition()[1], text.getPosition()[2]));
	msg->setLocalTranslation(Ogre::Vector3(text.getPosition()[0], text.getPosition()[1], text.getPosition()[2]));
	msg->setColor(Ogre::ColourValue(text.getColor()[0], text.getColor()[1], text.getColor()[2]));
	msg->setCharacterHeight(0.1);
	msg->setSpaceWidth(0.1);

	Ogre::SceneNode *sceneNode = this->mSceneMgr->getRootSceneNode()->createChildSceneNode("txt" + name);

	sceneNode->setVisible(true);
	sceneNode->attachObject(msg);
}

void OgreVisualizer::renderScene()
{

	pthread_mutex_lock(&sceneLock);

    if (scene != NULL) {

		std::string str;
		ostringstream convert;

        for (int i = 0; i < (int) scene->lines.size(); i++)
		{
		convert << "line" << lineCount;
		str = convert.str();
			drawLine(scene->lines[i], str);
			lineCount++;
        }
        for (int i = 0; i < (int) scene->sceneText.size(); i++)
		{
		convert << "line" << lineCount;
		str = convert.str();

            drawText(scene->sceneText[i], str);
			lineCount++;
		}
        for (int i = 0; i < (int) scene->drawnMeshes.size(); i++)
            scene->drawnMeshes[i].draw();
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
		{
            scene->solidMeshes[i].draw();
		}
//        vector<pair<float, int> > order(scene->transparentMeshes.size());
//        for (int i = 0; i < (int) order.size(); i++)
//            order[i] = make_pair((float)(~X_GC.R()*scene->transparentMeshes[i].getTransform().p())[2], i);
//        sort(order.begin(), order.end());
//        for (int i = 0; i < (int) order.size(); i++)
//            scene->transparentMeshes[order[i].second].draw();

        scene->sceneHasBeenDrawn = true;
    }

	pthread_cond_signal(&sceneHasBeenDrawn);  
	pthread_mutex_unlock(&sceneLock);

}

void OgreVisualizer::go()
{

	if(!setup())
	{
		return;
	}

	createScene();

	while(true)
	{
		renderScene();
		// Do stuff
		mRoot->renderOneFrame();
	}
}

int main(int argc, char *argv[])
{

	bool talkingToSimulator = false;
	int inPipe, outPipe;

	if (argc >= 3) {
		std::stringstream(argv[1]) >> inPipe;
		std::stringstream(argv[2]) >> outPipe;
        talkingToSimulator = true; // presumably those were the pipes
  	} else {
        printf("\n**** VISUALIZER HAS NO SIMULATOR TO TALK TO ****\n");
        printf("The Simbody VisualizerGUI was invoked directly with no simulator\n");
        printf("process to talk to. Will attempt to bring up the display anyway\n");
        printf("in case you want to look at the About message.\n");
        printf("The VisualizerGUI is intended to be invoked programmatically.\n");
    }

	// Create application object
	OgreVisualizer app;

	std::cout << "-=-=-=-[Initializing] inPipe: " << inPipe << std::endl;
	if(talkingToSimulator) 
	{
		app.shakeHandsWithSimulator(inPipe, outPipe);
		pthread_t thread;

		VisualizerBase::VisualizerCom vcom;
		vcom.arg = &app;
		vcom.inPipe = inPipe;
		vcom.outPipe = outPipe;
		int rc = pthread_create(&thread, NULL, VisualizerBase::listenForInputHelper, &vcom);
	}
//	app.dumpMessage();

	try {
		app.go();

	} catch( Ogre::Exception& e ) {
		std::cerr << "An exception has occured: " <<
			e.getFullDescription().c_str() << std::endl;
	}

	return 0;

}
