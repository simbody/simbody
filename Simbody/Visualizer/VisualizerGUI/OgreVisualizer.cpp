/*
 * OgreVisualizer.cpp
 * Taken from Ogre Wiki
*/
#include "OgreVisualizer.h"
#include "DynamicLines.h"
#include <OgreSceneNode.h>
#include <pthread.h>
#include <sstream>

//-------------------------------------------------------------------------------------
OgreVisualizer::OgreVisualizer()
{
//	CreateContext();
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
	Ogre::Plane skyPlane;
	skyPlane.d = 1000;
	skyPlane.normal = Ogre::Vector3::NEGATIVE_UNIT_Y;

	mSceneMgr->setSkyPlane(true, skyPlane, "Examples/CloudySky", 1500, 40, true, 1.5f, 150, 150);

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
 
    Ogre::Light* spotLight = mSceneMgr->createLight("spotLight");
    spotLight->setType(Ogre::Light::LT_SPOTLIGHT);
    spotLight->setDiffuseColour(0, 0, 1.0);
    spotLight->setSpecularColour(0, 0, 1.0);
 
    spotLight->setDirection(-1, -1, 0);
    spotLight->setPosition(Ogre::Vector3(300, 300, 0));
 
    spotLight->setSpotlightRange(Ogre::Degree(35), Ogre::Degree(50));
}

void OgreVisualizer::renderScene()
{

	Ogre::ResourcePtr lineMat = Ogre::MaterialManager::getSingleton().create("line", "General");
	Ogre::MaterialPtr mat = Ogre::MaterialPtr(lineMat);

    if (scene != NULL) {

        for (int i = 0; i < (int) scene->lines.size(); i++)
		{
			DynamicLines *line = new DynamicLines(Ogre::RenderOperation::OT_LINE_LIST);

			for (int j = 0; j < scene->lines[i].getLines().size(); j+=3 )
			{
//				std::cout << "Lines points: " << scene->lines[i].getLines().at(j) << std::endl;
				Ogre::Vector3 vec(scene->lines[i].getLines().at(j),
						scene->lines[i].getLines().at(j+1),
						scene->lines[i].getLines().at(j+2));
				line->addPoint(vec);
			}
			
			mat.get()->setAmbient(scene->lines[i].getColor()[0], 
			scene->lines[i].getColor()[1], 
			scene->lines[i].getColor()[2]);

			line->setMaterial("line");
			
			line->update();
			Ogre::SceneNode *linesNode = mSceneMgr->getRootSceneNode()->createChildSceneNode();
			linesNode->attachObject(line);

//        glLineWidth(2);
        }
        for (int i = 0; i < (int) scene->sceneText.size(); i++)
            scene->sceneText[i].draw();
//        glLineWidth(1);
//        glEnableClientState(GL_NORMAL_ARRAY);
        for (int i = 0; i < (int) scene->drawnMeshes.size(); i++)
            scene->drawnMeshes[i].draw();
//        glEnable(GL_LIGHTING);
        for (int i = 0; i < (int) scene->solidMeshes.size(); i++)
            scene->solidMeshes[i].draw();

//        vector<pair<float, int> > order(scene->transparentMeshes.size());
//        for (int i = 0; i < (int) order.size(); i++)
//            order[i] = make_pair((float)(~X_GC.R()*scene->transparentMeshes[i].getTransform().p())[2], i);
//        sort(order.begin(), order.end());
//        for (int i = 0; i < (int) order.size(); i++)
//            scene->transparentMeshes[order[i].second].draw();


        scene->sceneHasBeenDrawn = true;
    }

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
