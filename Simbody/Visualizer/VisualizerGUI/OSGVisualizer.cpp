#include "OSGVisualizer.h"

using namespace osg;

OSGVisualizer::OSGVisualizer()
{
    root = new osg::Group();
}

osg::Node* OSGVisualizer::createScene()
{
	sceneGeode = new Geode();

	return sceneGeode;

}

void OSGVisualizer::drawLine(RenderedLine& line)
{

	osg::Geometry* lineGeom = new osg::Geometry();
	osg::Vec3Array* vertices = new osg::Vec3Array();

	for(int j = 0; j < line.getLines().size(); j+=3)
	{
		vertices->push_back( osg::Vec3(line.getLines().at(j),
										line.getLines().at(j+1),
										line.getLines().at(j+2)
										)
										);
	}

	lineGeom->setVertexArray(vertices);
	

	osg::Vec4Array * colours = new osg::Vec4Array;
	colours->push_back(osg::Vec4(line.getColor()[0], 
				line.getColor()[1],
				line.getColor()[2],
				0.0f
				)
				);
	lineGeom->setColorArray(colours);

	lineGeom->setColorBinding(osg::Geometry::BIND_OVERALL);

	// set the normal in the same way color.
	osg::Vec3Array* normals = new osg::Vec3Array;
	normals->push_back(osg::Vec3(0.0f,-1.0f,0.0f));
	lineGeom->setNormalArray(normals);
	lineGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);

	lineGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,
									line.getLines().size()));

	// Add Line to Geode
	sceneGeode->addDrawable( lineGeom );
}

void OSGVisualizer::renderScene()
{
//	this->root->removeChildren(0, this->root->getNumChildren());

	pthread_mutex_lock(&sceneLock);
	if( scene != NULL) 
	{
		for (int i = 0; i < scene->lines.size(); i++)
		{
			drawLine(scene->lines[i]);
		}
	}

	pthread_cond_signal(&sceneHasBeenDrawn);
	pthread_mutex_unlock(&sceneLock);
}

void OSGVisualizer::go()
{
	root->addChild(createScene());

	root->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

	viewer.setUpViewInWindow(100, 100, 800, 600);
    viewer.setSceneData( root );

    viewer.setCameraManipulator(new osgGA::TrackballManipulator());

	viewer.realize();
	while( !viewer.done() )
	{
		renderScene();
		viewer.frame();
	}

}
int main( int argc, char** argv) 
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

	OSGVisualizer app;

	// Create application object

	if(talkingToSimulator) 
	{
		app.shakeHandsWithSimulator(inPipe, outPipe);
		pthread_t thread;

		VisualizerBase::VisualizerCom vcom;
		vcom.arg = &app;
		vcom.inPipe = inPipe;
		vcom.outPipe = outPipe;
		int rc = pthread_create(&thread, NULL, VisualizerBase::listenForInputHelper, &vcom);
		//int rc = pthread_create(&thread, NULL, VisualizerBase::listenForInputHelper, &app);
	}

	app.createScene();
	app.go();

}
