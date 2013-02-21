#include "OSGVisualizer.h"

using namespace osg;

osg::Node* OSGVisualizer::makeSky()
{
    int i, j;
    float lev[] = { -5, -1.0, 1.0, 15.0, 30.0, 60.0, 90.0  };
    float cc[][4] =
    {
        { 0.0, 0.0, 0.15 },
        { 0.0, 0.0, 0.15 },
        { 0.4, 0.4, 0.7 },
        { 0.2, 0.2, 0.6 },
        { 0.1, 0.1, 0.6 },
        { 0.1, 0.1, 0.6 },
        { 0.1, 0.1, 0.6 },
    };
    float x, y, z;
    float alpha, theta;
    float radius = 20.0f;
    int nlev = sizeof( lev )/sizeof(float);

    Geometry *geom = new Geometry;

    Vec3Array& coords = *(new Vec3Array(19*nlev));
    Vec4Array& colors = *(new Vec4Array(19*nlev));
    Vec2Array& tcoords = *(new Vec2Array(19*nlev));
    
    
    int ci = 0;

    for( i = 0; i < nlev; i++ )
    {
        for( j = 0; j <= 18; j++ )
        {
            alpha = osg::DegreesToRadians(lev[i]);
            theta = osg::DegreesToRadians((float)(j*20));

            x = radius * cosf( alpha ) * cosf( theta );
            y = radius * cosf( alpha ) * -sinf( theta );
            z = radius * sinf( alpha );

            coords[ci][0] = x;
            coords[ci][1] = y;
            coords[ci][2] = z;

            colors[ci][0] = cc[i][0];
            colors[ci][1] = cc[i][1];
            colors[ci][2] = cc[i][2];
            colors[ci][3] = 1.0;

            tcoords[ci][0] = (float)j/18.0;
            tcoords[ci][1] = (float)i/(float)(nlev-1);

            ci++;
        }


    }

    for( i = 0; i < nlev-1; i++ )
    {
        DrawElementsUShort* drawElements = new DrawElementsUShort(PrimitiveSet::TRIANGLE_STRIP);
        drawElements->reserve(38);

        for( j = 0; j <= 18; j++ )
        {
            drawElements->push_back((i+1)*19+j);
            drawElements->push_back((i+0)*19+j);
        }

        geom->addPrimitiveSet(drawElements);
    }
    
    geom->setVertexArray( &coords );
    geom->setTexCoordArray( 0, &tcoords );

    geom->setColorArray( &colors );
    geom->setColorBinding( Geometry::BIND_PER_VERTEX );


    Texture2D *tex = new Texture2D;
    tex->setImage(osgDB::readImageFile("Images/white.rgb"));

    StateSet *dstate = new StateSet;

    dstate->setTextureAttributeAndModes(0, tex, StateAttribute::OFF );
    dstate->setTextureAttribute(0, new TexEnv );
    dstate->setMode( GL_LIGHTING, StateAttribute::OFF );
    dstate->setMode( GL_CULL_FACE, StateAttribute::ON );
    

    // clear the depth to the far plane.
    osg::Depth* depth = new osg::Depth;
    depth->setFunction(osg::Depth::ALWAYS);
    depth->setRange(1.0,1.0);   
    dstate->setAttributeAndModes(depth,StateAttribute::ON );

    dstate->setRenderBinDetails(-2,"RenderBin");

    geom->setStateSet( dstate );

    Geode *geode = new Geode;
    geode->addDrawable( geom );

    geode->setName( "Sky" );

    return geode;
}

OSGVisualizer::OSGVisualizer()
{
    root = new osg::Group();

	geode = new Geode();
	lineGeom = new osg::Geometry();
	vertices = new osg::Vec3Array();
	geode->addDrawable( lineGeom );
	root->addChild(geode);

	root->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

	viewer.setUpViewInWindow(100, 100, 800, 600);
    viewer.setSceneData( root );

    viewer.setCameraManipulator(new osgGA::TrackballManipulator());
}

void OSGVisualizer::renderScene()
{

	pthread_mutex_lock(&sceneLock);
	vertices = new osg::Vec3Array();

	if( scene != NULL) 
	{
		for (int i = 0; i < (int) scene->lines.size(); i++)
		{
			
			for(int j = 0; j < scene->lines[i].getLines().size(); j+=3)
			{
				vertices->push_back( osg::Vec3(scene->lines[i].getLines().at(j),
												scene->lines[i].getLines().at(j+1),
												scene->lines[i].getLines().at(j+2)
												)
												);
			}

			lineGeom->setVertexArray(vertices);
		/*	
			osg::Vec4Array * colours = new osg::Vec4Array;
			colours->push_back(osg::Vec4(scene->lines[i].getColor()[0], 
						scene->lines[i].getColor()[1],
						scene->lines[i].getColor()[2],
						0.0f
						)
						);
	        lineGeom->setColorArray(colours);

			lineGeom->setColorBinding(osg::Geometry::BIND_PER_PRIMITIVE);

			// set the normal in the same way color.
			osg::Vec3Array* normals = new osg::Vec3Array;
			normals->push_back(osg::Vec3(0.0f,-1.0f,0.0f));
			lineGeom->setNormalArray(normals);
			lineGeom->setNormalBinding(osg::Geometry::BIND_PER_PRIMITIVE);
*/
			lineGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,
											scene->lines[i].getLines().size()));

		}
	}

	pthread_cond_signal(&sceneHasBeenDrawn);
	pthread_mutex_unlock(&sceneLock);
}

void OSGVisualizer::go()
{
	this->viewer.realize();
//	this->root->addChild(makeSky());
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

	app.go();

}
