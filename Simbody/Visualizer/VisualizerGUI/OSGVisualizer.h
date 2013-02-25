#ifndef _OSGVISUALIZER_H_
#define _OSGVISUALIZER_H_

#include <osg/Node>
#include <osg/Group>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Texture2D>
#include <osg/TexEnv>
#include <osg/Depth>
#include <osg/StateSet>
#include <osgDB/ReadFile> 
#include <osgViewer/Viewer>
#include <osg/PositionAttitudeTransform>
#include <osgGA/TrackballManipulator>

#include <osg/MatrixTransform>

#include "VisualizerBase.h"

using namespace osg;

class OSGVisualizer : public VisualizerBase
{
public:
	OSGVisualizer();
	void go();

	void createScene();

	osgViewer::Viewer viewer;
	osg::Group* root;
	osg::Geode* sceneGeode;

protected:
	//	Draw functions
	virtual void drawLine(RenderedLine& line);
	virtual void drawBox() {}
	virtual void drawCylinder() {}
	virtual void drawSphere() {}
	virtual void drawCircle() {}

	void renderScene();

};

#endif // _OSGVISUALIZER_H_
