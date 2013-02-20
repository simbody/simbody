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

#include "VisualizerBase.h"

using namespace osg;

class OSGVisualizer : public VisualizerBase
{
public:
	OSGVisualizer();
	void go();
	osg::Node * makeSky();

	osgViewer::Viewer viewer;
	osg::Group* root;
	osg::Geode* geode;
	osg::Geometry* lineGeom;
	osg::Vec3Array* vertices;

protected:
	//	Draw functions
	virtual void drawBox() {}
	virtual void drawCylinder() {}
	virtual void drawSphere() {}
	virtual void drawCircle() {}

	void renderScene();

};

#endif // _OSGVISUALIZER_H_
