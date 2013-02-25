/*
 * OgreVisualizer.h
 * Taken from Ogre wiki
*/
#ifndef __OgreVisualizer_h_
#define __OgreVisualizer_h_

#include "OgreVisualizerBase.h"
#include "VisualizerBase.h"

class OgreVisualizer : public OgreVisualizerBase, public VisualizerBase
{
public:
    OgreVisualizer();
    virtual ~OgreVisualizer();

	void go();

protected:
    virtual void createScene();

//	Draw functions
	virtual void drawLine(RenderedLine& line, const std::string& name);
	virtual void drawText(RenderedText& text, const std::string& name);
	virtual void drawMesh(RenderedMesh& mesh, const std::string& name);
	virtual void drawBox() {}
	virtual void drawCylinder() {}
	virtual void drawSphere() {}
	virtual void drawCircle() {}

	virtual void renderScene();
	
	Ogre::SceneNode * backgroundNode;
	Ogre::SceneNode * mainNode;

	int lineCount;

private:
};

#endif // #ifndef __OgreVisualizer_h_
