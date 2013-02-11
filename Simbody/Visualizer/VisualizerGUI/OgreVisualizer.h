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
	virtual void drawBox() {}
	virtual void drawCylinder() {}
	virtual void drawSphere() {}
	virtual void drawCircle() {}

	virtual void renderScene();
	
private:
};

#endif // #ifndef __OgreVisualizer_h_
