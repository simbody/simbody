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

protected:
    virtual void createScene();

	virtual void preRender();
//	Draw functions
	virtual void drawBox() {}
	virtual void drawCylinder() {}
	virtual void drawSphere() {}
	virtual void drawCircle() {}
	virtual void drawGroundAndSky();

	virtual void startRendering();
	virtual void finishRendering();
	
private:
};

#endif // #ifndef __OgreVisualizer_h_
