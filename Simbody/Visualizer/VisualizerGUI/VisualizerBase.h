#ifndef _VISUALIZER_BASE_H_
#define _VISUALIZER_BASE_H_

#include <string.h>
#include <sstream>
#include "SimTKcommon.h"
#include "../src/VisualizerProtocol.h"
#include <errno.h>
#include <stdio.h>
#include <pthread.h>
#include "Scene.h"

// Next, get the functions necessary for reading from and writing to pipes.
#ifdef _WIN32
    #include <io.h>
    #define READ _read
#else
    #include <unistd.h>
    #define READ read
#endif

// gcc 4.4.3 complains bitterly if you don't check the return
// status from the write() system call. This avoids those
// warnings and maybe, someday, will catch an error.
#define WRITE(pipeno, buf, len) \
	{ int status=write((pipeno), (buf), (len)); \
    SimTK_ERRCHK4_ALWAYS(status!=0, "VisualizerGUI",  \
    "An attempt to write() %d bytes to pipe %d failed with errno=%d (%s).", \
    (len),(pipeno),errno,strerror(errno));}

class VisualizerBase 
{

public:
//	class for communication
	struct VisualizerCom
	{
		VisualizerBase* arg;
		int inPipe;
		int outPipe;
	};

	Scene * scene;

public:
	VisualizerBase();
	~VisualizerBase();
	
	void shakeHandsWithSimulator(int inPipe, int outPipe);
	void readDataFromPipe(int srcPipe, unsigned char* buffer, int bytes);

	// Hack job atm
	static void* listenForInputHelper(void* args) 
	{
		VisualizerCom *vcom = (VisualizerCom*) args;
		vcom->arg->listenForInput(vcom->inPipe, vcom->outPipe);
	}

	void dumpMessage();
protected:
	std::string getName();
	void setName(std::string appName);
   	std::string	getVersion();
	std::string getExecutableName();

//	Drawing functions
	virtual void drawBox() = 0;
	virtual void drawCylinder() = 0;
	virtual void drawSphere() = 0;
	virtual void drawCircle() = 0;

	void listenForInput(int inPipe, int outPipe);
	Scene* readNewScene(int inPipe);
	virtual void renderScene() = 0;

//	Dummy Display pointers (Taken from gazebo)
	void * dummyDisplay;
	unsigned long long dummyWindowId;
	void * dummyContext;

private:
	std::string name;
	std::string simbodyVersionStr;
	std::string simulatorExecutableName;

}; // VisualizerBase


#endif //_VISUALIZER_BASE_H
