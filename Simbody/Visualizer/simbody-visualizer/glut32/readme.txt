On Windows we supply our own glut includes and libraries. We expect Windows to
provide only GL/gl.h and GL/glu.h; we'll provide glext.h and glut.h, as well as glut32.lib & dll binaries in 32 and 64 bits.

Note "glut32" is generically the name of glut libraries on Windows; it does
not imply 32 bit binaries. Sadly, even the 64 bit binary is named glut32.

The version of glut we're providing was built from glut 3.7.6 source, patched to
support the mouse wheel using a patch obtained here: 
    http://www.realmtech.net/opengl/glut.php.
We also fixed a bug that prevented glutTimerFunc() from working correctly.

During a Windows build (32 or 64 bit), we will copy the 32 or 64 bit dll into
the local binary directory in which the VisualizationGUI.exe executable is
begin built. That ensures that it will use our glut32.dll regardless of what's
on the path. When the INSTALL target is built, we'll copy that dll into the
SimTK_INSTALL_DIR/bin directory along with VisualizationGUI.exe. 

For Mac and Linux we expect to find an installed glut or freeglut that we can
use to build the Simbody Visualizer.

Sherm 11/16/2010

