# Uncomment DEBUG to use SimTK debug libraries instead of release libraries
# DEBUG=_d

# Default install directory
# SimTK_INSTALL_DIR=/usr/local/SimTK
#
ifeq "${SimTK_INSTALL_DIR}" ""
SimTK_HOME=/usr/local/SimTK
else
SimTK_HOME=${SimTK_INSTALL_DIR}
endif

LIB_DIR=$(SimTK_HOME)/lib
INCLUDE_DIR=$(SimTK_HOME)/include

LIBS= -lSimTKsimbody$(DEBUG)  -lpthread -lm -lSimTKmolmodel$(DEBUG) -lSimTKsimbody_aux$(DEBUG) -lSimTKsimbody$(DEBUG) -lSimTKmath$(DEBUG) -lSimTKcpodes$(DEBUG) -lSimTKcommon$(DEBUG) -lSimTKlapack \
      -lvtkCommon -lvtkHybrid -lvtkGenericFiltering -lvtkVolumeRendering -lvtksys -lvtkFiltering -lvtkRendering \
      -lvtkIO -lvtkGraphics -lvtkexoIIc -lvtkNetCDF -lvtkImaging -lvtkftgl -lvtkfreetype -lvtkDICOMParser -lvtkpng  \
      -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat 

# Note that we're compiling the examples with optimization here "-O2";
# if you want debug, use the debug flag "-g" instead of "-O2".
# They will run much slower in debug mode.
CFLAGS = -O2 -m32
#CFLAGS = -g -m32

ALL_PROGS = \
ExampleAdenylateMobilitiesVMD \
ExampleAdenylateMobilitiesVTK \
ExampleAmysIKProblem \
ExampleAssemblerPlayground \
ExampleBricardMechanism \
ExampleChain \
ExampleContactPlayground \
ExampleCreatePropane \
ExampleEventHandler  \
ExampleEventReporter  \
ExampleGears  \
ExampleKneeJoint \
ExampleLoadPdb  \
ExampleLoadPdb2  \
ExamplePendulum \
ExampleRigidProtein  \
ExampleSimpleProtein  \
ExampleSimpleRNA \
ExampleTwoArgonAtoms  \
ExampleTwoEthanes 

default: ExamplePendulum
	@echo Expecting $(SimTK_HOME) to be the SimTK 2.0 installation directory.


all : $(ALL_PROGS)
	@echo Expecting $(SimTK_HOME) to be the SimTK 2.0 installation directory.

# Treat all .cpp source files the same way
.cpp : 
	g++ $(CFLAGS) $< -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -o $*

clean : 
	rm $(ALL_PROGS)

