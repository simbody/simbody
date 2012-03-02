# Makefile for building Simbody examples on Linux or Mac.
# See README_Makefile.txt for more information.
 
# Comment this out if you're building 64 bit examples (Linux only).
M32FLAG = -m32

# Uncomment DEBUG to use Simbody debug libraries instead of release libraries
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

# You would need this whole list to link with static libraries; for 
# dynamic libraries (shared objects) the others are picked up as needed.
#       -llapack -lblas -lpthread -lm -lrt
LIBS= -lSimTKsimbody$(DEBUG) -lSimTKmath$(DEBUG) -lSimTKcommon$(DEBUG)

ifeq "${DEBUG}" ""
CFLAGS = -O2 ${M32FLAG}
else
CFLAGS = -g ${M32FLAG}
endif

ALL_PROGS = \
SimbodyInstallTestNoViz \
SimbodyInstallTest \
ChainExample \
ExampleAmysIKProblem \
ExampleAssemblerPlayground \
ExampleBricardMechanism \
ExampleContactPlayground \
ExampleEventHandler  \
ExampleEventReporter  \
ExampleGears  \
ExampleKneeJoint \
ExampleLongPendulum \
ExamplePendulum \
ExampleSampleAndHold \
JaredsDude \
Rattleback

default: SimbodyInstallTestNoViz SimbodyInstallTest
	@echo Using $(SimTK_HOME) as the SimTK installation directory.


all : $(ALL_PROGS)
	@echo Using $(SimTK_HOME) as the SimTK installation directory.

# Treat all .cpp source files the same way
.cpp : 
	g++ $(CFLAGS) $< -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -o $*

clean : 
	rm $(ALL_PROGS)
