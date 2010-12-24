# Uncomment DEBUG to use SimTK debug libraries instead of release libraries
#DEBUG=_d

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

LIBS= -lSimTKsimbody$(DEBUG) -lSimTKmath$(DEBUG) -lSimTKcommon$(DEBUG) -lSimTKlapack -lpthread -lm -lrt

ifeq "${DEBUG}" ""
CFLAGS = -O2
else
CFLAGS = -g
endif

ALL_PROGS = \
ExampleAmysIKProblem \
ExampleAssemblerPlayground \
ExampleBricardMechanism \
ExampleChain \
ExampleContactPlayground \
ExampleEventHandler  \
ExampleEventReporter  \
ExampleGears  \
ExampleKneeJoint \
ExamplePendulum 

default: ExamplePendulum
	@echo Using $(SimTK_HOME) as the SimTK installation directory.


all : $(ALL_PROGS)
	@echo Using $(SimTK_HOME) as the SimTK installation directory.

# Treat all .cpp source files the same way
.cpp : 
	g++ $(CFLAGS) $< -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -o $*

clean : 
	rm $(ALL_PROGS)
