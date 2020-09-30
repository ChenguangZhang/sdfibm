CC=g++
CFLAGS =-m64 -Dlinux32 -DWM_DP

C_DEBUG_FLAGS = -ggdb -O0
C_RELEASE_FLAGS = -s -O3

CFLAGS+=-Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor
CFLAGS+=-DNoRepository -ftemplate-depth-100 -fPIC -std=c++11
CFLAGS+=-DWM_ARCH_OPTION=64 -DWM_LABEL_SIZE=32

INFLAGS =-I.
INFLAGS+=-I../src
INFLAGS+=-IlnInclude
INFLAGS+=-I$(FOAM_SRC)/finiteVolume/lnInclude
INFLAGS+=-I$(FOAM_SRC)/OpenFOAM/lnInclude
INFLAGS+=-I$(FOAM_SRC)/dynamicFvMesh/lnInclude
INFLAGS+=-I$(FOAM_SRC)/OSspecific/POSIX/lnInclude
INFLAGS+=-I$(FOAM_SRC)/meshTools/lnInclude

LIBMPI=-L$(FOAM_LIBBIN)/openmpi-system -lPstream

LDFLAGS =-lm
LDFLAGS+=-L$(FOAM_LIBBIN) -lfiniteVolume -lOpenFOAM -ldynamicFvMesh -lmeshTools


APPSOURCES=main.cpp solidcloud.cpp
APPOBJECTS=$(APPSOURCES:.cpp=.o) ../src/meshinfo.o ../src/cellenumerator.o ../src/geometrictools.o

APPTARGET=smoothvof
LIBTARGETLINK=-L. -L../src/libmotion -lmotion -L../src/libshape -lshape -L../src/libshape -lshape -L../src/libcollision -lcollision

# release
.PHONY	: RELEASE
RELEASE	: CFLAGS+=$(C_RELEASE_FLAGS)
RELEASE	: all
# debug
.PHONY	: DEBUG
DEBUG 	: CFLAGS+=$(C_DEBUG_FLAGS)
DEBUG 	: all

.PHONY: $(APPTARGET)

all: $(APPTARGET)

$(APPTARGET): $(APPOBJECTS)
	@echo Building application $(APPTARGET)
	$(CC) -o $@ $(APPOBJECTS) $(LIBMPI) $(LIBTARGETLINK) $(LDFLAGS)

clean:
	rm $(LIBOBJECTS) $(APPSOURCES:.cpp=.o) $(APPTARGET)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(INFLAGS) -c
