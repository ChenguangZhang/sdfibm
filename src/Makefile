CC=g++
CFLAGS =-m64 -Dlinux32 -DWM_DP

C_DEBUG_FLAGS = -ggdb -O0
C_RELEASE_FLAGS = -s -O3

CFLAGS+=-Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor
CFLAGS+=-DNoRepository -ftemplate-depth-100 -fPIC -std=c++11
CFLAGS+=-DWM_ARCH_OPTION=64 -DWM_LABEL_SIZE=32

INFLAGS =-I.
INFLAGS+=-IlnInclude
INFLAGS+=-I$(FOAM_SRC)/finiteVolume/lnInclude
INFLAGS+=-I$(FOAM_SRC)/OpenFOAM/lnInclude
INFLAGS+=-I$(FOAM_SRC)/dynamicFvMesh/lnInclude
INFLAGS+=-I$(FOAM_SRC)/OSspecific/POSIX/lnInclude
INFLAGS+=-I$(FOAM_SRC)/meshTools/lnInclude

LIBMPI=-L$(FOAM_LIBBIN)/$(FOAM_MPI) -lPstream

LDFLAGS =-lm
LDFLAGS+=-L$(FOAM_LIBBIN) -lfiniteVolume -lOpenFOAM -lmeshTools


APPSOURCES=main.cpp solid.cpp solidcloud.cpp logger.cpp meshinfo.cpp cellenumerator.cpp geometrictools.cpp
APPOBJECTS=$(APPSOURCES:.cpp=.o)

APPTARGET=sdfibm
LIBTARGETLINK=-L. -L./libmotion -lmotion -L./libshape -lshape -L./libshape -lshape -L./libcollision -lcollision

# release
.PHONY	: RELEASE
RELEASE	: CFLAGS+=$(C_RELEASE_FLAGS)
RELEASE	: all
# debug
.PHONY	: DEBUG
DEBUG 	: CFLAGS+=$(C_DEBUG_FLAGS)
DEBUG 	: all

SUBDIRS= ./libcollision ./libmotion ./libshape
.PHONY: $(SUBDIRS) $(APPTARGET)

all: $(SUBDIRS) $(APPTARGET)

$(SUBDIRS):
	$(MAKE) -C $@

$(APPTARGET): $(APPOBJECTS)
	@echo Building application $(APPTARGET)
	$(CC) -o $@ $(APPOBJECTS) $(LIBMPI) $(LIBTARGETLINK) $(LDFLAGS)

clean:
	rm $(LIBOBJECTS) $(APPOBJECTS) $(APPTARGET)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(INFLAGS) -c
