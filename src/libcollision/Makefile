EXTDIR=$(CURDIR)/../../external

CC=g++
CFLAGS =-m64 -Dlinux32 -DWM_DP -O2 # remove O0 when no bug
CFLAGS+=-Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor
CFLAGS+=-DNoRepository -ftemplate-depth-100 -fPIC -std=c++14
CFLAGS+=-DWM_ARCH_OPTION=64 -DWM_LABEL_SIZE=32

INFLAGS =-I.
INFLAGS+=-I$(EXTDIR)/
INFLAGS+=-IlnInclude
INFLAGS+=-I$(FOAM_SRC)/finiteVolume/lnInclude
INFLAGS+=-I$(FOAM_SRC)/OpenFOAM/lnInclude
INFLAGS+=-I$(FOAM_SRC)/dynamicFvMesh/lnInclude
INFLAGS+=-I$(FOAM_SRC)/OSspecific/POSIX/lnInclude
INFLAGS+=-I$(FOAM_SRC)/meshTools/lnInclude

LIBMPI=$(FOAM_LIBBIN)/openmpi-system/libPstream.so

LDFLAGS =-lm
LDFLAGS+=-L$(FOAM_LIBBIN) -lfiniteVolume -lOpenFOAM -ldynamicFvMesh -L$(EXTDIR) -lyaml-cpp

LIBSOURCES= ugrid.cpp collision.cpp
LIBOBJECTS=$(LIBSOURCES:.cpp=.o)


LIBTARGET=./libcollision.a

all: $(LIBTARGET)
lib: $(LIBTARGET)


$(LIBTARGET): $(LIBOBJECTS) 
	# $(CC) $(LIBOBJECTS) -c -o $@ $(LDFLAGS)
	ar rcs $@ $^

clean:
	rm $(LIBOBJECTS) $(LIBTARGET)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(INFLAGS) -c
