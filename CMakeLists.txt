cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(sdfibm LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(DEFINED ENV{WM_PROJECT_DIR})
	MESSAGE(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
    message(FATAL_ERROR "Cannot find OpenFOAM")
endif(DEFINED ENV{WM_PROJECT_DIR})

set(FOAMVERSION $ENV{WM_PROJECT_VERSION}) 
set(FOAMDIR $ENV{WM_PROJECT_DIR})
set(FOAMLIB $ENV{FOAM_LIBBIN})
set(FOAMSRC $ENV{FOAM_SRC})
set(FOAMMPI $ENV{FOAM_MPI})

set(FOAMFLAG "-std=c++17 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wold-style-cast
-Wnon-virtual-dtor -Wno-unused-variable -Wno-unused-local-typedef -Wno-invalid-offsetof -Wno-deprecated-register -Wno-undefined-var-template -DNoRepository -ftemplate-depth-100 -fPIC")

set(FOAMINC
    ${FOAMSRC}/finiteVolume/lnInclude 
    ${FOAMSRC}/OpenFOAM/lnInclude  
    ${FOAMSRC}/OSspecific/POSIX/lnInclude 
    ${FOAMSRC}/meshTools/lnInclude 
    ${FOAMSRC}/dynamicFvMesh/lnInclude 
)

link_directories(
    ${FOAMLIB}
    ${FOAMLIB}/${FOAMMPI}
)

add_subdirectory(src)
