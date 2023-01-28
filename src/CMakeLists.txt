add_subdirectory(libshape)
add_subdirectory(libmotion)
add_subdirectory(libmaterial)
add_subdirectory(libcollision)

add_executable(${PROJECT_NAME} main.cpp cellenumerator.cpp geometrictools.cpp logger.cpp meshinfo.cpp solid.cpp solidcloud.cpp)
add_definitions(${FOAMFLAG})

target_include_directories(${PROJECT_NAME} PUBLIC $FOAMINC)
target_link_libraries(${PROJECT_NAME} PUBLIC
    OpenFOAM finiteVolume meshTools Pstream
    dl m motion shape collision)

