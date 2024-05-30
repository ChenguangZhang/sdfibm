add_subdirectory(libshape)
add_subdirectory(libmotion)
add_subdirectory(libforcer)
add_subdirectory(libmaterial)
add_subdirectory(libcollision)
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")
    add_subdirectory(test)
endif()

add_library(sdfibm_util STATIC cellenumerator.cpp geometrictools.cpp meshinfo.cpp)
target_include_directories(sdfibm_util PUBLIC ${FOAMINC})
target_link_libraries(sdfibm_util PUBLIC m)

add_executable(${PROJECT_NAME} main.cpp logger.cpp solid.cpp solidcloud.cpp)
add_definitions(${FOAMFLAG})

target_include_directories(${PROJECT_NAME} PUBLIC $FOAMINC)
target_link_libraries(${PROJECT_NAME} PUBLIC
    OpenFOAM finiteVolume meshTools Pstream
    dl m motion shape collision sdfibm_util)

