set(libName collision)

set(libSrc
    ugrid.cpp
    collision.cpp
)

set(libInc
    bbox.h
    ugrid.h
    collision.h
)

add_library(
    ${libName}
    STATIC
    ${libSrc}
    ${libInc}
)

add_definitions(${FOAMFLAG})
target_include_directories(${libName} PUBLIC ${FOAMINC})
target_link_libraries(${libName} PUBLIC m)
