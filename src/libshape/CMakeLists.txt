set(libName shape)

set(
    libSrc
    shapefactory.cpp
)

set(
    libInc
    ishape.h
    shapefactory.h
    box.h
    circle.h
    circle_tail.h
    circle_twotail.h
    ellipse.h
    ellipsoid.h
    plane.h
    rectangle.h
    sphere.h
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
