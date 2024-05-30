set(libName force)

set(
    libSrc
    forcerfactory.cpp
)

set(
    libInc
    iforcer.h
    constant.h
    spring.h
    forcerfactory.h
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
