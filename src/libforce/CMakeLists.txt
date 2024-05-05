set(libName force)

set(
    libSrc
    forcefactory.cpp
)

set(
    libInc
    iforce.h
    constant.h
    spring.h
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
