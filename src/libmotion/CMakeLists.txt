set(libName motion)

set(
    libSrc
    motionfactory.cpp
)

set(
    libInc
    imotion.h
    motion01mask.h
    motion000002.h
    motion110002.h
    motion222000.h
    motionfactory.h
    motionopenclose.h
    motionrotor.h
    motionsinedirectional.h
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
