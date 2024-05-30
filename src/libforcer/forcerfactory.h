#pragma once

#include "../types.h"
#include "../genericfactory.h"
#include "iforcer.h"
namespace sdfibm::force {
    MAKESPECIALFACTORY(Forcer, IForcer, Foam::dictionary);
    #define REGISTERFORCE(m) \
        bool sdfibm::m::added = ForcerFactory::add(sdfibm::m::typeName(), sdfibm::m::create);

    #include "constant.h"
    REGISTERFORCE(force::Constant);
    #include "spring.h"
    REGISTERFORCE(force::Spring);
}