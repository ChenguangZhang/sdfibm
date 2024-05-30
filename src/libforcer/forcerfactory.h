#pragma once

#include "../types.h"
#include "../genericfactory.h"
#include "iforcer.h"
namespace sdfibm::forcer {
    MAKESPECIALFACTORY(Forcer, IForcer, Foam::dictionary);
    #define REGISTERFORCE(m) \
        bool sdfibm::m::added = ForcerFactory::add(sdfibm::m::typeName(), sdfibm::m::create);

    #include "constant.h"
    REGISTERFORCE(forcer::Constant);
    #include "spring.h"
    REGISTERFORCE(forcer::Spring);
}