#pragma once

#include "../types.h"
#include <memory>

namespace sdfibm::forcer {

#define TYPENAME(name) \
    static std::string typeName() {return name;} \
    static bool added;

class IForcer;
template <typename T>
class _creator
{
public:
    static std::unique_ptr<IForcer> create(const dictionary& para)
    {
        return std::make_unique<T>(para);
    }
};

class IForcer
{
public:
    using Force = std::pair<vector, vector>;
    IForcer() = default;
    virtual ~IForcer() = default;

    virtual Force generate(
        const scalar& time,
        const vector& position,
        const vector& velocity,
        const quaternion& orientation,
        const vector& omega
        ) = 0;
    virtual std::string description() const = 0;
};

}
