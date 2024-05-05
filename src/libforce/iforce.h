#ifndef IMOTION_H_INCLUDED
#define IMOTION_H_INCLUDED

#include <utility>
#include "../types.h"

namespace sdfibm::force {

#define TYPENAME(name) \
    static std::string typeName() {return name;} \
    static bool added;

class IForce;
template <typename T>
class _creator
{
public:
    static IForce* create(const dictionary& node)
    {
        return new T(node);
    }
};

class IForce
{
public:
    using Force = std::pair<vector, vector>;
    IForce() = default;
    virtual ~IForce() = default;

    virtual Force generate(const scalar& time, const vector& position, const vector& velocity) = 0;
    virtual std::string description() const = 0; // more detailed information
};

}
#endif
