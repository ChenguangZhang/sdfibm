#ifndef IMOTION_H_INCLUDED
#define IMOTION_H_INCLUDED

#include "../types.h"

namespace sdfibm{
// nameing convention: 0 completely frozen, 1 totally free, 2 user specified

#define TYPENAME(name) \
    static std::string typeName() {return name;} \
    static bool added;

// define templated base to insert creator function into child class
// equivalent to manually adding the following line into each child class
//static IMotion* create(const Node& node){return new Motion000000(node);}
class IMotion;
template <typename T>
class _creator
{
public:
    static IMotion* create(const dictionary& node)
    {
        return new T(node);
    }
};

class IMotion
{
public:
    IMotion() = default;
    virtual ~IMotion() = default;

    // the work horse
    virtual void constraint(const scalar& time, vector& velocity, vector& omega) = 0;
    virtual std::string description() const = 0; // more detailed information
};

}
#endif
