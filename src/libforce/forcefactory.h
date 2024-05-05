#ifndef FORCEFACTORY_H
#define FORCEFACTORY_H

#include <map>
#include <memory>
#include "../types.h"
namespace sdfibm::force {

class IForce;

class ForceFactory
{
public:
    using TCreateMethod = IForce* (*)(const dictionary&);

public:
    ForceFactory() = delete;

    static bool add(const std::string& name, TCreateMethod create_method);
    static IForce* create(const std::string& name, const dictionary&);

    static void report(std::ostream& os = std::cout)
    {
        // displays all the forces the factory can "produce"
        int i = 1;
        for(auto it : m_methods)
        {
            os << '[' << i << "] " << it.first << std::endl;
            ++i;
        }
    }

private:
    static std::map<std::string, TCreateMethod> m_methods;
};

}
#endif
