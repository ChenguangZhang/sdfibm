#include "iforce.h"
#include "forcefactory.h"

#include <iostream>
namespace sdfibm::force {

#define REGISTERFORCE(m) \
    bool sdfibm::m::added = ForceFactory::add(sdfibm::m::typeName(), sdfibm::m::create);

bool ForceFactory::add(const std::string& name, TCreateMethod create_method)
{
    auto it = m_methods.find(name);
    // std::cout << "Registering..." << name << std::endl;
    if(it == m_methods.end())
    {
        m_methods[name] = create_method;
        return true;
    }
    return false;
}

IForce* ForceFactory::create(const std::string& name, const dictionary& node)
{
    auto it = m_methods.find(name);
    if(it != m_methods.end())
    {
        return it->second(node);
    }
    return nullptr;
}

std::map<std::string, ForceFactory::TCreateMethod> ForceFactory::m_methods;

}
// add forces
#include "constant.h"
REGISTERFORCE(force::Constant);
#include "spring.h"
REGISTERFORCE(force::Spring);
