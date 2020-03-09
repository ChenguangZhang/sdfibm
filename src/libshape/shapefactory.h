#ifndef SHAPEFACTORY_H
#define SHAPEFACTORY_H

#include <map>
#include <memory>
#include "../types.h"
namespace sdfibm{

class IShape;

class ShapeFactory
{
public:
    using TCreateMethod = IShape* (*)(const dictionary&);

public:
    // could be a singleton, for now just disable ctor
    ShapeFactory() = delete;

    static bool add(const string& name, TCreateMethod create_method);
    static IShape* create(const string& name, const dictionary&);

    static void report(std::ostream& os = std::cout)
    {
        // displays all the motions the factory can "produce"
        int i = 1;
        for(auto it : m_methods)
        {
            os << '[' << i << "] " << it.first << std::endl;
            ++i;
        }
    }

private:
    static std::map<string, TCreateMethod> m_methods;
};

}
#endif
