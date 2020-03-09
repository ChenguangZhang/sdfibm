#include "ishape.h"
#include "shapefactory.h"

#include <iostream>
namespace sdfibm{

#define REGISTERSHAPE(m) \
    bool sdfibm::m::added = ShapeFactory::add(sdfibm::m::typeName(), sdfibm::m::create);

bool ShapeFactory::add(const string& name, TCreateMethod create_method)
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

IShape* ShapeFactory::create(const string& name, const dictionary& para)
{
    auto it = m_methods.find(name);
    if(it != m_methods.end())
    {
        return it->second(para);
    }
    return nullptr;
}

std::map<string, ShapeFactory::TCreateMethod> ShapeFactory::m_methods;

}

// add shapes
#include "circle.h"
REGISTERSHAPE(Circle);

#include "sphere.h"
REGISTERSHAPE(Sphere);

#include "ellipse.h"
REGISTERSHAPE(Ellipse);

#include "ellipsoid.h"
REGISTERSHAPE(Ellipsoid);

#include "rectangle.h"
REGISTERSHAPE(Rectangle);

#include "circle_tail.h"
REGISTERSHAPE(Circle_Tail);

#include "circle_twotail.h"
REGISTERSHAPE(Circle_TwoTail);

#include "box.h"
REGISTERSHAPE(Box);

#include "plane.h"
REGISTERSHAPE(Plane);

