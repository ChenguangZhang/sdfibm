#include "imotion.h"
#include "motionfactory.h"

#include <iostream>

#define REGISTERMOTION(m) \
    bool m::added = MotionFactory::add(m::typeName(), m::create);

bool MotionFactory::add(const string& name, TCreateMethod create_method)
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

IMotion* MotionFactory::create(const string& name, const dictionary& node)
{
    auto it = m_methods.find(name);
    if(it != m_methods.end())
    {
        return it->second(node);
    }
    return nullptr;
}

std::map<string, MotionFactory::TCreateMethod> MotionFactory::m_methods;

// add motions
#include "motion000000.h"
REGISTERMOTION(Motion000000);
#include "motion000100.h"
REGISTERMOTION(Motion000100);
#include "motion000010.h"
REGISTERMOTION(Motion000010);
#include "motion000001.h"
REGISTERMOTION(Motion000001);
#include "motion000111.h"
REGISTERMOTION(Motion000111);
#include "motion000002.h"
REGISTERMOTION(Motion000002);
#include "motion111111.h"
REGISTERMOTION(Motion111111);
#include "motion222000.h"
REGISTERMOTION(Motion222000);
#include "motionsinedirectional.h"
REGISTERMOTION(MotionSineDirectional);

#include      "motion001001.h"
REGISTERMOTION(Motion001001);
#include      "motion010010.h"
REGISTERMOTION(Motion010010);
#include      "motion100100.h"
REGISTERMOTION(Motion100100);
