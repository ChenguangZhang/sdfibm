#pragma once

#include "../types.h"
#include "../genericfactory.h"
#include "ishape.h"
#include <functional>
#include <cassert>
namespace sdfibm {
    MAKESPECIALFACTORY(Shape, IShape, Foam::dictionary);

    #define REGISTERSHAPE(m) \
        bool sdfibm::m::added = ShapeFactory::add(sdfibm::m::typeName(), sdfibm::m::create);

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
}
