#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "types.hpp"
#include "utilsccd.hpp"
#include "utils.hpp"
#include "solid.hpp"
#include "../shapelib/shapes.hpp"

real circleCircleCollision(Solid& c1, Solid& c2, vector& collision_location, vector& collision_direction);
real circlePlaneCollision (Solid& c,  Solid& p,  vector& collision_location, vector& collision_direction);
real spherePlaneCollision (Solid& c,  Solid& p,  vector& collision_location, vector& collision_direction);
real ellipsePlaneCollision(Solid& c,  Solid& p,  vector& collision_location, vector& collision_direction);
real sphereSphereCollision(Solid& c1,  Solid& c2,  vector& collision_location, vector& collision_direction);

extern real (*collisionFunctionTable[10][10])(
           Solid& s1,
           Solid& s2,
           vector& contact_location,
           vector& contact_normal
        );

void collisionFunctionTableInit();
#endif // COLLISION_HPP
