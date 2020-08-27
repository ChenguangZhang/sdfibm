#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "../types.h"
#include "../utils.h"
#include "../solid.h"
#include <map>
namespace sdfibm{

const int kNUM_COL_FUNC = 10;
static std::map<std::string, int> SHAPE2ID = {
    {"Plane",  0},
    {"Circle", 1},
    {"Sphere", 2}
};

// 2d collision
scalar circlePlaneCollision (const Solid&  c, const Solid&  p, vector& cP, vector& cN);
scalar circleCircleCollision(const Solid& c1, const Solid& c2, vector& cP, vector& cN);

// 3d collision
scalar spherePlaneCollision (const Solid&  c, const Solid&  p, vector& cP, vector& cN);
scalar sphereSphereCollision(const Solid& c1, const Solid& c2, vector& cP, vector& cN);

typedef scalar (*collisionFunc)(const Solid& s1, const Solid& s2, vector& cP, vector& cN);
extern collisionFunc collisionFuncTable[kNUM_COL_FUNC][kNUM_COL_FUNC]; // table of collsion handling functions
// extern scalar (*collisionFuncTable[kNUM_COL_FUNC][kNUM_COL_FUNC])(Solid& s1, Solid& s2, vector& cP, vector& cN);
void InitCollisionFuncTable();
collisionFunc getCollisionFunc(const std::string& name1, const std::string& name2);

}
#endif // COLLISION_HPP
