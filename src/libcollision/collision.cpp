#include "collision.h"
namespace sdfibm{
/* 0. shapeShapeCollision takes two solids, and two references of the collision position (cP) and normal (cN)
   1. it returns the overlap between two shapes -- collision occur when overlap > 0
   2. the normal points from the first towards the second shape */

real sphereSphereCollision(const Solid& s1, const Solid& s2, vector& cP, vector& cN)
{
    Foam::vector s2s = s2.getCenter() - s1.getCenter();
    cN = Foam::normalised(s2s);

    real r1 = s1.getShape()->getRadiusB();
    real r2 = s2.getShape()->getRadiusB();

    cP = 0.5*(s1.getCenter() + s2.getCenter() + (r1 - r2)*cN);
    return r1 + r2 - Foam::mag(s2s);
}
real circleCircleCollision(const Solid& s1, const Solid& s2, vector& cP, vector& cN) {
    return sphereSphereCollision(s1, s2, cP, cN); }

real planeSphereCollision(const Solid& p, const Solid& s, vector& cP, vector& cN)
{
    vector s_center = Foam::conjugate(p.getOrientation()).transform(s.getCenter() - p.getCenter());
    real p2s = s_center.y();
    cN = p.getOrientation().transform(Foam::vector(0,1,0));
    cP = s.getCenter() - s.getRadiusB()*cN;
    return s.getRadiusB() - p2s;
}

real planeCircleCollision(const Solid& p, const Solid& s, vector& cP, vector& cN) {
    return planeSphereCollision(p, s, cP, cN); }
real spherePlaneCollision(const Solid& s, const Solid& p, vector& cP, vector& cN) {
    return planeSphereCollision(p, s, cP, cN); }
real circlePlaneCollision(const Solid& s, const Solid& p, vector& cP, vector& cN) {
    return spherePlaneCollision(s, p, cP, cN); }

collisionFunc collisionFuncTable[kNUM_COL_FUNC][kNUM_COL_FUNC];

void InitCollisionFuncTable()
{
    for(label i = 0; i < kNUM_COL_FUNC; i++){
        for(label j = 0; j < kNUM_COL_FUNC; j++) {
            collisionFuncTable[i][j] = nullptr;
    } }

    collisionFuncTable[SHAPE2ID["Circle"]][SHAPE2ID["Plane"] ] = circlePlaneCollision;
    collisionFuncTable[SHAPE2ID["Plane"] ][SHAPE2ID["Circle"]] = planeCircleCollision;
    collisionFuncTable[SHAPE2ID["Sphere"]][SHAPE2ID["Plane"] ] = spherePlaneCollision;
    collisionFuncTable[SHAPE2ID["Plane"] ][SHAPE2ID["Sphere"]] = planeSphereCollision;
    collisionFuncTable[SHAPE2ID["Circle"]][SHAPE2ID["Circle"]] = circleCircleCollision;
    collisionFuncTable[SHAPE2ID["Sphere"]][SHAPE2ID["Sphere"]] = sphereSphereCollision;
}

collisionFunc getCollisionFunc(const std::string& name1, const std::string& name2)
{
    return collisionFuncTable[SHAPE2ID[name1]][SHAPE2ID[name2]];
}


}
