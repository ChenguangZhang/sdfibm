#include "collision.h"

const int MAXCOLFUN = 10;

real (*collisionFunctionTable[MAXCOLFUN][MAXCOLFUN])(
           Solid& s1,
           Solid& s2,
           vector& contact_location,
           vector& contact_normal
        );

real circleCircleCollision(Solid& c1, Solid& c2,
                           vector& collision_location,
                           vector& collision_normal)
{
    collision_location = 0.5*(c1.getCenter() + c2.getCenter());
    collision_normal = Foam::normalised(c1.getCenter() - c2.getCenter());
    real c2c = Foam::mag(c1.getCenter() - c2.getCenter());
    // doc: depth = touch distance - real distance, so when collide there is a
    // postive signal of depth
    return (c1.getShape()->getBoundingRadius()+
            c2.getShape()->getBoundingRadius() ) - c2c;
}

real sphereSphereCollision(Solid& c1, Solid& c2,
                           vector& collision_location,
                           vector& collision_normal)
{
    collision_location = 0.5*(c1.getCenter() + c2.getCenter());
    collision_normal = Foam::normalised(c1.getCenter() - c2.getCenter());
    real c2c = Foam::mag(c1.getCenter() - c2.getCenter());
    return (c1.getShape()->getBoundingRadius()+
            c2.getShape()->getBoundingRadius() ) - c2c;
}

real circlePlaneCollision(Solid& c, Solid& p,
                          vector& collision_location,
                          vector& collision_normal)
{
    // doc: the collision normal all points from the second to the first
    vector newcenter = Foam::conjugate(p.getOrientation()).transform(c.getCenter() - p.getCenter());
    real c2p = newcenter.y();
    vector world_plane_normal = p.getOrientation() * vector(0,1,0);
    collision_location = c.getCenter() - c.getBoundingRadius()*world_plane_normal;
    collision_normal = world_plane_normal;
    return c.getBoundingRadius() - c2p;
}

real spherePlaneCollision(Solid& c, Solid& p,
                          vector& collision_location,
                          vector& collision_normal)
{
    // doc: the collision normal all points from the second to the first
    vector newcenter = Foam::conjugate(p.getOrientation()).transform(c.getCenter() - p.getCenter());
    real c2p = newcenter.y();
    vector world_plane_normal = p.getOrientation() * vector(0,1,0);
    collision_location = c.getCenter() - c.getBoundingRadius()*world_plane_normal;
    collision_normal = world_plane_normal;

    return c.getBoundingRadius() - c2p;
}

real ellipsePlaneCollision(Solid& e, Solid& p,
                          vector& collision_location,
                          vector& collision_normal)
{
    // transform to the ellipse cooridnate
    vector plane_center_t = glm::conjugate(e.getOrientation()) * (p.getCenter() - e.getCenter());
    vector plane_normal_t = glm::conjugate(e.getOrientation()) * p.getOrientation() * Vector(0,1,0);
    plane_normal_t = glm::normalize(plane_normal_t);


    collision_normal = p.getOrientation() * Vector3(0,1,0);

    Ellipse* eshape = dynamic_cast<Ellipse*>(e.getShape());
    real a = eshape->getRadiusa();
    real b = eshape->getRadiusb();
    real nx = plane_normal_t.x;
    real ny = plane_normal_t.y;

    real gammap = 1.0/std::sqrt(a*a*nx*nx + b*b*ny*ny);
    vector pp = gammap * vector(a*a*nx, b*b*ny, 0.0);
    vector pn = -pp;


    real distp = -(pp - plane_center_t) & plane_normal_t;
    real distn = -(pn - plane_center_t) & plane_normal_t;

    if(distn < distp)
    {
        collision_location = e.getOrientation()*pp + e.getCenter();
        return distp;
    }
    else
    {
        collision_location = e.getOrientation()*pn + e.getCenter();
        return distn;
    }
}

void collisionFunctionTableInit()
{
    for(label i = 0; i < 10; i++)
        for(label j = 0; j < 10; j++)
        {
            collisionFunctionTable[i][j] = NULL;
        }
    collisionFunctionTable[SHAPE::CIRC][SHAPE::PLAN] = circlePlaneCollision;
    //collisionFunctionTable[SHAPE::PLAN][SHAPE::CIRC] = circlePlaneCollision;
    //collisionFunctionTable[SHAPE::PLAN][SHAPE::CIRC] = circlePlaneCollision;

    collisionFunctionTable[SHAPE::ELLI][SHAPE::PLAN] = ellipsePlaneCollision;
    //collisionFunctionTable[SHAPE::PLAN][SHAPE::ELLI] = ellipsePlaneCollision;

    collisionFunctionTable[SHAPE::SPHE][SHAPE::PLAN] = spherePlaneCollision;
    //collisionFunctionTable[SHAPE::PLAN][SHAPE::SPHE] = spherePlaneCollision;

    collisionFunctionTable[SHAPE::CIRC][SHAPE::CIRC] = circleCircleCollision;

    collisionFunctionTable[SHAPE::SPHE][SHAPE::SPHE] = sphereSphereCollision;
}
