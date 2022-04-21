#ifndef SOLID_H_INCLUDED
#define SOLID_H_INCLUDED

#include "types.h"
//#include "./libshape/vof.h"
#include "./libshape/ishape.h"
#include "./libmotion/imotion.h"
#include "./libmaterial/imaterial.h"
#include "utils.h"
#include "fvc.H"

namespace sdfibm {

class Solid
{
protected:
    label id;     // integer label
    label hostid; // id of the mesh cell that host the solid center

    // state
    vector center;
    quaternion orientation;

    // state's 1st time derivative
    vector velocity;
    vector omega;

    // state's 2nd time derivative
    vector force;
    vector torque;

    // fluid only force
    vector fluid_force;
    vector fluid_torque;

    // fluid only force from last time step
    vector fluid_force_old;
    vector fluid_torque_old;

    // property pointers
    IMotion*   ptr_motion;
    IShape*    ptr_shape;
    IMaterial* ptr_material;

    // object properties = shape x material
    scalar mass;
    scalar mass_inv;
    tensor moi_inv;

public:
    // getters
    inline         label getID()       const {return id;}
    inline const vector& getCenter()   const {return center;}
    inline const vector& getVelocity() const {return velocity;}
    inline const vector& getOmega()    const {return omega;}
    inline const vector& getForce()    const {return force;}
    inline const vector& getTorque()   const {return torque;}
    inline const quaternion& getOrientation() const {return orientation;}

    inline const vector& getFluidForce()  {return fluid_force; }
    inline const vector& getFluidTorque() {return fluid_torque; }

    // setters
    inline void setCenter  (const vector& c){ center = c; }
    inline void setVelocity(const vector& v){ velocity = v; }
    inline void setOmega   (const vector& o){ omega = o; }
    inline void setForce   (const vector& f){ force = f; }
    inline void setTorque  (const vector& t){ torque= t; }
    inline void setOrientation(const vector& angles)
    {
        // from a fixed order of Euler angles
        orientation = quaternion(quaternion::XYZ, angles);
    }
    inline void setFluidForceOld(const vector& f)  {fluid_force_old = f; }
    inline void setFluidTorqueOld(const vector& t) {fluid_torque_old= t; }

    inline IMotion*   getMotion()   const {return ptr_motion;  }
    inline IShape*    getShape()    const {return ptr_shape;   }
    inline IMaterial* getMaterial() const {return ptr_material;}
    inline scalar     getRadiusB()  const {return ptr_shape->getRadiusB();}

    // done getters and setters
    // set solid material and shape
    inline void setMaterialAndShape(IMaterial* material,
                                    IShape*    shape)
    {
        ptr_material = material;
        ptr_shape    = shape;

        // update composite properties
        scalar rho = ptr_material->getRho();
        mass     = ptr_shape->m_volume    * rho;
        mass_inv = ptr_shape->m_volumeINV / rho;
        moi_inv  = ptr_shape->m_moiINV    / rho;
    }

    // control of solid motion
    inline void setMotion(IMotion* solid_motion) { ptr_motion = solid_motion; }
    inline void unsetMotion()                    { ptr_motion = nullptr; }

    /* Check if point p is inside the solid
     *
     * 20/02/2021 All geometric calculations are delegated to shape component of the solid. Before
     * that however, the solid should transfer form the point to the object space of its shape component.
     * TODO: for circle/sphere shapes, the quaternion transformation can be omitted.
     *
     * @param p World coordinate of a point (e.g., a mesh vertex or cell center
     * @return true if point inside solid, false if outside of the solid
     */
    inline bool isInside(const vector& p) const
    {
        return ptr_shape->isInside(
            Foam::conjugate(orientation).transform(p - center)
        );
    }

    /* Find point p's sdf w.r.t. the solid
     *
     * Note: see that of isInside function
     *
     * @param p World coordinate of a point (e.g., a mesh vertex or cell center
     * @return phi the sdf of the point
     */
    inline scalar signedDistance(const vector& p) const
    {
        return ptr_shape->signedDistance(
            Foam::conjugate(orientation).transform(p - center)
        );
    }

    inline vector evalPointVelocity(const vector& p) const
    {
        return velocity + (omega^(p - center));
    }

    inline void addAcceleration(const vector& acc)
    {
        // dedicated to body force which needs the (priviate) mass info
        force += mass*acc;
    }
    inline void clearForceAndTorque()
    {
        force  = vector::zero;
        torque = vector::zero;
    }
    inline void setFluidForceAndTorque(const vector& fluid_force_in, const vector& fluid_torque_in)
    {
        fluid_force  = fluid_force_in;
        fluid_torque = fluid_torque_in;
    }
    inline void storeOldFluidForce()
    {
        fluid_force_old  = fluid_force;
        fluid_torque_old = fluid_torque;
    }
    inline void addMidFluidForceAndTorque()
    {
        force  += (1.5*fluid_force  - 0.5*fluid_force_old);
        torque += (1.5*fluid_torque - 0.5*fluid_torque_old);
    }
    inline void addForceAndTorque(const vector& external_force, const vector& external_torque)
    {
        force  += external_force;
        torque += external_torque;
    }

    void move(const scalar& time, const scalar& dt)
    {
        // motion = velocity & omega
        // temporarily store motion at time n
        vector velocity_old = velocity;
        vector omega_old    = omega;

        // update motion to n+1 using force at time n+1/2
        velocity += force*mass_inv*dt;        // velocity updated to t + dt
        tensor R = orientation.R();
        tensor moi_inv_world = R & moi_inv & R.T();
        omega += (moi_inv_world & torque)*dt; // omega updated to t + dt

        // constrain motion
        if(ptr_motion != nullptr)
            ptr_motion->constraint(time, velocity, omega);

        // position & orientation updated AFTER constraint
        center      += 0.5*(velocity + velocity_old)*dt;
        orientation += 0.5*quaternion(0.5*(omega + omega_old))*orientation*dt;
        orientation.normalize(); // no need to normalize every step, but cheap anyway
    }

    Solid(label solid_id,
          const vector& solid_center,
          const quaternion& solid_quaternion):
        id      (solid_id),
        center  (solid_center),
        orientation (solid_quaternion)
    {
        velocity = vector::zero;
        omega    = vector::zero;
        force    = vector::zero;
        torque   = vector::zero;

        fluid_force  = vector::zero;
        fluid_torque = vector::zero;
        fluid_force_old  = vector::zero;
        fluid_torque_old = vector::zero;

        ptr_motion   = nullptr;
        ptr_shape    = nullptr;
        ptr_material = nullptr;
        mass = 0;
        mass_inv = 0;
        hostid =-1;
    }

    ~Solid(){}
    friend std::ostream& operator<<(std::ostream& os, const Solid& s);
    friend void write2D(std::ostream& os, const Solid& s);
};

}
#endif
