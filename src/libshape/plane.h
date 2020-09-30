#ifndef PLANE_H
#define PLANE_H

#include "ishape.h"
namespace sdfibm {

class Plane : public IShape, _shapecreator<Plane>
{
private:
    scalar m_radius;
    scalar m_radiusSQR;

public:
    // const static int shape_id = SHAPE::CIRC;
    Plane(const dictionary& para) {}

    // typename and description
    SHAPETYPENAME("Plane")
    virtual std::string description() const override {return "Plane (x-z plane)";}

    // implement interface
    virtual inline bool isInside(
            const vector& pworld,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        vector p = transform(pworld, shape_center, shape_orientation);
        return p.y() < 0;
    }

    virtual inline scalar signedDistance(
            const vector& pworld,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        vector p = transform(pworld, shape_center, shape_orientation);
        return p.y();
    }
};

}
#endif
