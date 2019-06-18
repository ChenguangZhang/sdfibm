#ifndef SPHERE_H
#define SPHERE_H

#include "ishape.h"

class Sphere : public IShape, _shapecreator<Sphere>
{
private:
    real m_radius;
    real m_radiusSQR;

public:
    Sphere(const dictionary& para)
    {
        m_radius = Foam::readScalar(para.lookup("radius"));
        m_com = para.lookupOrDefault("com", vector::zero);

        m_radiusSQR = m_radius * m_radius;

        m_volume = 4.0/3.0*M_PI*m_radiusSQR*m_radius;
        m_volumeINV = 1.0/m_volume;

        real tmp = 0.4*m_volume*m_radiusSQR;
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;

        m_moiINV = Foam::inv(m_moi);
        m_radiusB = m_radius;
    }

    inline real getRadius() const { return m_radius;}
    inline real getVolume() const { return m_volume;}

    SHAPETYPENAME("Sphere")
    virtual std::string description() const override {return "sphere, r = " + std::to_string(m_radius);}

    virtual inline bool isInside(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const
    {
        return _sdf_circle_bool_fast(m_com + p - shape_center, m_radiusSQR);
    }

    virtual inline real signedDistance(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const
    {
        return _sdf_filter(_sdf_circle_real(m_com + p - shape_center, m_radius));
    }

};

#endif
