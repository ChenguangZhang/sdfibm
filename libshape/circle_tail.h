#ifndef CIRCLE_TAIL_H
#define CIRCLE_TAIL_H

#include "ishape.h"

class Circle_Tail : public IShape, _shapecreator<Circle_Tail>
{
private:
    real m_radius;
    real m_ratio;   // tail-length/circle-radius to derive tail length
    real m_radiusb; // tail-thickness

    real m_radiusSQR;
    real m_radiusa;

public:
    Circle_Tail(const dictionary& para)
    {
        m_radius = Foam::readScalar(para.lookup("radius"));
        m_ratio  = Foam::readScalar(para.lookup("ratio"));
        m_radiusb = Foam::readScalar(para.lookup("thickness"));
        m_com = para.lookupOrDefault("com", vector::zero);

        m_radiusSQR = m_radius * m_radius;
        m_radiusa = (m_ratio + 1)*0.5*m_radius;
        // set inherited variables
        m_volume = M_PI*m_radiusSQR;
        m_volumeINV = 1.0/m_volume;
        real tmp = 0.5*m_volume*m_radiusSQR;
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = 2*m_radiusa;
    }
    inline real getRadius() const { return m_radius;}
    inline real getVolume() const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Circle_Tail")
    virtual std::string description() const override {return "Circle_Tail (x-y plane), r = " + std::to_string(m_radius);}

    virtual inline bool isInside(
            const vector& pworld,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        vector p = m_com + transform(pworld, shape_center, shape_orientation);
        p.z() = 0.0;

        bool b1 = _sdf_rectangle_bool(_sdf_offset(p, vector(m_radiusa, 0.0, 0.0)), m_radiusa, m_radiusb);
        bool b2 = _sdf_circle_bool_fast(vector(p.x(), p.y(), 0.0), m_radiusSQR);
        return _sdf_union({b1, b2});
    }

    virtual inline real signedDistance(
            const vector& pworld,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        vector p = m_com + transform(pworld, shape_center, shape_orientation);
        p.z() = 0.0;

        real d1 = _sdf_circle_real(p, m_radius);
        real d2 = _sdf_rectangle_real(_sdf_offset(p, vector(m_radiusa, 0.0, 0.0)),
                            m_radiusa, m_radiusb);

        return _sdf_filter(_sdf_union({d1, d2}));
    }
};
#endif
