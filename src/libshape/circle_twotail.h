#ifndef CIRCLE_TWOTAIL_H
#define CIRCLE_TWOTAIL_H

#include "ishape.h"
namespace sdfibm{

class Circle_TwoTail : public IShape, _shapecreator<Circle_TwoTail>
{
private:
    real m_radius;
    real m_ratio;   // tail-length/circle-radius to derive tail length
    real m_radiusb; // tail-thickness

    real m_radiusSQR;
    real m_radiusa;

public:
    Circle_TwoTail(const dictionary& para)
    {
        m_radius = Foam::readScalar(para.lookup("radius"));
        m_ratio  = Foam::readScalar(para.lookup("ratio"));
        m_radiusb = Foam::readScalar(para.lookup("thickness"))*0.5;
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
    SHAPETYPENAME("Circle_TwoTail")
    virtual std::string description() const override {return "Circle_TwoTail (x-y plane), r = " + std::to_string(m_radius);}

    virtual inline bool isInside(
            const vector& pworld,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        vector p = m_com + transform(pworld, shape_center, shape_orientation);
        p.z() = 0.0;

        bool dc = _sdf_circle_bool_fast(p, m_radiusSQR);
        bool d1 = _sdf_rectangle_bool(_sdf_offset(_sdf_rot30(p), vector(m_radiusa, 0.0, 0.0)),
                    m_radiusa, m_radiusb);
        bool d2 = _sdf_rectangle_bool(_sdf_offset(_sdf_rot30(_sdf_flipy(p)), vector(m_radiusa, 0.0, 0.0)),
                    m_radiusa, m_radiusb);

        return _sdf_union({dc, d1, d2});
    }

    virtual inline real signedDistance(
            const vector& pworld,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        vector p = m_com + transform(pworld, shape_center, shape_orientation);
        p.z() = 0.0;

        real dc = _sdf_circle_real(p, m_radius);
        real d1 = _sdf_rectangle_real(_sdf_offset(_sdf_rot30(p), vector(m_radiusa, 0.0, 0.0)),
                    m_radiusa, m_radiusb);
        real d2 = _sdf_rectangle_real(_sdf_offset(_sdf_rot30(_sdf_flipy(p)), vector(m_radiusa, 0.0, 0.0)),
                    m_radiusa, m_radiusb);

        return _sdf_filter(
                    _sdf_union({dc, d1, d2})
        );
    }
};

}
#endif
