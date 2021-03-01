#ifndef CIRCLE_TWOTAIL_H
#define CIRCLE_TWOTAIL_H

#include "ishape.h"
namespace sdfibm{

class Circle_TwoTail : public IShape, _shapecreator<Circle_TwoTail>
{
private:
    scalar m_radius;
    scalar m_ratio;   // tail-length/circle-radius to derive tail length
    scalar m_radiusb; // tail-thickness

    scalar m_radiusSQR;
    scalar m_radiusa;

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
        scalar tmp = 0.5*m_volume*m_radiusSQR;
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = 2*m_radiusa;
    }
    inline scalar getRadius() const { return m_radius;}
    inline scalar getVolume() const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Circle_TwoTail")
    virtual std::string description() const override {return "Circle_TwoTail (x-y plane), r = " + std::to_string(m_radius);}

    virtual inline bool isInside(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;

        bool dc = sdf::circle_bool_fast(p2d, m_radiusSQR);
        bool d1 = sdf::rectangle_bool(sdf::offset(sdf::rot30(            p2d), vector(m_radiusa, 0, 0)),
                    m_radiusa, m_radiusb);
        bool d2 = sdf::rectangle_bool(sdf::offset(sdf::rot30(sdf::flipy(p2d)), vector(m_radiusa, 0, 0)),
                    m_radiusa, m_radiusb);

        return sdf::U({dc, d1, d2});
    }

    virtual inline scalar signedDistance(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;

        scalar dc = sdf::circle(p2d, m_radius);
        scalar d1 = sdf::rectangle(sdf::offset(sdf::rot30(            p2d), vector(m_radiusa, 0, 0)),
                    m_radiusa, m_radiusb);
        scalar d2 = sdf::rectangle(sdf::offset(sdf::rot30(sdf::flipy(p2d)), vector(m_radiusa, 0, 0)),
                    m_radiusa, m_radiusb);

        return sdf::filter(sdf::U({dc, d1, d2}));
    }
};

}
#endif
