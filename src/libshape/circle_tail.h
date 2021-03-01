#ifndef CIRCLE_TAIL_H
#define CIRCLE_TAIL_H

#include "ishape.h"
namespace sdfibm{

class Circle_Tail : public IShape, _shapecreator<Circle_Tail>
{
private:
    scalar m_radius;
    scalar m_ratio;   // tail-length/circle-radius to derive tail length
    scalar m_radiusb; // tail-thickness

    scalar m_radiusSQR;
    scalar m_radiusa;

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
        scalar tmp = 0.5*m_volume*m_radiusSQR;
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = 2*m_radiusa;
    }
    inline scalar getRadius() const { return m_radius;}
    inline scalar getVolume() const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Circle_Tail")
    virtual std::string description() const override {return "Circle_Tail (x-y plane), r = " + std::to_string(m_radius);}

    virtual inline bool isInside(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;

        bool b1 = sdf::rectangle_bool(sdf::offset(p2d, vector(m_radiusa, 0.0, 0.0)), m_radiusa, m_radiusb);
        bool b2 = sdf::circle_bool_fast(vector(p2d.x(), p2d.y(), 0.0), m_radiusSQR);

        return sdf::U({b1, b2});
    }

    virtual inline scalar signedDistance(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;

        scalar d1 = sdf::circle(p2d, m_radius);
        scalar d2 = sdf::rectangle(sdf::offset(p2d, vector(m_radiusa, 0.0, 0.0)),
                            m_radiusa, m_radiusb);

        return sdf::filter(sdf::U({d1, d2}));
    }
};

}
#endif
