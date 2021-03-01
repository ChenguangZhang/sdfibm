#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "ishape.h"
namespace sdfibm{

class Rectangle : public IShape, _shapecreator<Rectangle>
{
private:
    scalar m_radiusa, m_radiusb;
    scalar m_radiusaSQR;
    scalar m_radiusbSQR;

public:
    // const static int shape_id = SHAPE::CIRC;
    Rectangle(const dictionary& para)
    {
        m_radiusa =  Foam::readScalar(para.lookup("radiusa"));
        m_radiusb =  Foam::readScalar(para.lookup("radiusb"));
        m_com = para.lookupOrDefault("com", vector::zero);
        m_radiusaSQR = m_radiusa * m_radiusa;
        m_radiusbSQR = m_radiusb * m_radiusb;

        // set inherited variables
        m_volume = 4.0*m_radiusa*m_radiusb;
        m_volumeINV = 1.0/m_volume;

        scalar tmp = 1.0/3.0*m_volume*(m_radiusaSQR + m_radiusbSQR);
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = std::max(m_radiusa, m_radiusb);
    }
    inline scalar getRadiusa() const { return m_radiusa;}
    inline scalar getRadiusb() const { return m_radiusb;}
    inline scalar getVolume() const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Rectangle")
    virtual std::string description() const override {return "rectangle (x-y plane), [ra, rb] = " + std::to_string(m_radiusa) + ", " + std::to_string(m_radiusb);}

    virtual inline bool isInside(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;
        return sdf::rectangle_bool(
                p2d,
                m_radiusa,
                m_radiusb);
    }

    virtual inline scalar signedDistance(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;
        return sdf::filter(sdf::rectangle(
                    p2d,
                    m_radiusa,
                    m_radiusb)
        );
    }
};

}
#endif
