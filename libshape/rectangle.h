#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "ishape.h"

class Rectangle : public IShape, _shapecreator<Rectangle>
{
private:
    real m_radiusa, m_radiusb;
    real m_radiusaSQR;
    real m_radiusbSQR;

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

        real tmp = 1.0/3.0*m_volume*(m_radiusaSQR + m_radiusbSQR);
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = std::max(m_radiusa, m_radiusb);
    }
    inline real getRadiusa() const { return m_radiusa;}
    inline real getRadiusb() const { return m_radiusb;}
    inline real getVolume() const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Rectangle")
    virtual std::string description() const override {return "rectangle (x-y plane), [ra, rb] = " + std::to_string(m_radiusa) + ", " + std::to_string(m_radiusb);}

    virtual inline bool isInside(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        return _sdf_rectangle_bool(
                m_com + transform(p, shape_center, shape_orientation),
                m_radiusa,
                m_radiusb);
    }

    virtual inline real signedDistance(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        return _sdf_filter(_sdf_rectangle_real(
                    m_com + transform(p, shape_center, shape_orientation),
                    m_radiusa,
                    m_radiusb)
        );
    }
};
#endif
