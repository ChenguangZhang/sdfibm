#ifndef BOX_H
#define BOX_H

#include "ishape.h"
namespace sdfibm{

class Box : public IShape, _shapecreator<Box>
{
private:
    real m_radiusa, m_radiusb, m_radiusc;

public:
    // const static int shape_id = SHAPE::CIRC;
    Box(const dictionary& para)
    {
        m_radiusa =  Foam::readScalar(para.lookup("radiusa"));
        m_radiusb =  Foam::readScalar(para.lookup("radiusb"));
        m_radiusc =  Foam::readScalar(para.lookup("radiusc"));
        m_com = para.lookupOrDefault("com", vector::zero);

        // set inherited variables
        m_volume = 8.0*m_radiusa*m_radiusb*m_radiusc;
        m_volumeINV = 1.0/m_volume;

        m_moi[0] = 1.0/3.0*m_volume*(m_radiusb*m_radiusb + m_radiusc*m_radiusc);
        m_moi[4] = 1.0/3.0*m_volume*(m_radiusa*m_radiusa + m_radiusc*m_radiusc);
        m_moi[8] = 1.0/3.0*m_volume*(m_radiusb*m_radiusb + m_radiusa*m_radiusa);
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = std::max(std::max(m_radiusa, m_radiusb),m_radiusc);
    }
    inline real getRadiusa() const { return m_radiusa;}
    inline real getRadiusb() const { return m_radiusb;}
    inline real getRadiusc() const { return m_radiusc;}
    inline real getVolume()  const { return m_volume;}

    SHAPETYPENAME("Box")
    virtual std::string description() const override {return "Box, [ra, rb, rc] = " + std::to_string(m_radiusa) + ", " + std::to_string(m_radiusb) + ", " + std::to_string(m_radiusc);}

    // implement interface
    virtual inline bool isInside(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        return _sdf_box_bool(
                m_com + transform(p, shape_center, shape_orientation),
                m_radiusa,
                m_radiusb,
                m_radiusc);
    }

    virtual inline real signedDistance(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        return _sdf_filter(_sdf_box_real(
                    m_com + transform(p, shape_center, shape_orientation),
                    m_radiusa,
                    m_radiusb,
                    m_radiusc)
        );
    }
};

}
#endif
