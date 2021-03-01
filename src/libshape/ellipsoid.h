#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "ishape.h"
namespace sdfibm {

class Ellipsoid : public IShape, _shapecreator<Ellipsoid>
{
private:
    scalar m_radiusa, m_radiusb, m_radiusc;
    scalar m_radiusaSQRINV, m_radiusbSQRINV, m_radiuscSQRINV;

public:
    // const static int shape_id = SHAPE::CIRC;
    Ellipsoid(const dictionary& para)
    {
        m_radiusa =  Foam::readScalar(para.lookup("radiusa"));
        m_radiusb =  Foam::readScalar(para.lookup("radiusb"));
        m_radiusc =  Foam::readScalar(para.lookup("radiusc"));
        m_radiusaSQRINV = 1.0/(m_radiusa * m_radiusa);
        m_radiusbSQRINV = 1.0/(m_radiusb * m_radiusb);
        m_radiuscSQRINV = 1.0/(m_radiusc * m_radiusc);

        // set inherited variables
        m_volume = 4.0/3.0*M_PI*m_radiusa*m_radiusb*m_radiusc;
        m_volumeINV = 1.0/m_volume;

        m_moi[0] = 0.2*m_volume*(m_radiusb*m_radiusb + m_radiusc*m_radiusc);
        m_moi[4] = 0.2*m_volume*(m_radiusa*m_radiusa + m_radiusc*m_radiusc);
        m_moi[8] = 0.2*m_volume*(m_radiusa*m_radiusa + m_radiusb*m_radiusb);
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = std::max(std::max(m_radiusa, m_radiusb),m_radiusc);
    }
    inline scalar getRadiusa() const { return m_radiusa;}
    inline scalar getRadiusb() const { return m_radiusb;}
    inline scalar getRadiusc() const { return m_radiusc;}
    inline scalar getVolume() const { return m_volume;}

    SHAPETYPENAME("Ellipsoid")
    virtual std::string description() const override {return "ellipsoid, [ra, rb, rc] = " + std::to_string(m_radiusa) + ", " + std::to_string(m_radiusb) + ", " + std::to_string(m_radiusc);}

    // implement interface
    virtual inline bool isInside(const vector& p) const override
    {
        return sdf::ellipsoid_bool_fast(p, m_radiusaSQRINV, m_radiusbSQRINV, m_radiuscSQRINV);
    }

    virtual inline scalar signedDistance(const vector& p) const override
    {
        return sdf::filter(
            sdf::ellipsoid( p, m_radiusaSQRINV, m_radiusbSQRINV, m_radiuscSQRINV)
        );
    }
};

}
#endif
