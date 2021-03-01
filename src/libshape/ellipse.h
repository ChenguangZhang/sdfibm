#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "ishape.h"
namespace sdfibm {

class Ellipse : public IShape, _shapecreator<Ellipse>
{
private:
    scalar m_radiusa, m_radiusb;
    scalar m_radiusaSQRINV, m_radiusbSQRINV;

public:
    Ellipse(const dictionary& para)
    {
        m_radiusa = Foam::readScalar(para.lookup("radiusa"));
        m_radiusb = Foam::readScalar(para.lookup("radiusb"));
        m_com = para.lookupOrDefault("com", vector::zero);

        m_radiusaSQRINV = 1.0/(m_radiusa * m_radiusa);
        m_radiusbSQRINV = 1.0/(m_radiusb * m_radiusb);

        // set inherited variables
        m_volume = M_PI*m_radiusa*m_radiusb;
        m_volumeINV = 1.0/m_volume;

        scalar tmp = 0.25*m_volume*(m_radiusa*m_radiusa + m_radiusb*m_radiusb);
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = std::max(m_radiusa, m_radiusb);
    }
    inline scalar getRadiusa() const { return m_radiusa;}
    inline scalar getRadiusb() const { return m_radiusb;}
    inline scalar getVolume()  const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Ellipse")
    virtual std::string description() const override {return "ellipse (x-y plane), [ra, rb] = " + std::to_string(m_radiusa) + ", " + std::to_string(m_radiusb);}

    virtual inline bool isInside(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;

        return sdf::ellipse_bool_fast(p2d, m_radiusaSQRINV, m_radiusbSQRINV);
    }

    virtual inline scalar signedDistance(const vector& p) const override
    {
        vector p2d = m_com + p; p2d.z() = 0.0;

        return sdf::filter(
            sdf::ellipse(p2d, m_radiusaSQRINV, m_radiusbSQRINV)
        );
    }
};

}
#endif
