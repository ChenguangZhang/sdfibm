#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "ishape.h"
namespace sdfibm{

class Ellipse : public IShape, _shapecreator<Ellipse>
{
private:
    real m_radiusa, m_radiusb;
    real m_radiusaSQRINV, m_radiusbSQRINV;

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

        real tmp = 0.25*m_volume*(m_radiusa*m_radiusa + m_radiusb*m_radiusb);
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp;
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = std::max(m_radiusa, m_radiusb);
    }
    inline real getRadiusa() const { return m_radiusa;}
    inline real getRadiusb() const { return m_radiusb;}
    inline real getVolume()  const { return m_volume;}

    // implement interface
    SHAPETYPENAME("Ellipse")
    virtual std::string description() const override {return "ellipse (x-y plane), [ra, rb] = " + std::to_string(m_radiusa) + ", " + std::to_string(m_radiusb);}

    virtual inline bool isInside(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        return _sdf_ellipse_bool_fast(
                m_com + transform(p, shape_center, shape_orientation),
                m_radiusaSQRINV,
                m_radiusbSQRINV
        );
    }

    virtual inline real signedDistance(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        return _sdf_filter(
            _sdf_ellipse_real_fast(
                m_com + transform(p,shape_center,shape_orientation),
                m_radiusaSQRINV,
                m_radiusbSQRINV
           )
        );
    }
};

}
#endif
