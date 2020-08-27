#ifndef CHANGE_H
#define CHANGE_H

#include "ishape.h"
namespace sdfibm{

class CHANGE : public IShape, _shapecreator<CHANGE>
{
private:
    scalar m_radius; // define properties of your shape CHANGE

public:
    // const static int shape_id = SHAPE::CIRC;
    CHANGE(const dictionary& para)
    {
        // fetch needed property from para
        m_radius =  Foam::readScalar(para.lookup("radius")); // CHANGE
        m_com = para.lookupOrDefault("com", vector::zero);

        m_volume = M_PI*m_radiusSQR; // volume of your shape CHANGE
        m_volumeINV = 1.0/m_volume;
        scalar tmp = 0.5*m_volume*m_radiusSQR; // moi of your shape CHANGE
        m_moi[0] = tmp; m_moi[4] = tmp; m_moi[8] = tmp; // moi of your shape CHANGE
        m_moiINV = Foam::inv(m_moi);
        m_radiusB = m_radius; // bounding radius of your shape CHANGE
    }
    // how to get properties of your shape
    inline scalar getRadius() const { return m_radius;}
    inline scalar getVolume() const { return m_volume;}

    // typename and description
    SHAPETYPENAME("CHANGE")
    // describe your shape
    virtual std::string description() const override {return "CHANGE r = " + std::to_string(m_radius);} // CHANGE

    // implement interface
    virtual inline bool isInside(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        // tell if a point is inside or outside of your shape
        return _sdf_circle_bool_fast(m_com + vector(p.x(), p.y(),0.0)-shape_center, m_radiusSQR); // CHANGE
    }

    virtual inline scalar signedDistance(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const override
    {
        // sdf of your shape, final result must be wrapped by _sdf_filter to
        // avoid possible floating-point error
        return _sdf_filter(_sdf_circle_real(m_com + vector(p.x(), p.y(),0.0)-shape_center, m_radius)); // CHANGE
    }
};

}
#endif
