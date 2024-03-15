#ifndef ISHAPE_H
#define ISHAPE_H

#include "../types.h"
#include "./sdf/sdf.h"
#include <algorithm>
namespace sdfibm {

#define SHAPETYPENAME(name) \
    static std::string typeName() {return name;} \
    static bool added; \
    virtual std::string getTypeName() const {return name;}

class IShape;

template <typename T>
class _shapecreator
{
public:
    static IShape* create(const dictionary& para)
    {
        return new T(para);
    }
};

class IShape
{
public:
    struct Transformation {const quaternion& q; const vector& t;};
    const static int m_id = -1;
    scalar m_radiusB;             // radius of bounding sphere
    scalar m_volume, m_volumeINV; // volume
    vector m_com;                 // center of mass
    tensor m_moi, m_moiINV;       // moi in principal frame (nonzero only along diagonal)
    bool finite;

public:
    IShape()
    {
        m_radiusB   = 0.0;
        m_volume    = 0.0;
        m_volumeINV = 0.0;
        m_com       = vector(0.0, 0.0, 0.0);

        m_moi    = tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        m_moiINV = tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        finite = true;
    }

    SHAPETYPENAME("IShape");
    
    inline vector world2local(const vector& p, const Transformation& tr) const
    {
        return Foam::conjugate(tr.q).transform(p - tr.t);
    }

    virtual int    getShapeID() const {return m_id;}
    virtual scalar getRadiusB() const {return m_radiusB;}

    virtual bool phi01(const vector& p, const Transformation& tr) { return isInside      (world2local(p, tr)); }
    virtual scalar phi(const vector& p, const Transformation& tr) { return signedDistance(world2local(p, tr)); }

    virtual std::string description   ()                const = 0;

    virtual ~IShape(){}

private:
    virtual bool   isInside      (const vector& p) const = 0; // local coordinate
    virtual scalar signedDistance(const vector& p) const = 0; // local coordinate
};

}
#endif
