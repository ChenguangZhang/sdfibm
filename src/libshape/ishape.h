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
    struct Transformation {vector t; quaternion q;};
    const static int m_id = -1;
    scalar m_radiusB   {0.0}; // radius of bounding sphere
    scalar m_volume    {0.0};
    scalar m_volumeINV {0.0}; 
    vector m_com    {vector::zero};
    tensor m_moi    {tensor::I}; // in principal frame, diagonal
    tensor m_moiINV {tensor::I};
    bool finite {true};

public:
    SHAPETYPENAME("IShape");
    
    inline static vector world2local(const vector& p, const Transformation& tr)
    {
        return Foam::conjugate(tr.q).transform(p - tr.t);
    }

    virtual int    getShapeID() const {return m_id;}
    virtual scalar getRadiusB() const {return m_radiusB;}

    bool phi01(const vector& p, const Transformation& tr) { return isInside      (world2local(p, tr)); }
    scalar phi(const vector& p, const Transformation& tr) { return signedDistance(world2local(p, tr)); }

    virtual std::string description() const = 0;

    virtual ~IShape(){}

private:
    virtual bool   isInside      (const vector& p) const = 0; // local coordinate
    virtual scalar signedDistance(const vector& p) const = 0; // local coordinate
};

}
#endif
