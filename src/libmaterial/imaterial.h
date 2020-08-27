#ifndef IMATERIAL_HPP
#define IMATERIAL_HPP

#include "../types.h"
namespace sdfibm{

class IMaterial
{
private:
    // density
    scalar m_rho;
public:
    IMaterial(scalar rho = 1.0)
    {
        m_rho = rho;
    }
    const scalar& getRho() const {return m_rho;}
};

class MaterialDefault : public IMaterial
{
public:
    MaterialDefault():IMaterial() {}
};

}
#endif
