#ifndef IMATERIAL_HPP
#define IMATERIAL_HPP

#include "types.h"

class IMaterial
{
private:
    // density
    real m_rho;
public:
    IMaterial(real rho = 1.0)
    {
        m_rho = rho;
    }
    const real& getRho() const {return m_rho;}
};

class MaterialDefault : public IMaterial
{
public:
    MaterialDefault():IMaterial() {}
};

#endif
