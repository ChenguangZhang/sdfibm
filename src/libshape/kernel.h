#include "../types.h"

#define THRESHOLD (1.0e-6)

inline int signt(const real& val) {
    if(std::fabs(val) < 1.0e-6)
        return 1;
    return (0.0 < val) - (val < 0.0);
}

real squareFraction(const real& ll_in, const real& lr_in, const real& ur_in, const real& ul_in)
{
    // accept four dimensional distances
    real ur = std::fabs(ur_in) < THRESHOLD ? THRESHOLD : ur_in;
    real ul = std::fabs(ul_in) < THRESHOLD ? THRESHOLD : ul_in;
    real lr = std::fabs(lr_in) < THRESHOLD ? THRESHOLD : lr_in;
    real ll = std::fabs(ll_in) < THRESHOLD ? THRESHOLD : ll_in;
    int sign_sum = signt(ur) + signt(ul) + signt(lr) + signt(ll);
    // case 1, 4 is out sign_sum =  4
    if(sign_sum == 4) return 0;
    // case 5, 0 is out sign_sum = -4
    if(sign_sum ==-4) return 1;

    // case 2, 3 is out sign_sum =  2
    if(sign_sum == 2)
    {
        real epsilon1(0.0), epsilon2(0.0);
        if(ul < 0) // ul is outside
        {
            epsilon1 = ul/(ul - ur);
            epsilon2 = ul/(ul - ll);
        }
        if(ur < 0) // ur is outside
        {
            epsilon1 = ur/(ur - ul);
            epsilon2 = ur/(ur - lr);
        }
        if(ll < 0) // ur is outside
        {
            epsilon1 = ll/(ll - ul);
            epsilon2 = ll/(ll - lr);
        }
        if(lr < 0) // ur is outside
        {
            epsilon1 = lr/(lr - ll);
            epsilon2 = lr/(lr - ur);
        }
        return 0.5*(epsilon1*epsilon2);
    }
    
    // case 4, 1 is out sign_sum = -2
    if(sign_sum == -2)
    {
        real epsilon1(0.0), epsilon2(0.0);
        if(ul > 0) // ul is outside
        {
            epsilon1 = ul/(ul - ur);
            epsilon2 = ul/(ul - ll);
        }
        if(ur > 0) // ur is outside
        {
            epsilon1 = ur/(ur - ul);
            epsilon2 = ur/(ur - lr);
        }
        if(ll > 0) // ur is outside
        {
            epsilon1 = ll/(ll - ul);
            epsilon2 = ll/(ll - lr);
        }
        if(lr > 0) // ur is outside
        {
            epsilon1 = lr/(lr - ll);
            epsilon2 = lr/(lr - ur);
        }
        return 1.0 - 0.5*(epsilon1*epsilon2);
    }
    // case 3, 2 is out sign_sum =  0
    if(sign_sum == 0)
    {
        real positive_sum = 0.0;
        if(ur > 0) positive_sum += ur;
        if(ul > 0) positive_sum += ul;
        if(lr > 0) positive_sum += lr;
        if(ll > 0) positive_sum += ll;
        return 1.0 - positive_sum/(std::fabs(ur)+
                                   std::fabs(ul)+
                                   std::fabs(lr)+
                                   std::fabs(ll));
    }
    return -1;
}

real cubeFraction(const double (&d)[8])
{
    static const vector CUBE_FNORMAL[6] = {
        vector( 0, 0, 1),
        vector( 0,-1, 0),
        vector(-1, 0, 0),
        vector( 0, 1, 0),
        vector( 1, 0, 0),
        vector( 0, 0,-1)
    };
    // unlike in vof.cpp, here unitx/y/z are dimensionless. I did some math Sept. 11 which eventually removes delta and
    // cell center from the calculation
    static const vector unitx(1, 0, 0);
    static const vector unity(0, 1, 0);
    static const vector unitz(0, 0, 1);

    real area_fraction[6];
    area_fraction[0] = squareFraction(d[0], d[1], d[2], d[3]);
    area_fraction[1] = squareFraction(d[7], d[6], d[2], d[3]);
    area_fraction[2] = squareFraction(d[5], d[1], d[2], d[6]);
    area_fraction[3] = squareFraction(d[4], d[5], d[1], d[0]);
    area_fraction[4] = squareFraction(d[4], d[0], d[3], d[7]);
    area_fraction[5] = squareFraction(d[4], d[5], d[6], d[7]);

    // find and use intersection point as the top vertex of all the prisms
    real r;
    vector p_I;
    if(d[0] * d[6] < 0)
    {
       r = std::min(1.0, std::fabs(d[0])/(std::fabs(d[0]) + std::fabs(d[6])));
       p_I = r * (unitx + unity + unitz);
    }
    else if(d[1] * d[7] < 0)
    {
       r = std::min(1.0, std::fabs(d[1])/(std::fabs(d[1]) + std::fabs(d[7])));
       p_I = (1 - r) * (unitx) + r * (unity + unitz);
    }
    else if(d[2] * d[4] < 0)
    {
       r = std::min(1.0, std::fabs(d[2])/(std::fabs(d[2]) + std::fabs(d[4])));
       p_I = (1 - r) * (unitx + unity) + r * unitz;
    }
    else if(d[3] * d[5] < 0)
    {
       r = std::min(1.0, std::fabs(d[3])/(std::fabs(d[3]) + std::fabs(d[5])));
       p_I = (1 - r) * unity + r * (unitx + unitz);
    }

    // project the cell center on the plane, with this new point
    real vof = 0.0;
    for(label i = 0; i < 6; ++i)
    {
       vof = vof + area_fraction[i] * (CUBE_FNORMAL[i] & (p_I - 0.5*vector::one + 0.5*CUBE_FNORMAL[i]) );
    }
    return vof/3.0;
//        Real nxdeltaabs = std::fabs(d[1] - d[0]);
//        Real nydeltaabs = std::fabs(d[3] - d[0]);
//        Real nzdeltaabs = std::fabs(d[4] - d[0]);
//        Real Acut = 0.0;
//        if(nxdeltaabs > nydeltaabs && nxdeltaabs > nzdeltaabs)
//            Acut = std::fabs(area_fraction[5] - area_fraction[0])/nxdeltaabs;
//        if(nydeltaabs > nxdeltaabs && nydeltaabs > nzdeltaabs)
//            Acut = std::fabs(area_fraction[3] - area_fraction[1])/nydeltaabs;
//        else
//            Acut = std::fabs(area_fraction[4] - area_fraction[2])/nzdeltaabs;
    //    return Acut;

}


