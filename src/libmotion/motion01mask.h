#ifndef MOTION01MASK_H
#define MOTION01MASK_H

#include "imotion.h"
namespace sdfibm{

class Motion01Mask:public IMotion, _creator<Motion01Mask>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    virtual ~Motion01Mask() override final {}
    TYPENAME("Motion01Mask")
    virtual string description() const override {return "general (0|1){6} motion mask";}

    Motion01Mask(const dictionary& para)
    {
        vmask = vector::zero;
        omask = vector::zero;

        std::string mask = Foam::word(para.lookup("mask"));

        // sanity check
        if(mask[0]!='b')
        {
            std::cout << "Mask starts with b\n";
            std::exit(1);
        }
        if(mask.size() != 7)
        {
            std::cout << "Mask size != 7\n";
            std::exit(1);
        }
        for(int i = 1; i<=6; i++)
        {
            if(mask[i]!='0' && mask[i]!='1')
            {
                std::cout << "Mask ele != 0 or 1\n";
                std::exit(1);
            }
        }

        vmask[0] = (mask[1]=='0') ? 0 : 1;
        vmask[1] = (mask[2]=='0') ? 0 : 1;
        vmask[2] = (mask[3]=='0') ? 0 : 1;
        omask[0] = (mask[4]=='0') ? 0 : 1;
        omask[1] = (mask[5]=='0') ? 0 : 1;
        omask[2] = (mask[6]=='0') ? 0 : 1;
    }

private:
    vector vmask, omask;
};

void Motion01Mask::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = cmptMultiply(velocity, vmask);
    omega    = cmptMultiply(omega, omask);
}

}
#endif // MOTION01MASK
