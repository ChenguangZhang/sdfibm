#ifndef MOTION000100_H
#define MOTION000100_H

#include "imotion.h"
namespace sdfibm{

class Motion000100:public IMotion, _creator<Motion000100>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion000100(const dictionary&){}
    virtual ~Motion000100() override final {}
    TYPENAME("Motion000100")
    virtual string description() const override {return "only free to rotate in x";}
};

void Motion000100::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
    omega[1] = 0; omega[2] = 0;
}

}
#endif // MOTION000100_H
