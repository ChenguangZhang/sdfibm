#ifndef MOTION000000_H
#define MOTION000000_H

#include "imotion.h"
namespace sdfibm{

class Motion000000:public IMotion, _creator<Motion000000>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion000000(const dictionary&){}
    virtual ~Motion000000() override final {}
    TYPENAME("Motion000000")
    virtual string description() const override {return "stationary";}
};

void Motion000000::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
    omega    = vector::zero;
}

}
#endif // MOTION000000_H
