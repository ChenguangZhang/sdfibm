#ifndef MOTION000001_H
#define MOTION000001_H

#include "imotion.h"
namespace sdfibm{

class Motion000001:public IMotion, _creator<Motion000001>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion000001(const dictionary&){}
    virtual ~Motion000001() override final {}
    TYPENAME("Motion000001")
    virtual string description() const override {return "only free to rotate in z";}
};

void Motion000001::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
    omega[0] = 0; omega[1] = 0;
}

}
#endif // MOTION000001_H
