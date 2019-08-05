#ifndef MOTION000010_H
#define MOTION000010_H

#include "imotion.h"

class Motion000010:public IMotion, _creator<Motion000010>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion000010(const dictionary&){}
    virtual ~Motion000010() override final {}
    TYPENAME("Motion000010")
    virtual string description() const override {return "only free to rotate in y";}
};

void Motion000010::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
    omega[0] = 0; omega[2] = 0;
}

#endif // MOTION000010_H
