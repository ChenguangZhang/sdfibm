#ifndef MOTION000111_H
#define MOTION000111_H

#include "imotion.h"

class Motion000111:public IMotion, _creator<Motion000111>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion000111(const dictionary&){}
    virtual ~Motion000111() override final {}
    TYPENAME("Motion000111")
    virtual string description() const override {return "only free to rotate in xyz";}
};

void Motion000111::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
}

#endif // MOTION000111_H
