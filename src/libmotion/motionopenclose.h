#ifndef MOTIONOPENCLOSE_H
#define MOTIONOPENCLOSE_H

#include "imotion.h"
namespace sdfibm {

class MotionOpenClose:public IMotion, _creator<MotionOpenClose>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    MotionOpenClose(const dictionary& para)
    {
    }
    virtual ~MotionOpenClose() override final {}
    TYPENAME("MotionOpenClose")
    virtual string description() const override {return "model the open-close operation of a gate";}
private:
    scalar getOpenCloseVelocity(const scalar& time)
    {
        scalar t = fmod(time, 5.0);
        scalar v = 0.0;
        if (t > 1 && t < 2)
            v = -1.0;
        if (t > 3 && t < 4)
            v =  1.0;
        return v;
    }
};

void MotionOpenClose::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity = vector(0,getOpenCloseVelocity(time),0);
    omega    = vector::zero;
}

}
#endif // MOTIONOPENCLOSE_H
