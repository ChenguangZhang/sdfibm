#ifndef MOTION001001_H
#define MOTION001001_H

#include "imotion.h"
namespace sdfibm{

class Motion001001:public IMotion, _creator<Motion001001>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion001001(const dictionary&){}
    virtual ~Motion001001() override final {}
    TYPENAME("Motion001001")
    virtual string description() const override {return "only free to translate & rotate in z";}
};

void Motion001001::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity[0] = 0; velocity[1] = 0;
    omega[0] = 0; omega[1] = 0;
}

}
#endif // MOTION001001_H
