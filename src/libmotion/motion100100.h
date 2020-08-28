#ifndef MOTION100100_H
#define MOTION100100_H

#include "imotion.h"
namespace sdfibm{

class Motion100100:public IMotion, _creator<Motion100100>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion100100(const dictionary&){}
    virtual ~Motion100100() override final {}
    TYPENAME("Motion100100")
    virtual string description() const override {return "only free to translate & rotate in x";}
};

void Motion100100::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity[1] = 0; velocity[2] = 0;
    omega[1] = 0; omega[2] = 0;
}

}
#endif // MOTION100100_H
