#ifndef MOTION010010_H
#define MOTION010010_H

#include "imotion.h"
namespace sdfibm{

class Motion010010:public IMotion, _creator<Motion010010>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion010010(const dictionary&){}
    virtual ~Motion010010() override final {}
    TYPENAME("Motion010010")
    virtual string description() const override {return "only free to translate & rotate in y";}
};

void Motion010010::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity[0] = 0; velocity[2] = 0;
    omega[0] = 0; omega[2] = 0;
}

}
#endif // MOTION010010_H
