#ifndef MOTION111111_H
#define MOTION111111_H

#include "imotion.h"
namespace sdfibm{

class Motion111111:public IMotion, _creator<Motion111111>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion111111(const dictionary&){}
    virtual ~Motion111111() override final {}
    TYPENAME("Motion111111")
    virtual string description() const override {return "free to move and rotate";}
};

void Motion111111::constraint(const scalar &time, vector &velocity, vector &omega)
{
    // free motion, do nothing
}

}
#endif // MOTION111111_H
