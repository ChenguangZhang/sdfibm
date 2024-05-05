#ifndef MOTION01MASK_H
#define MOTION01MASK_H

#include "iforce.h"
namespace sdfibm::force {

class Constant : public IForce, _creator<Constant>
{
public:
    virtual Force generate(const scalar &time, const vector& position, const vector& velocity) override final;
    virtual ~Constant() override final {}
    TYPENAME("Constant")
    virtual std::string description() const override {return "const force and torque";}

    Constant(const dictionary& para)
    {
        force = para.lookup("force");
        torque = para.lookup("torque");
    }

private:
    vector force, torque;
};

IForce::Force Constant::generate(const scalar &time, const vector& position, const vector& velocity)
{
    return {force, torque};
}

}
#endif
