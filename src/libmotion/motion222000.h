#ifndef MOTION222000_H
#define MOTION222000_H

#include "imotion.h"
namespace sdfibm{

class Motion222000:public IMotion, _creator<Motion222000>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion222000(const dictionary& para)
    {
        try {
            m_u = Foam::readScalar(para.lookup("u"));
            m_v = Foam::readScalar(para.lookup("v"));
            m_w = Foam::readScalar(para.lookup("w"));
        } catch (const std::exception& e) {
            std::cout << "Problem in creating Motion222000 (aka linear motion)"
                      << e.what() << std::endl;
        }
    }
    virtual ~Motion222000() override final {}
    TYPENAME("Motion222000") // a linear motion
    virtual string description() const override {return "move at prescribed (u,v,w), no rotation";}
private:
    real m_u, m_v, m_w;
};

void Motion222000::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = vector(m_u, m_v, m_w);
    omega    = vector::zero;
}

}
#endif // MOTION222000_H
