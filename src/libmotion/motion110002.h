#ifndef MOTION110002_H
#define MOTION110002_H

#include "imotion.h"
namespace sdfibm{

class Motion110002:public IMotion, _creator<Motion110002>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion110002(const dictionary& para)
    {
        try {
            m_period = Foam::readScalar(para.lookup("period"));
        } catch (const std::exception& e) {
            std::cout << "Problem in creating Motion110002 (aka const rotation)" << ' '
                      << e.what() << std::endl;
            std::exit(1);
        }
        m_omega = 2*M_PI/m_period;
    }
    virtual ~Motion110002() override final {}
    TYPENAME("Motion110002") // a linear motion
    virtual string description() const override {return "only rotate in z with prescribed omega";}
private:
    scalar m_period, m_omega; // m_omegaz is only the z component
};

void Motion110002::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity[2] = 0.0;
    omega    = vector::zero;
    omega[2] = m_omega;
}

}
#endif // MOTION110002_H
