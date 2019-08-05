#ifndef MOTION000002_H
#define MOTION000002_H

#include "imotion.h"

class Motion000002:public IMotion, _creator<Motion000002>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    Motion000002(const dictionary& para)
    {
        try {
            m_period = Foam::readScalar(para.lookup("period"));
        } catch (const std::exception& e) {
            std::cout << "Problem in creating Motion000002 (aka const rotation)" << ' '
                      << e.what() << std::endl;
            std::exit(1);
        }
        m_omega = 2*M_PI/m_period;
    }
    virtual ~Motion000002() override final {}
    TYPENAME("Motion000002") // a linear motion
    virtual string description() const override {return "only rotate in z with prescribed omega";}
private:
    real m_period, m_omega; // m_omegaz is only the z component
};

void Motion000002::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
    omega    = vector::zero;
    omega[2] = m_omega;
    std::cout << "set motion Motion000002 to" << m_omega << std::endl;
}

#endif // MOTION000002_H
