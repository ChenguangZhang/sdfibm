#ifndef MOTIONROTOR_H
#define MOTIONROTOR_H

#include "imotion.h"
namespace sdfibm{

class MotionRotor:public IMotion, _creator<MotionRotor>
{
public:
    // same signature for all motions
    virtual void constraint(
            const scalar& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    MotionRotor(const dictionary& para)
    {
        try {
            m_period = Foam::readScalar(para.lookup("period"));
            m_radius = Foam::readScalar(para.lookup("radius"));
            m_theta0 = Foam::readScalar(para.lookup("theta0"));
            m_selfom = Foam::readScalar(para.lookup("selfom"));
        } catch (const std::exception& e) {
            std::cout << "Problem in creating MotionRotor" << ' '
                      << e.what() << std::endl;
            std::exit(1);
        }
        m_omega = 2*M_PI/m_period;
    }
    virtual ~MotionRotor() override final {}
    TYPENAME("MotionRotor")
    virtual string description() const override {return "rotate around origin with prescribed omega";}
private:
    scalar m_period, m_omega;
    scalar m_radius, m_theta0;
    scalar m_selfom;
};

void MotionRotor::constraint(const scalar &time, vector &velocity, vector &omega)
{
    velocity = vector::zero;
    velocity[0] = -m_radius*m_omega*std::sin(m_omega*time + m_theta0);
    velocity[1] =  m_radius*m_omega*std::cos(m_omega*time + m_theta0);
    omega    = vector::zero;
    omega[2] = m_selfom;
}

}
#endif // MOTIONROTOR_H
