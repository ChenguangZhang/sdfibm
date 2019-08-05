#ifndef MOTIONSINEDIRECTIONAL_H
#define MOTIONSINEDIRECTIONAL_H

#include "imotion.h"

class MotionSineDirectional:public IMotion, _creator<MotionSineDirectional>
{
public:
    // same signature for all motions
    virtual void constraint(
            const real& time,
            vector& velocity,
            vector& omega) override final;

    // update below
    MotionSineDirectional(const dictionary& para)
    {
        try {
            m_amplitude = Foam::readScalar(para.lookup("amplitude"));
            m_period    = Foam::readScalar(para.lookup("period"));
            m_omega = 2*M_PI/m_period;

            vector direction = para.lookup("direction");
            m_direction = direction;

            if(Foam::magSqr(m_direction) > 1.001)
                throw("direction vector not normalized!");
        } catch (const std::exception& e) {
            std::cout << "Problem in creating MotionSineDirectional (aka linear oscillation)"
                      << e.what() << std::endl;
        }
    }
    virtual ~MotionSineDirectional() override final {}
    TYPENAME("MotionSineDirectional") // a linear motion
    virtual string description() const override {return "oscillate as vec(x) = vec(x0) + a*sin(omega t)*vec(n), no rotation";}
private:
    real m_amplitude, m_period;
    vector m_direction;

    real m_omega;
};

void MotionSineDirectional::constraint(const real &time, vector &velocity, vector &omega)
{
    velocity = m_amplitude * m_omega * std::cos(m_omega* time)*m_direction;
    omega    = vector::zero;
}


#endif // MOTIONSINEDIRECTIONAL_H
