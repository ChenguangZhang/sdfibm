#pragma once

#include "iforcer.h"

namespace sdfibm::forcer {

class Magnetic : public IForcer, _creator<Magnetic>
{
public:
    virtual Force generate(const scalar&, const vector&, const vector&, const quaternion&, const vector&) override final;

    Magnetic(const dictionary& para)
    {
        direction = para.lookup("direction");
        // normalize the direction
        try {
            A = Foam::readScalar(para.lookup("A"));
            w = Foam::readScalar(para.lookup("w"));
        } catch (const std::exception& e) {
            std::cout << "Problem in creating Magnetic force" << ' '
                      << e.what() << std::endl;
            std::exit(1);
        }
    }
    virtual ~Magnetic() override final {}
    TYPENAME("Magnetic")
    virtual std::string description() const override {return "Magnetic forcer: A*cos(omega*t)*direction";}
private:
    vector direction;  // direction, should be a unit vector
    scalar A; // A
    scalar w; // omega, set to zero to have a constant force
};

IForcer::Force Magnetic::generate(
        const scalar& time,
        const vector& position,
        const vector& velocity,
        const quaternion& orientation,
        const vector& omega
)
{
    // hardcoded magnetic moment as magnitude * unit direction vector
    static vector m_original = (1.0) * (vector(0.0, 0.0, 1.0));
    vector B = A * std::cos(w*time)*direction; // magnetic field
    vector m = Foam::conjugate(orientation).transform(m_original);
    return {vector::zero, m ^ B};
}

}