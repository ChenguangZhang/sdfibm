#ifndef FORCESPRING_H
#define FORCESPRING_H

#include "iforcer.h"

namespace sdfibm::force {

class Spring : public IForcer, _creator<Spring>
{
public:
    virtual Force generate(const scalar &time, const vector& position, const vector& velocity) override final;
    // update below
    Spring(const dictionary& para)
    {
        pivot = para.lookup("pivot");
        try {
            k = Foam::readScalar(para.lookup("k"));
            l = Foam::readScalar(para.lookup("l"));
        } catch (const std::exception& e) {
            std::cout << "Problem in creating Spring force" << ' '
                      << e.what() << std::endl;
            std::exit(1);
        }
    }
    virtual ~Spring() override final {}
    TYPENAME("Spring") // a linear motion
    virtual std::string description() const override {return "Spring forcer with a pivot, stiffness (k), and rest length (l)";}
private:
    vector pivot; 
    scalar k;
    scalar l;
};

IForcer::Force Spring::generate(const scalar &time, const vector &position, const vector &velocity)
{
    vector r = position - pivot;
    return {-k * r * (1.0 - l/mag(r)), vector::zero};
}

}
#endif
