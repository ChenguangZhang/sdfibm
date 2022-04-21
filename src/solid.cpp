#include "solid.h"

namespace sdfibm {

std::ostream& operator<<(std::ostream& os, const Solid& s)
{
    vector v;
    v = s.getCenter();  os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    v = s.getVelocity();os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    v = s.getForce();   os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';

    v = s.getOrientation().eulerAngles(quaternion::XYZ);
                            os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    v = s.getOmega();   os << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    v = s.getTorque();  os << v.x() << ' ' << v.y() << ' ' << v.z();
    return os;
}

void write2D(std::ostream& os, const Solid& s)
{
    vector v;
    v = s.getCenter();  os << v.x() << ' ' << v.y() << ' ';
    v = s.getVelocity();os << v.x() << ' ' << v.y() << ' ';
    v = s.getForce();   os << v.x() << ' ' << v.y() << ' ';
    v = s.getOrientation().eulerAngles(quaternion::XYZ);
                            os << v.z() << ' ';
    v = s.getOmega();   os << v.z() << ' ';
    v = s.getTorque();  os << v.z();
}

}