#ifndef ISHAPE_H
#define ISHAPE_H

#include "../types.h"
#include <algorithm>

#define SHAPETYPENAME(name) \
    static string typeName() {return name;} \
    static bool added;

class IShape;
template <typename T>
class _shapecreator
{
public:
    static IShape* create(const dictionary& para)
    {
        return new T(para);
    }
};

class IShape
{
public:
    const static int m_id = -1;
    real   m_radiusB; // radius of bounding sphere
    real   m_volume, m_volumeINV; // volume
    vector m_com; // center of mass
    // moi always in principal frame (nonzero only along diagonal)
    tensor m_moi, m_moiINV;

public:
    IShape()
    {
        // set default values
        m_radiusB   = 0.0; // bounding raidus
        m_volume    = 0.0; // volume
        m_volumeINV = 0.0; // 1/volume
        m_com = vector(0, 0, 0); // center of mass

        m_moi    = tensor(0, 0, 0, 0, 0, 0, 0, 0, 0);
        m_moiINV = tensor(0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    virtual int  getShapeID() const {return m_id;}
    virtual real getRadiusB() const {return m_radiusB;}

    // collection of inline sdf utilities, all are pure geometric functions
    // transformtion functions
    // because we are transforming the field, the points transformation are all "reversed"
    // for example rotate the field 30\degree counter-clockwise is equivalent to rotate the point 30\degree clockwise.
    static inline vector _sdf_rot30(const vector& p) {
        return vector( 0.866025404*p.x()+0.5*p.y(), 0.866025404*p.y()-0.5*p.x(),0.0);
    }
    static inline vector _sdf_rot45(const vector& p) {
        return 0.707106781*vector(p.x()+p.y(),-p.x()+p.y(), 0.0);
    }
    static inline vector _sdf_rot60(const vector& p) {
        return vector(0.866025404*p.y()+0.5*p.x(),-0.866025404*p.x()+0.5*p.y(),0.0);
    }
    static inline vector _sdf_rot90(const vector& p) {return vector(p.y(), -p.x(), 0.0);}

    static inline vector _sdf_rotth(const vector& p, const real& th) {
        real s = std::sin(th);
        real c = std::cos(th);
        return vector(p.x()*c+p.y()*s,-p.x()*s+p.y()*c, 0.0);
    }

    // mirror regarding to x-axis and y-axis
    static inline vector _sdf_flipy(const vector& p) {return vector( p.x(), -p.y(), p.z());}
    static inline vector _sdf_flipx(const vector& p) {return vector(-p.x(),  p.y(), p.z());}
    static inline vector _sdf_offset(const vector& p, const vector& offset) {return p - offset;}

    // boolean operations
    // diff(erence) is binary function
    static inline real _sdf_diff(const real& d1, const real& d2) { return std::max(d1,-d2);}
    static inline bool _sdf_diff(bool d1, bool d2) { return d1 && (!d2);}
    template<typename T> static inline T _sdf_union(const std::initializer_list<T>& sdfs) {return std::min(sdfs);}
    template<typename T> static inline T _sdf_inter(const std::initializer_list<T>& sdfs) {return std::max(sdfs);}
    static inline bool _sdf_union(const std::initializer_list<bool>& sdfs) {return std::max(sdfs);}
    static inline bool _sdf_inter(const std::initializer_list<bool>& sdfs) {return std::min(sdfs);}

    // safeguard function
    static inline real _sdf_filter(const real& sdf) {return (std::fabs(sdf)<1.0e-8) ? -1e-8 : sdf;}

    // all w.r.t the origin
    static inline bool _sdf_circle_bool(const vector& p, const real& r) {
        return Foam::magSqr(p) < r*r;
    }
    static inline bool _sdf_circle_bool_fast(const vector& p, const real& rSQR) {
        return Foam::magSqr(p) < rSQR;
    }
    static inline real _sdf_circle_real(const vector& p, const real& r) {
        return Foam::mag(p) - r;
    }
    static inline real _sdf_rectangle_bool(const vector& p, const real& ra, const real& rb) {
        return std::fabs(p.x()) < ra && std::fabs(p.y()) < rb;
    }
    static inline real _sdf_rectangle_real(const vector& p, const real& ra, const real& rb) {
        real dx = std::fabs(p.x()) - ra;
        real dy = std::fabs(p.y()) - rb;
        real dxp = std::max(0.0, dx);
        real dyp = std::max(0.0, dy);
        return std::sqrt(dxp*dxp + dyp*dyp) + std::min(0.0, std::max(dx, dy));
    }
    static inline real _sdf_box_bool(const vector& p, const real& ra, const real& rb, const real& rc) {
        return std::fabs(p.x()) < ra && std::fabs(p.y()) < rb && std::fabs(p.z()) < rc;
    }
    static inline real _sdf_box_real(const vector& p, const real& ra, const real& rb, const real& rc) {
        real dx = std::fabs(p.x()) - ra;
        real dy = std::fabs(p.y()) - rb;
        real dz = std::fabs(p.z()) - rc;
        real dxp = std::max(0.0, dx);
        real dyp = std::max(0.0, dy);
        real dzp = std::max(0.0, dz);
        return std::sqrt(dxp*dxp+dyp*dyp+dzp*dzp)+std::min(0.0,std::max(dz,std::max(dx, dy)));
    }
    static inline real _sdf_ellipse_bool_fast(const vector& p, const real& raSQRINV, const real& rbSQRINV) {
        return p.x()*p.x()*raSQRINV + p.y()*p.y()*rbSQRINV < 1.0;
    }
    static inline real _sdf_ellipse_real_fast(const vector& p, const real& raSQRINV, const real& rbSQRINV) {
        // raSQRINV = 1/(a*a), and rbSQRINV = 1/(b*b)
        real xbyaSQR = p.x()*p.x()*raSQRINV; // (x/a)^2
        real ybybSQR = p.y()*p.y()*rbSQRINV; // (y/b)^2
        return 0.5*(xbyaSQR+ybybSQR-1.0)/(std::sqrt(xbyaSQR*raSQRINV+ybybSQR*rbSQRINV));
    }

    static inline real _sdf_ellipsoid_bool_fast(const vector& p, const real& raSQRINV, const real& rbSQRINV, const real& rcSQRINV) {
        return p.x()*p.x()*raSQRINV + p.y()*p.y()*rbSQRINV + p.z()*p.z()*rcSQRINV < 1.0;
    }
    static inline real _sdf_ellipsoid_real_fast(const vector& p,
            const real& raSQRINV,
            const real& rbSQRINV,
            const real& rcSQRINV) {
        // raSQRINV = 1/(a*a), and rbSQRINV = 1/(b*b)
        real xbyaSQR = p.x()*p.x()*raSQRINV; // (x/a)^2
        real ybybSQR = p.y()*p.y()*rbSQRINV; // (y/b)^2
        real zbycSQR = p.z()*p.z()*rcSQRINV; // (z/c)^2
        return 0.5*(xbyaSQR+ybybSQR+zbycSQR-1.0)/(std::sqrt(xbyaSQR*raSQRINV+ybybSQR*rbSQRINV+zbycSQR*rcSQRINV));
    }

    // transfer from world to body frame, note the use of conjugate
    static vector transform(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_quaternion)
    {
        return Foam::conjugate(shape_quaternion).transform(p-shape_center);
    }
    // above content are commont to all shapes -- DO NOT MODIFY
    // below content to be overwritten by child classes
    virtual bool isInside(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const = 0;

    virtual real signedDistance(
            const vector& p,
            const vector& shape_center,
            const quaternion& shape_orientation) const = 0;
    virtual std::string description() const = 0; // more detailed information
    virtual ~IShape(){}
};

#endif
