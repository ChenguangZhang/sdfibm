#include "vector.H"

/* This namespace collects sdf-related functions (mostly inlined for efficency). */
namespace sdf
{
    using Foam::scalar;
    using Foam::vector;

    const scalar TOL = 1e-8;
    /*********************************************************************************************
     *                                     Shape primitives                                      *
     * Note: all sdf are calculated in the object space                                          *
     ********************************************************************************************/
    // circle
    inline bool circle_bool(const vector& p, const scalar& r)
    {
        return Foam::magSqr(p) < r*r;
    }
    inline bool circle_bool_fast(const vector& p, const scalar& rSQR)
    {
        return Foam::magSqr(p) < rSQR;
    }
    inline scalar circle(const vector& p, const scalar& r)
    {
        return Foam::mag(p) - r;
    }

    // rectangle
    inline scalar rectangle_bool(const vector& p, const scalar& ra, const scalar& rb)
    {
        return std::fabs(p.x()) < ra && std::fabs(p.y()) < rb;
    }
    inline scalar rectangle(const vector& p, const scalar& ra, const scalar& rb)
    {
        scalar dx = std::fabs(p.x()) - ra;
        scalar dy = std::fabs(p.y()) - rb;
        scalar dxp = std::max(0.0, dx);
        scalar dyp = std::max(0.0, dy);
        return std::sqrt(dxp*dxp + dyp*dyp) + std::min(0.0, std::max(dx, dy));
    }

    // box
    inline scalar box_bool(const vector& p, const scalar& ra, const scalar& rb, const scalar& rc)
    {
        return std::fabs(p.x()) < ra && std::fabs(p.y()) < rb && std::fabs(p.z()) < rc;
    }
    inline scalar box(const vector& p, const scalar& ra, const scalar& rb, const scalar& rc)
    {
        scalar dx = std::fabs(p.x()) - ra;
        scalar dy = std::fabs(p.y()) - rb;
        scalar dz = std::fabs(p.z()) - rc;
        scalar dxp = std::max(0.0, dx);
        scalar dyp = std::max(0.0, dy);
        scalar dzp = std::max(0.0, dz);
        return std::sqrt(dxp*dxp+dyp*dyp+dzp*dzp)+std::min(0.0,std::max(dz,std::max(dx, dy)));
    }

    // ellipse
    inline scalar ellipse_bool_fast(const vector& p, const scalar& raSQRINV, const scalar& rbSQRINV)
    {
        return p.x()*p.x()*raSQRINV + p.y()*p.y()*rbSQRINV < 1.0;
    }
    inline scalar ellipse(const vector& p, const scalar& raSQRINV, const scalar& rbSQRINV)
    {
        // raSQRINV = 1/(a*a) and so on
        scalar xbyaSQR = p.x()*p.x()*raSQRINV; // (x/a)^2
        scalar ybybSQR = p.y()*p.y()*rbSQRINV; // (y/b)^2
        return 0.5*(xbyaSQR+ybybSQR-1.0)/(std::sqrt(xbyaSQR*raSQRINV+ybybSQR*rbSQRINV));
    }

    // ellipsoid
    inline scalar ellipsoid_bool_fast(const vector& p, const scalar& raSQRINV, const scalar& rbSQRINV, const scalar& rcSQRINV)
    {
        return p.x()*p.x()*raSQRINV + p.y()*p.y()*rbSQRINV + p.z()*p.z()*rcSQRINV < 1.0;
    }
    inline scalar ellipsoid(const vector& p, const scalar& raSQRINV, const scalar& rbSQRINV, const scalar& rcSQRINV)
    {
        // raSQRINV = 1/(a*a), and so on
        scalar xbyaSQR = p.x()*p.x()*raSQRINV; // (x/a)^2
        scalar ybybSQR = p.y()*p.y()*rbSQRINV; // (y/b)^2
        scalar zbycSQR = p.z()*p.z()*rcSQRINV; // (z/c)^2
        return 0.5*(xbyaSQR+ybybSQR+zbycSQR-1.0)/(std::sqrt(xbyaSQR*raSQRINV+ybybSQR*rbSQRINV+zbycSQR*rcSQRINV));
    }

    /*********************************************************************************************
     *                                     Transformations                                       *
     * Note: to transform the field, the point transformation is reversed.                       *
     * For example, to rotate the field +30\degree, hte point is rotated -30\degree.             *
     ********************************************************************************************/
    // rotation
    inline vector rot30(const vector& p)
    {
        return vector(0.866025404*p.x()+0.5*p.y(), 0.866025404*p.y()-0.5*p.x(), 0.0);
    }
    inline vector rot45(const vector& p)
    {
        return 0.707106781*vector(p.x()+p.y(),-p.x()+p.y(), 0.0);
    }
    inline vector rot60(const vector& p)
    {
        return vector(0.866025404*p.y()+0.5*p.x(),-0.866025404*p.x()+0.5*p.y(), 0.0);
    }
    inline vector rot90(const vector& p)
    {
        return vector(p.y(), -p.x(), 0.0);
    }
    inline vector rotth(const vector& p, const scalar& th)
    {
        scalar s = std::sin(th);
        scalar c = std::cos(th);
        return vector(p.x()*c+p.y()*s,-p.x()*s+p.y()*c, 0.0);
    }

    // reflection
    inline vector flipy(const vector& p)
    {
        return vector( p.x(), -p.y(), p.z());
    }
    inline vector flipx(const vector& p)
    {
        return vector(-p.x(),  p.y(), p.z());
    }
    inline vector offset(const vector& p, const vector& offset)
    {
        return p - offset;
    }

    /*********************************************************************************************
     *                                    Boolean Operations                                     *
     * Note: the operation on Boolean and floating point values are different                    *
     ********************************************************************************************/
    // difference (a binary operation)
    inline scalar D(const scalar& d1, const scalar& d2) { return std::max(d1,-d2);}
    inline bool   D(bool d1,          bool d2         ) { return      d1 && (!d2);}

    // union
    inline scalar U(const std::initializer_list<scalar>& phis) {return std::min(phis);}
    inline bool   U(const std::initializer_list<bool>  & phis) {return std::max(phis);}

    // intersection
    inline scalar I(const std::initializer_list<scalar>& phis) {return std::max(phis);}
    inline bool   I(const std::initializer_list<bool>  & phis) {return std::min(phis);}

    /*********************************************************************************************
     *                                     Helper functions                                      *
     ********************************************************************************************/
    inline scalar filter(const scalar& phi)
    {
        return (std::fabs(phi)<TOL) ? -TOL : phi;
    }

    inline scalar clamp(scalar x, scalar a, scalar b)
    {
        return (x>b) ? b: ((x<a) ? a : x);
    }
}

