#include "vof.h"
#include "kernel.h"


// functions that map a geometry into a volume fraction field on a Cartesian mesh mesh
real cellFraction2D(const IShape* shape,
                    const vector& shape_center,
                    const quaternion& shape_orientation,
                    const Foam::scalarField& cs,
                    const Foam::vectorField& cc,
                    label celli)
{
    real cellsize = 0.5*cs[celli];
    vector c = cc[celli];
    c[2] = 0;

    vector shift1(cellsize, cellsize, 0.0);
    vector shift2(cellsize,-cellsize, 0.0);
    real ur = shape->signedDistance(c + shift1, shape_center, shape_orientation);
    real ul = shape->signedDistance(c - shift2, shape_center, shape_orientation);
    real lr = shape->signedDistance(c + shift2, shape_center, shape_orientation);
    real ll = shape->signedDistance(c - shift1, shape_center, shape_orientation);
    return squareFraction(ll, lr, ur, ul);
}

real cellFraction3D(const IShape* shape,
                    const vector& shape_center,
                    const quaternion& shape_orientation,
                    const Foam::scalarField& cs,
                    const Foam::vectorField& cc,
                    label celli)
{
        real cellsize = cs[celli];
        vector c = cc[celli];

        vector origin = c - 0.5*vector(cellsize, cellsize, cellsize);
        vector unitx(cellsize, 0, 0);
        vector unity(0, cellsize, 0);
        vector unitz(0, 0, cellsize);

        real d[8]; // SDF at 8 vertices
        d[0] = shape->signedDistance(origin,                 shape_center, shape_orientation);
        d[1] = shape->signedDistance(origin + unitx,         shape_center, shape_orientation);
        d[2] = shape->signedDistance(origin + unitx + unity, shape_center, shape_orientation);
        d[3] = shape->signedDistance(origin + unity,         shape_center, shape_orientation);
        d[4] = shape->signedDistance(origin + unitz,         shape_center, shape_orientation);
        d[5] = shape->signedDistance(origin + unitz + unitx, shape_center, shape_orientation);
        d[6] = shape->signedDistance(origin + unitz + unitx + unity, shape_center, shape_orientation);
        d[7] = shape->signedDistance(origin + unitz + unity, shape_center, shape_orientation);
        return cubeFraction(d);
}
