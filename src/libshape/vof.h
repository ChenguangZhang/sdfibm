#ifndef VOF_H
#define VOF_H

#include "ishape.h"
#include "fvc.H"

real cellFraction2D(const IShape* shape,
                    const vector& shape_center,
                    const quaternion& shape_orientation,
                    const Foam::scalarField& cs,
                    const Foam::vectorField& cc,
                    label celli);
real cellFraction3D(const IShape* shape,
                    const vector& shape_center,
                    const quaternion& shape_orientation,
                    const Foam::scalarField& cs,
                    const Foam::vectorField& cc,
                    label celli);
#endif
