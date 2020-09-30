#ifndef SOLIDCLOUD_H
#define SOLIDCLOUD_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "fvm.H"
#include "types.h"
#include "solid.h"
#include "./libmaterial/imaterial.h"
#include "utils.h"
#include "geometrictools.h"
#include "cellenumerator.h"

namespace sdfibm {

class SolidCloud
{
private:
    bool m_ON_TWOD;

private:
    std::vector<Solid> m_solids;
    std::vector<Solid> m_planes;

    const Foam::fvMesh& m_mesh;
    Foam::volScalarField* m_ptr_As;

private:
    void solidFluidInteract(Solid& s,  const scalar& dt);

    std::map<std::string, IMotion*  > m_libmotion;
    std::map<std::string, IMaterial*> m_libmat;
    std::map<std::string, IShape*   > m_libshape;
    std::ofstream logfile;

    GeometricTools m_geotools;
    CellEnumerator m_cellenum;

public:
    SolidCloud(const Foam::word& dictfile, const Foam::fvMesh& mesh);
    ~SolidCloud(){};

    inline void addSolid(const Solid& solid) { m_solids.push_back(solid); }
    inline void addPlane(const Solid& solid) { m_planes.push_back(solid); }
    void writeVOF(const Foam::word& field_name);
};

}
#endif

