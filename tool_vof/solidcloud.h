#ifndef SOLIDCLOUD_H
#define SOLIDCLOUD_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <functional>

#include "fvm.H"
#include "types.h"
#include "solid.h"
#include "./libmaterial/imaterial.h"
#include "utils.h"
#include "meshSearch.H"

#include "./libcollision/ugrid.h"
#include "./libcollision/bbox.h"
#include "./libcollision/collision.h"
namespace sdfibm{

class SolidCloud
{
private:
    // control variables (note that solid is always on)
    bool m_ON_TWOD;  // two dimensional
    std::string m_field_name;

private:
    std::vector<Solid> m_solids; // finite solids
    std::vector<Solid> m_planes; // infinite planes (also solids)

    // collision handling
    real m_radiusB; // maximum bounding radius of shapes

    // fluid
    Foam::fvMesh*         m_ptr_Mesh;

    Foam::scalarField* m_ptr_cs; // cell size
    Foam::volScalarField* m_ptr_ct; // cell type
    // solid's projected fields on fluid mesh
    Foam::volScalarField* m_ptr_As;
    Foam::meshSearch* m_ms;

private:
    void solidFluidInteract(Solid& s,  const real& dt);

    void resolveCollisionPairs();

    std::map<std::string, IMotion*  > m_libmotion;
    std::map<std::string, IMaterial*> m_libmat;
    std::map<std::string, IShape*   > m_libshape;
    std::ofstream logfile;

    // helper functions to print info to the console
    inline const std::string GenBanner(const std::string& title) const
    {
        unsigned int nside = (78 - title.length())/2;
        return std::string(nside, '*') + ' ' + title +  ' ' + std::string(nside, '*') + '\n';
    }

public:
    SolidCloud();
    ~SolidCloud();

    // related to initial setup
    inline void addSolid(const Solid& solid) { m_solids.push_back(solid); }
    inline void addPlane(const Solid& solid) { m_planes.push_back(solid); }
    void addBoundingBox(const BBox& particle_bbox);
    void linkFluid(const Foam::volScalarField& U);
    void initialCorrect();

    void setFieldName(const std::string& field_name){m_field_name = field_name;}

    const Solid& operator[](label i) const;
    void checkAlpha() const;
    real totalSolidVolume() const;
    bool inline isOnTwoD() const  {return m_ON_TWOD;}

    // related to time stepping
    void interact(const real& time, const real& dt);
};

}
#endif

