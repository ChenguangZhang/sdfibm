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
#include "logger.h"
#include "geometrictools.h"
#include "cellenumerator.h"

#include "./libcollision/ugrid.h"
#include "./libcollision/bbox.h"
#include "./libcollision/collision.h"

namespace sdfibm {

class SolidCloud
{
private:
    bool m_ON_FLUID; // is fsi?
    bool m_ON_TWOD;  // is 2-d?
    bool m_ON_VOFONLY;
    bool m_ON_RESTART;
    SolidCloud (const SolidCloud&) = delete;
    SolidCloud& operator=(const SolidCloud&) = delete;

private:
    std::vector<Solid> m_solids; // finite solids
    std::vector<Solid> m_planes; // infinite planes (also solids)
    dictionary m_solidDict;

    // collision handling
    BBox* m_ptr_bbox;
    UGrid* m_ptr_ugrid;
    std::vector<CollisionPair> collision_pairs;
    scalar m_radiusB; // max bounding radius of shapes

    // environmental info
    vector m_gravity;
    scalar m_rhof; // fluid density

    // variables on mesh
    const Foam::fvMesh&   m_mesh;
    Foam::volVectorField& m_Uf;
    Foam::volScalarField& m_ct;
    Foam::volScalarField& m_As;
    Foam::volVectorField& m_Fs;
    Foam::volScalarField& m_Ts;

    GeometricTools m_geotools;
    CellEnumerator m_cellenum;

    std::map<std::string, IMotion*  > m_libmotion;
    std::map<std::string, IMaterial*> m_libmat;
    std::map<std::string, IShape*   > m_libshape;
    std::ofstream statefile;

private:
    void solidSolidInteract();
    void solidSolidCollision(Solid& s1, Solid& s2);
    void resolveCollisionPairs();

    void solidFluidInteract(Solid& s,  const scalar& dt);
    void solidFluidCorrect(Solid& s, const scalar& dt);

public:
    SolidCloud(const Foam::word& dictfile, Foam::volVectorField& U);
    ~SolidCloud();

    // setup
    inline void addSolid(const Solid& solid) { m_solids.push_back(solid); }
    inline void addPlane(const Solid& solid) { m_planes.push_back(solid); }
    void addBoundingBox(const BBox& particle_bbox);

    // io
    void saveState(const scalar& time);
    void initFromDictionary(const Foam::word& dictname);
    void saveRestart(const std::string& filename);
    const Solid& operator[](label i) const;

    // info
    void checkAlpha() const;
    scalar totalSolidVolume() const;
    bool inline isOnFluid() const {return m_ON_FLUID;}
    bool inline isOnTwoD() const  {return m_ON_TWOD;}

    // time stepping
    void storeOld();
    void restoreOld();
    void evolve  (const scalar& time, const scalar& dt);
    void interact(const scalar& time, const scalar& dt);
    void addMidEnvironment();
    void fixInternal(const scalar& dt);
    void initialCorrect();
};

}
#endif

