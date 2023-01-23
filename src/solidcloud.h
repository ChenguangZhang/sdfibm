#ifndef SOLIDCLOUD_H
#define SOLIDCLOUD_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "dimensionedScalar.H"
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
    // disallow copy constructor and assignment operator
    SolidCloud (const SolidCloud&) = delete;
    SolidCloud& operator=(const SolidCloud&) = delete;

private:
    bool m_ON_FLUID; // is fsi?
    bool m_ON_TWOD;  // is 2-d?
    bool m_ON_VOFONLY;
    bool m_ON_RESTART;

    unsigned int m_timeStepCounter;
    unsigned int m_writeFrequency;
    scalar m_time;

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

    void solidFluidInteract(Solid& s, scalar dt);
    void solidFluidCorrect (Solid& s, scalar dt);

public:
    SolidCloud(const Foam::word& dictfile, Foam::volVectorField& U, scalar time);
    ~SolidCloud();

    // setup
    inline void addSolid(Solid&& solid) { m_solids.emplace_back(solid); }
    inline void addPlane(Solid&& solid) { m_planes.emplace_back(solid); }
    void addBoundingBox(const BBox& particle_bbox);

    // io
    void saveState();
    void initFromDictionary(const Foam::word& dictname);
    void saveRestart(const std::string& filename);
    const Solid& operator[](label i) const;

    // info
    void checkAlpha() const;
    scalar totalSolidVolume() const;
    bool inline isOnFluid() const {return m_ON_FLUID;}
    bool inline isOnTwoD()  const {return m_ON_TWOD;}

    // time stepping
    void storeOld();
    void restoreOld();
    void evolve  (scalar time, scalar dt);
    void interact(scalar time, scalar dt);
    void addMidEnvironment();
    void fixInternal(scalar dt);
    void initialCorrect();

    friend std::ostream& operator<<(std::ostream& os, const SolidCloud& sc);
};

}
#endif
