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

#include "utils.h"
#include "logger.h"
#include "geometrictools.h"
#include "cellenumerator.h"

#include "entitylibrary.h"

#include "./libcollision/ugrid.h"
#include "./libcollision/bbox.h"
#include "./libcollision/collision.h"

class IMotion;
class IMaterial;
class IShape;

namespace sdfibm::force {
    class IForcer;
}

namespace sdfibm {

class SolidCloud
{
private:
    bool m_ON_FLUID;
    bool m_ON_TWOD;
    bool m_ON_VOFONLY;
    bool m_ON_RESTART;
    bool m_ON_MEANFIELD {false};

    unsigned int m_timeStepCounter;
    unsigned int m_writeFrequency;
    scalar m_time;
    std::string m_sampler;

    std::vector<Solid> m_solids;
    dictionary m_solidDict;

    // collision handling
    BBox* m_ptr_bbox;
    UGrid* m_ptr_ugrid;
    std::vector<CollisionPair> collision_pairs;
    scalar m_radiusB;

    // environmental info
    vector m_gravity;
    scalar m_rhof;

    // variables on mesh
    const Foam::fvMesh&   m_mesh;
    Foam::volVectorField& m_Uf;
    Foam::volScalarField& m_ct;
    Foam::volScalarField& m_As;
    Foam::volVectorField& m_Fs;
    Foam::volScalarField& m_Ts;

    GeometricTools m_geotools;
    std::unique_ptr<Foam::meshSearch> m_ms;

    std::map<std::string, IMotion*  > m_libmotion;
    std::map<std::string, IMaterial*> m_libmat;
    EntityLibrary<IShape> m_libshape;
    EntityLibrary<force::IForcer> m_libforcer;
    std::ofstream statefile;
    std::ofstream meanFieldFile;

    void solidSolidInteract();
    void solidSolidCollision(Solid& s1, Solid& s2);
    void resolveCollisionPairs();

    void solidFluidInteract(Solid& s, scalar dt);
    void solidFluidCorrect (Solid& s, scalar dt);
    template<class Type>
    Type calcMeanField(Solid& s, IShape* shape, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& field);

public:
    SolidCloud(const Foam::word& dictfile, Foam::volVectorField& U, scalar time);
    ~SolidCloud();

    // setup
    inline void addSolid(Solid&& solid) { m_solids.emplace_back(solid); }
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
    void writeMeanField();
    void evolve  (scalar time, scalar dt);
    void interact(scalar time, scalar dt);
    void addMidEnvironment();
    void fixInternal(scalar dt);
    void initialCorrect();

    friend std::ostream& operator<<(std::ostream& os, const SolidCloud& sc);
    
    // deleted methods
    SolidCloud (const SolidCloud&) = delete;
    SolidCloud& operator=(const SolidCloud&) = delete;
};

}
#endif
