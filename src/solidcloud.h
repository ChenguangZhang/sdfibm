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
#include "logger.h"


#include "./libcollision/ugrid.h"
#include "./libcollision/bbox.h"
#include "./libcollision/collision.h"
namespace sdfibm{

class SolidCloud
{
private:
    // control variables (note that solid is always on)
    bool m_ON_FLUID; // fluid is on by default, but could be disabled
    bool m_ON_TWOD;  // two dimensional
    bool m_ON_VOFONLY;
    bool m_ON_RESTART;

private:
    std::vector<Solid> m_solids; // finite solids
    std::vector<Solid> m_planes; // infinite planes (also solids)
    dictionary m_solidDict;

    // function binding of 2D or 3D geometric function to handle intersected cells
    /* std::function<real (const IShape*, const vector&, const quaternion&, */
    /*                     const Foam::scalarField&, const Foam::vectorField&, label)> m_vofCalculator; */

    // collision handling
    BBox* m_ptr_bbox;
    UGrid* m_ptr_ugrid;
    std::vector<CollisionPair> collision_pairs;
    real m_radiusB; // maximum bounding radius of shapes

    // environmental variables
    vector m_gravity; // gravity
    real   m_rho;     // fluid density

    // fluid
    Foam::fvMesh*         m_ptr_Mesh;
    Foam::volVectorField* m_ptr_Uf;

    Foam::scalarField* m_ptr_cs; // cell size
    Foam::volScalarField* m_ptr_ct; // cell type
    // solid's projected fields on fluid mesh
    Foam::volScalarField* m_ptr_As;
    Foam::volVectorField* m_ptr_Fs;
    Foam::volScalarField* m_ptr_Ts;

    Foam::surfaceScalarField* m_ptr_Asf;
    Foam::meshSearch* m_ms;

private:
    void solidSolidInteract();
    void solidSolidCollision(Solid& s1, Solid& s2);
    void solidFluidInteract(Solid& s,  const real& dt);

    void resolveCollisionPairs();

    std::map<std::string, IMotion*  > m_libmotion;
    std::map<std::string, IMaterial*> m_libmat;
    std::map<std::string, IShape*   > m_libshape;
    std::ofstream statefile;

    inline const std::string GenBanner(const std::string& title) const
    {
        unsigned int nside = (78 - title.length())/2;
        return std::string(nside, '*') + ' ' + title +  ' ' + std::string(nside, '*') + '\n';
    }

public:
    SolidCloud(const Foam::word& dictfile);
    ~SolidCloud();

    // related to initial setup
    inline void addSolid(const Solid& solid) { m_solids.push_back(solid); }
    inline void addPlane(const Solid& solid) { m_planes.push_back(solid); }
    void addBoundingBox(const BBox& particle_bbox);
    void linkFluid(const Foam::volVectorField& U);
    void initialCorrect();

    // io
    void saveState(const real& time);
    void initFromDictionary(const Foam::word& dictname);// arg is dictionary name
    void saveRestart(const std::string& filename);

    const Solid& operator[](label i) const;
    void checkAlpha() const;
    real totalSolidVolume() const;
    bool inline isOnFluid() const {return m_ON_FLUID;}
    bool inline isOnTwoD() const  {return m_ON_TWOD;}

    // related to time stepping
    void storeOld();
    void restoreOld();
    void evolve  (const real& time, const real& dt);
    void interact(const real& time, const real& dt);
    void addMidEnvironment();
    void fixInternal();
};

}
#endif

