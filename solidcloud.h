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

class SolidCloud
{
private:
    // control variables (note that solid is always on)
    bool m_ON_FLUID; // fluid is on by default, but could be disabled
    bool m_ON_TWOD;  // two dimensional
    bool m_ON_RESTART;

private:
    std::vector<Solid> m_solids; // finite solids
    std::vector<Solid> m_planes; // infinite solids

    // bind the correct 2D/3D function to process intersected cells
    std::function<real (const IShape*,
                        const vector&,
                        const quaternion&,
                        const Foam::scalarField&,
                        const Foam::vectorField&,
                        label)> m_vofCalculator;

    // collision handling
    std::vector<CollisionPair> collision_pairs;
    UGrid* m_ptr_ugrid;
    BBox* m_ptr_bbox;

    // environments
    vector m_gravity;
    real   m_rho; // fluid density

    real m_radiusB;

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
    void solidFluidInteract(Solid& s,  const real& dt);

    //void resolveCollisionPairs();

    std::map<std::string, IMotion*  > m_libmotion;
    std::map<std::string, IMaterial*> m_libmat;
    std::map<std::string, IShape*   > m_libshape;
    std::ofstream statefile;
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
    //void addBoundingBox(const BoundingBox& particle_bounding_box);
    void linkFluid(const Foam::volVectorField& U);
    void linkMesh(const Foam::fvMesh& mesh);
    void initialCorrect();

    // io
    void saveState(const real& time);
    void loadRestart(const std::string& filename);
    void saveRestart(const std::string& filename);

    const Solid& operator[](label i) const;
    void checkAlpha() const;
    real totalSolidVolume() const;
    bool inline isOnFluid() const {return m_ON_FLUID;}
    bool inline isOnTwoD() const  {return m_ON_TWOD;}

    // related to time stepping
    void storeOld();
    void restoreOld();
    void evolve  (const real& time, const real& dt); // returns the max change in cloud
    void interact(const real& time, const real& dt);
    void addMidEnvironment();
    void fixInternal();
};

#endif
