#include <cstdlib>
#include <iomanip>
#include <limits>
#include <chrono>
#include "Pstream.H"
#include "solidcloud.h"
#include "sstream"
#include "./libmotion/motionfactory.h"
#include "./libshape/shapefactory.h"
namespace sdfibm {

void SolidCloud::initFromDictionary(const Foam::word& dictfile)
{
    Foam::IFstream ifstream = Foam::IFstream(dictfile);
    dictionary root(ifstream());
    if (Foam::Pstream::master())
    {
        std::ostringstream msg;
        msg << "Init from " << dictfile;
        LOG(msg.str());
    }

    // meta information
    const dictionary& meta = root.subDict("meta");
    m_ON_FLUID = Foam::readBool(meta.lookup("on_fluid"));
    m_ON_TWOD  = Foam::readBool(meta.lookup("on_twod"));
    if (meta.found("on_meanfield"))
    {
        m_ON_MEANFIELD = true;
        meanFieldFile.open("meanfield.out", std::fstream::out);
        meanFieldFile << std::scientific;
    }

    m_gravity  = meta.lookup("gravity");
    m_writeFrequency = meta.lookupOrDefault("writeFrequency", 1);
    meta.found("sampler") ? m_sampler = Foam::word(meta.lookup("sampler")) : "";

    // log meta information
    if (Foam::Pstream::master())
    {
        std::ostringstream msg;
        std::string dim = m_ON_TWOD ? "2D" : "3D";
        std::string type = m_ON_FLUID? "FSI": "DEM (fluid disabled)";

        msg << "Summary: "
            << dim << ' ' 
            << type << ", "
            << "g = " << "("<<m_gravity[0]<<' '<< m_gravity[1]<<' '<< m_gravity[2]<<"). ";
        msg << "Binary built at " << __DATE__ << ' ' << __TIME__ << '\n';
        LOG(msg.str());
    }

   // read and create shape, motion, and material
   try {
        // read and create shapes
        if (Foam::Pstream::master())
        {
            LOGF << GenBanner("CREATE: SHAPES");
            LOGF << "--> Available shape types:\n";
            ShapeFactory::report(LOGF);
            LOGF << "--> Used shapes:\n";
        }
        const dictionary& shapes = root.subDict("shapes");
        m_radiusB = -1.0;
        for (int i=0; i < shapes.size(); ++i)
        {
            const dictionary& para = shapes.subDict(shapes.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));
            std::string name = Foam::word(para.lookup("name"));

            m_libshape[name] = ShapeFactory::create(type, para);
            if (m_libshape[name] == nullptr)
                throw std::runtime_error(std::string("Unrecognized shape type " + type + '\n'));
            m_radiusB = std::max(m_radiusB, m_libshape[name]->getRadiusB());

            if (Foam::Pstream::master())
                LOGF << "[+] " << type << " as " << name << " (" << m_libshape[name]->description() << ")\n";
        }

        // find the maximum bounding radius, used to create the uniform grid for collision detection
        if (Foam::Pstream::master())
        {
            LOGF << "maximum bounding radius is " << m_radiusB << '\n';
        }

        // read and create motions
        if (Foam::Pstream::master())
        {
            LOGF << GenBanner("CREATE: MOTIONS");
            LOGF << "--> Available motion types:\n";
            MotionFactory::report(LOGF);
            LOGF << "--> Used motions:\n";
        }
        const dictionary& motions = root.subDict("motions");
        for (int i=0; i < motions.size(); ++i)
        {
            const dictionary& para = motions.subDict(motions.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));
            std::string name = Foam::word(para.lookup("name"));
            m_libmotion[name] = MotionFactory::create(type, para);
            if (m_libmotion[name] == nullptr)
               throw std::runtime_error(std::string("Unrecognized motion type " + type + '\n'));

            if (Foam::Pstream::master())
                LOGF << "[+] " << type << " as " << name << " (" << m_libmotion[name]->description() << ")\n";
        }

        // read and create materials
        if (Foam::Pstream::master())
            LOGF << GenBanner("CREATE: MATERIALS");
        const dictionary& materials = root.subDict("materials");
        for (int i=0; i < materials.size(); ++i)
        {
            const dictionary& para = materials.subDict(materials.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));

            if (type == "General")
            {
                scalar rho = readScalar(para.lookup("rho"));
                std::string name = Foam::word(para.lookup("name"));
                m_libmat[name] = new IMaterial(rho);
            }
            else
            {
                throw std::runtime_error("Unrecognizable material parameter!");
            }
       }
    }

    catch (const std::exception& e)
    {
        if (Foam::Pstream::master())
        {
            std::cout << e.what();
            LOGF << "Error when creating shape/motion/material!" << e.what() << '\n';
        }
        std::exit(1);
    }

    // create solids
    if (Foam::Pstream::master())
        LOGF << GenBanner("CREATE: SOLIDS");

    const dictionary &solids = root.subDict("solids");
    forAll(solids, i)
    {
        const dictionary &solid = solids.subDict(solids.toc()[i]);

        vector pos = solid.lookup("pos");
        if (m_ON_TWOD)
        {
            if (pos.z() !=  0)
                Quit("Solid must has z=0 in 2D simulation, violated by solid # " + std::to_string(i));
        }
        // create solid
        Solid s(i, pos, quaternion::I);
        // read velocity, euler angle, angular velocity
        vector vel = solid.lookupOrDefault("vel", vector::zero);
        s.setVelocity(vel);
        vector euler = solid.lookupOrDefault("euler", vector::zero);
        s.setOrientation(euler*M_PI/180.0);
        vector omega = solid.lookupOrDefault("omega", vector::zero);
        s.setOmega(omega);

        // add motion, material, and shape
        std::string mot_name = Foam::word(solid.lookup("mot_name"));
        std::string mat_name = Foam::word(solid.lookup("mat_name"));
        std::string shp_name = Foam::word(solid.lookup("shp_name"));
        if (mot_name!="free")
            s.setMotion(m_libmotion[mot_name]);
        s.setShape(m_libshape[shp_name]);
        s.setMaterial(m_libmat[mat_name]);
        this->addSolid(std::move(s));

        if (Foam::Pstream::master())
            LOGF << "Solid " << i << ":" << " motion = " << mot_name
                << ", material = " << mat_name << ", shape = " << shp_name << "\n";
    }
    m_solidDict = root;
}

// ctor
SolidCloud::SolidCloud(const Foam::word& dictfile, Foam::volVectorField& U, scalar time = 0.0) :
    m_mesh(U.mesh()),
    m_Uf (U),
    m_ct (const_cast<Foam::volScalarField&>(m_mesh.lookupObject<Foam::volScalarField>("Ct"))),
    m_As (const_cast<Foam::volScalarField&>(m_mesh.lookupObject<Foam::volScalarField>("As"))),
    m_Fs (const_cast<Foam::volVectorField&>(m_mesh.lookupObject<Foam::volVectorField>("Fs"))),
    m_Ts (const_cast<Foam::volScalarField&>(m_mesh.lookupObject<Foam::volScalarField>("Ts"))),
    m_geotools(GeometricTools(m_mesh)),
    m_cellenum(CellEnumerator(m_mesh))
{
    m_time = time;
    m_timeStepCounter = 0;
    m_writeFrequency = 1;
    m_solids.reserve(10);

    collision_pairs.reserve(100);
    m_ptr_ugrid = nullptr;
    InitCollisionFuncTable();

    initFromDictionary(Foam::word(dictfile));
    if (Foam::Pstream::master())
    {
        LOGF << "Totally [" << m_solids.size() << "] solids.\n";
    }

    // open output file, overwrite if file exists
    statefile.open("cloud.out", std::fstream::out);
    statefile << std::scientific;

    // create bounding box
    Foam::vector low  = m_mesh.bounds().min();
    Foam::vector high = m_mesh.bounds().max();
    scalar a[3] = {low.x(), low.y(), low.z()};
    scalar b[3] = {high.x(),high.y(),high.z()};
    m_ptr_bbox = new BBox(a, b);
    m_ptr_bbox->report(std::cout);
    m_ptr_ugrid = new UGrid(*m_ptr_bbox, 2*m_radiusB);

    const Foam::dictionary& transportProperties = m_mesh.lookupObject<Foam::IOdictionary>("transportProperties");
    const Foam::dimensionedScalar&          rho = transportProperties.lookup("rho");
    m_rhof = rho.value();
    if (!m_ON_FLUID)
        m_rhof = 0.0;

    // correct initial condition to reduce #iteration during 1st time step
    m_ON_RESTART = false;
    if (m_mesh.time().value() > 0)
        m_ON_RESTART = true;
    if (!m_ON_RESTART)
        initialCorrect();

    if (Foam::Pstream::master())
    {
        LOGF << GenBanner("END OF INIT");
        LOGF << '\n';
    }
}

SolidCloud::~SolidCloud()
{
    if (Foam::Pstream::master())
        LOG("Simulation finished! Congratulations!");
    delete m_ptr_ugrid, m_ptr_ugrid = nullptr;
    delete m_ptr_bbox , m_ptr_bbox  = nullptr;
    statefile.close();
}

void SolidCloud::initialCorrect()
{
    this->interact(0, 1);
    m_As.write();
    if (Foam::Pstream::master())
        LOG("Initial As written to 0 directory");
    for (Solid& solid : m_solids)
    {
        solid.clearForceAndTorque();
    }
}

void SolidCloud::fixInternal(scalar dt)
{
    // after updating solid velocity, fix fluid velocity in cells fully covered by solid
    const Foam::vectorField& cc = m_mesh.C().internalField();
    const Foam::volScalarField & ct = m_ct;
    forAll(cc, icell)
    {
        if (ct[icell] >= 4) // if cell totally within a solid
        {
            label id = ct[icell] - 4; // id of the solid that contains this cell
            m_Uf[icell] =  m_solids[id].evalPointVelocity(cc[icell]);
        }
    }
    m_Uf.correctBoundaryConditions();
}

void SolidCloud::writeMeanField()
{
    if (m_sampler == "")
        return;
    for (Solid& solid : m_solids)
    {
        vector vMean = calcMeanField<Foam::vector>(solid, m_libshape[m_sampler], m_Uf);
        meanFieldFile << vMean.x() << ' ' << vMean.y() << ' ' << vMean.z() << ' ';
    }
    meanFieldFile << '\n';
}

template<class Type>
Type SolidCloud::calcMeanField(Solid& solid, IShape* shape, const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& field)
{
    Solid tmpSolid = solid;
    tmpSolid.setShape(shape);

    m_geotools.clearCache();
    const Foam::scalarField& cv = m_mesh.V();

    Type meanField  = vector::zero;
    scalar volume = 0.0;

    m_cellenum.SetSolid(tmpSolid);
    int insideType = tmpSolid.getID() + 4;
    scalar alpha = 0.0;
    while (!m_cellenum.Empty())
    {
        int icur = m_cellenum.GetCurCellInd();
        if (m_cellenum.GetCurCellType() == CellEnumerator::ALL_INSIDE)
        {
            alpha = 1.0;
            m_ct[icur] = insideType;
        }
        else
        {
            m_ct[icur] = m_cellenum.GetCurCellType();
            alpha = m_geotools.calcCellVolume(icur, tmpSolid, m_ON_TWOD)/cv[icur];
        }
        scalar dV = alpha*cv[icur];
        volume    += dV;
        meanField += dV*field[icur];
        m_cellenum.Next();
    }

    if (Foam::Pstream::parRun())
    {
        Foam::reduce(meanField, Foam::sumOp<Type>());
        Foam::reduce(volume, Foam::sumOp<Foam::scalar>());
    }
    return meanField/volume;
}


void SolidCloud::solidFluidInteract(Solid& solid, scalar dt)
{
    m_geotools.clearCache();
    const Foam::vectorField& cc = m_mesh.cellCentres();
    const Foam::scalarField& cv = m_mesh.V();

    scalar dtINV = 1.0/dt;
    vector force  = vector::zero;
    vector torque = vector::zero;

    m_cellenum.SetSolid(solid);

    int numInsideCell = 0;
    int numBorderCell = 0;
    int insideType = solid.getID() + 4;
    scalar alpha = 0.0;
    while (!m_cellenum.Empty())
    {
        int icur = m_cellenum.GetCurCellInd();
        if (m_cellenum.GetCurCellType() == CellEnumerator::ALL_INSIDE)
        {
            ++numInsideCell;
            alpha = 1.0;
            m_ct[icur] = insideType;
        }
        else
        {
            m_ct[icur] = m_cellenum.GetCurCellType();
            ++numBorderCell;
            alpha = m_geotools.calcCellVolume(icur, solid, m_ON_TWOD)/cv[icur];
        }

        vector us = solid.evalPointVelocity(cc[icur]);
        vector uf = m_Uf[icur];
        vector localforce = alpha*cv[icur]*(uf - us)*dtINV;
        force    += localforce;
        torque   += (cc[icur]-solid.getCenter()) ^ localforce;
        m_Fs[icur] += localforce/cv[icur];
        m_As[icur] += alpha;
        m_Ts[icur] += alpha;

        m_cellenum.Next();
    }

    if (Foam::Pstream::parRun())
    {
        Foam::reduce(numInsideCell, Foam::sumOp<Foam::scalar>());
        Foam::reduce(numBorderCell, Foam::sumOp<Foam::scalar>());
    }

    Foam::Info << "Solid " << solid.getID() << " has " << numInsideCell << '/'
               << numBorderCell << " internal/boundary cells.\n";

    force  *= m_rhof;
    torque *= m_rhof;

    if (Foam::Pstream::parRun())
    {
        Foam::reduce(force,  Foam::sumOp<Foam::vector>());
        Foam::reduce(torque, Foam::sumOp<Foam::vector>());
    }
    solid.setFluidForceAndTorque(force, torque);
}

void SolidCloud::interact(scalar time, scalar dt)
{
    // reset solid field, which are source terms for the fluid solver
    m_ct = 0;
    m_As = 0.0;
    m_Fs = Foam::dimensionedVector("zero", Foam::dimAcceleration, Foam::vector::zero);
    m_Ts = Foam::dimensionedScalar("zero", Foam::dimTemperature, 0.0);
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    for (Solid& solid : m_solids)
        solidFluidInteract(solid, dt);

    checkAlpha();

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> t_elapse = duration_cast<duration<double>>(t2 - t1);
    
    if (Foam::Pstream::master())
    {
        std::ostringstream msg;
        msg << "t = " << std::setw(6) << time
            << " [FSI took " << std::left << std::setprecision(3) << std::setw(6)
            << 1000*t_elapse.count() << " ms]";
        LOG(msg.str());
    }

    m_As.correctBoundaryConditions();
    m_Fs.correctBoundaryConditions();
    m_Ts.correctBoundaryConditions();
}

void SolidCloud::addMidEnvironment()
{
    // add env effect (like gravity) at time = t + dt/2
    for (Solid& solid : m_solids)
    {
        scalar rhos = solid.getMaterial()->getRho();
        vector gprime = ((rhos - m_rhof)/rhos) * m_gravity;
        solid.addAcceleration(gprime);
    }
}

void SolidCloud::solidSolidInteract()
{
    m_ptr_ugrid->clear();
    for (const Solid& s : m_solids)
    {
        vector center = s.getCenter();
        m_ptr_ugrid->insert(center.x(), center.y(), center.z(), s.getID());
    }
    static std::vector<CollisionPair> pairs;
    pairs.clear();
    m_ptr_ugrid->generateCollisionPairs(pairs);
    for (CollisionPair& pair : pairs)
        solidSolidCollision(m_solids[pair.first], m_solids[pair.second]);
}


void SolidCloud::solidSolidCollision(Solid& s1, Solid& s2)
{
    scalar cd; vector cP, cN; // geometric values to be calculated
    collisionFunc cfunc = getCollisionFunc(s1.getShape()->getTypeName(), s2.getShape()->getTypeName());
    if (!cfunc)
    {
        // The collision between the two solid shapes is not implemented. This SHOULD BE issued as an warning message to the use. However, doing do greatly bloats the terminal (or log file). As a temporary solution, we simply skip the collision detection.
        // Mar. 1, 2022. Chenguang Zhang
        // Foam::Info << "Warning: collision between shapes " << s1.getShape()->getTypeName() << " and " << s2.getShape()->getTypeName() << " is not implemented. Skip!\n\n";
        return;
    }
    cd = cfunc(s1, s2, cP, cN); 

    vector force, torque;

    /* force law: f = f(cd, cN, material), subject to customization */
    {
        if (cd < 0) return;
        force  = 1e4*cd*cN;
        torque = vector::zero;
    }

    // add collision force and torque to solids
    s1.addForceAndTorque(-force,-torque);
    s2.addForceAndTorque( force, torque);
}

void SolidCloud::evolve(scalar time, scalar dt)
{
    m_time = time;
    static label N_SUBITER = 20;
    if (m_solids.size() == 1) N_SUBITER = 1;
    scalar dt_sub = dt / N_SUBITER;
    for (int i = 0; i < N_SUBITER; ++i)
    {
        // clear all forces
        for (Solid& solid : m_solids)
            solid.clearForceAndTorque();
        // NOW Fn = 0.0

        for (Solid& solid : m_solids)
        {
            solid.addMidFluidForceAndTorque();
        }
        // NOW Fn = F_f

        this->addMidEnvironment();
        // NOW Fn = F_f + m*g

        this->solidSolidInteract();
        // NOW Fn = F_f + F_c + m*g

        for (Solid& solid : m_solids)
        {
            solid.move(time, dt_sub);
        }
    }

    // appear plane has no env force and no solid-solid interactoin TODO

    for (Solid& solid : m_solids)
        solid.storeOldFluidForce();
}

void SolidCloud::checkAlpha() const
{
    forAll(m_mesh.cells(), icell)
    {
        m_As[icell] = std::min(m_As[icell], scalar(1.0));
    }
}

scalar SolidCloud::totalSolidVolume() const
{
    const Foam::scalarField& cell_vol = m_mesh.V();
    return Foam::gSum(m_As*cell_vol);
}

void SolidCloud::saveState()
{
    if (Foam::Pstream::master())
    {
        if (m_timeStepCounter % m_writeFrequency == 0)
        {
            statefile << (*this);
            statefile.flush();
        }
        if (m_ON_MEANFIELD)
        {
            this->writeMeanField();
        }
        ++m_timeStepCounter;
    }
}

std::ostream& operator<<(std::ostream& os, const SolidCloud& sc)
{
    // 3d cases save full  data: 1+18 columns
    // 2d cases save fewer data: 1+ 9 columns (time,xc,yc,ux,uy,fx,fy,ez,omegaz,tz).  
    if (sc.m_ON_TWOD)
    {
        for (const Solid& solid : sc.m_solids)
        {
            os << sc.m_time << ' ';
            write2D(os, solid);
            os << '\n';
        }
    }
    else
    {
        for (const Solid& solid : sc.m_solids)
            os << sc.m_time << ' ' << solid << '\n';
    }
    return os;
}

void SolidCloud::saveRestart(const std::string& filename)
{
    Foam::OFstream os(filename);
    {
        // update content to current time
        dictionary &solids = m_solidDict.subDict("solids");
        forAll(solids, i)
        {
            const Solid& s = m_solids[i];
            dictionary &solid = solids.subDict(solids.toc()[i]);

            vector tmp;

            tmp = s.getCenter(); if (m_ON_TWOD) tmp.z() = 0.0;
            solid.set("pos", tmp);
            tmp = s.getVelocity(); if (m_ON_TWOD) tmp.z() = 0.0;
            solid.set("vel", tmp);

            solid.set("euler", s.getOrientation().eulerAngles(quaternion::XYZ)*180.0/M_PI); // to degree
            tmp = s.getOmega(); if (m_ON_TWOD) {tmp.x() = 0.0;tmp.y()=0.0;}
            solid.set("omega", tmp);
        }
    }
    os << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
 "| =========                 |                                                 |\n"
 "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
 "|  \\\\    /   O peration     | Version:  6.0.0                                 |\n"
 "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n"
 "|    \\\\/     M anipulation  |                                                 |\n"
 "\\*---------------------------------------------------------------------------*/\n"
 "FoamFile\n"
 "{\n"
 "    version     2.0;\n"
 "    format      ascii;\n"
 "    class       dictionary;\n"
 "    object      solidDict;\n"
 "}\n"
 "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    m_solidDict.write(os, false);
    os << "\n// " << Foam::word(GetTimeString()) << "\n"; // append the time stamp
}

const Solid& SolidCloud::operator[](label i) const
{
    return m_solids[i];
}

}
