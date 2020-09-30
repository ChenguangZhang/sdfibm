#include <cstdlib>
#include <iomanip>
#include <limits>

#include "fvCFD.H"

#include "solidcloud.h"
#include "./libmotion/motionfactory.h"
#include "./libshape/shapefactory.h"
namespace sdfibm {

SolidCloud::SolidCloud(const Foam::word& dictfile, const Foam::fvMesh& mesh):
    m_mesh(mesh),
    m_geotools(GeometricTools(m_mesh)),
    m_cellenum(CellEnumerator(m_mesh))
{
    if(Foam::Pstream::master())
    {
        logfile.open("cloud.log");
        logfile << GenBanner("INIT SOLIDCLOUD");
    }
    m_solids.reserve(10);
    m_planes.reserve(10);

    Foam::IFstream ifstream = Foam::IFstream(dictfile);
    dictionary root(ifstream());

    if(Foam::Pstream::master())
        logfile << "reading solids from: " << root.name() << '\n';

    const dictionary& meta = root.subDict("meta");
    m_ON_TWOD  = Foam::readBool(meta.lookup("on_twod"));

    if(Foam::Pstream::master())
    {
        logfile << GenBanner("SUMMARY");
        string dim = m_ON_TWOD ? "2D" : "3D";
        logfile << "Binary was compiled at " << __DATE__ << ' ' << __TIME__ << '\n';
        logfile << "Simulation starts   at " << GetTimeString() << '\n';
    }

   // read and create shape, motion, and material
   try {
        // read and create shapes
        if(Foam::Pstream::master())
        {
            logfile << GenBanner("CREATE: SHAPES");
            logfile << "--> Available shape types:\n";
            ShapeFactory::report(logfile);
            logfile << "--> Used shapes:\n";
        }
        const dictionary& shapes = root.subDict("shapes");
        for (int i=0; i < shapes.size(); ++i)
        {
            const dictionary& para = shapes.subDict(shapes.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));
            std::string name = Foam::word(para.lookup("name"));

            m_libshape[name] = ShapeFactory::create(type, para);
            if(m_libshape[name] == nullptr)
                throw std::runtime_error(std::string("Unrecognized shape type " + type + '\n'));

            if(Foam::Pstream::master())
                logfile << "[+] " << type << " as " << name << " (" << m_libshape[name]->description() << ")\n";
        }

        m_libmat["mat"] = new IMaterial(1.0);
    }
    catch (const std::exception& e)
    {
        if(Foam::Pstream::master())
        {
            std::cout << e.what();
            logfile << "Error when creating shape/motion/material!" << e.what() << '\n';
        }
        std::exit(1);
    }

    // create solids
    if(Foam::Pstream::master())
        logfile << GenBanner("CREATE: SOLIDS & PLANES");

    const dictionary &solids = root.subDict("solids");
    forAll(solids, i)
    {
        const dictionary &solid = solids.subDict(solids.toc()[i]);

        vector pos = solid.lookup("pos");
        if (m_ON_TWOD)
        {
            // sanity check for 2d simulation: solid must have z = 0
            if(pos.z() !=  0)
                Quit("Solid must has z=0 in 2D simulation, violated by solid # " + std::to_string(i));
        }
        // create solid
        Solid s(i, pos, quaternion::I);
        // read velocity, euler angle, angular velocity
        vector euler = solid.lookupOrDefault("euler", vector::zero);
        s.setOrientation(euler*M_PI/180.0);

        std::string shp_name = Foam::word(solid.lookup("shp_name"));
        s.setMaterialAndShape(m_libmat["mat"], m_libshape[shp_name]);
        this->addSolid(s);

        if(Foam::Pstream::master())
            logfile << "Solid " << i << " shape = " << shp_name << "\n";
    }

    const dictionary &planes = root.subDict("planes");
    forAll(planes, i)
    {
        const dictionary &plane = planes.subDict(planes.toc()[i]);

        vector pos = plane.lookup("pos");
        if (m_ON_TWOD)
        {
            // sanity check for 2d simulation: plane must have z = 0
            if(pos.z() !=  0)
                Quit("Plane must has z=0 in 2D simulation, violated by plane # " + std::to_string(i));
        }
        // create plane
        Solid s(i, pos, quaternion::I);
        // read velocity, euler angle, angular velocity
        vector euler = plane.lookupOrDefault("euler", vector::zero);
        s.setOrientation(euler*M_PI/180.0);

        std::string shp_name = Foam::word(plane.lookup("shp_name"));
        s.setMaterialAndShape(m_libmat["mat"], m_libshape[shp_name]);
        this->addPlane(s);

        if(Foam::Pstream::master())
            logfile << "Plane " << i << " shape = " << shp_name << "\n";
    }

    if(Foam::Pstream::master())
    {
        logfile << "Totally [" << m_solids.size() << "] solids and [" << m_planes.size() << "] planes.\n";
        logfile << GenBanner("END OF INIT");
    }
}


void SolidCloud::solidFluidInteract(Solid& solid, const scalar& dt)
{
    m_geotools.clearCache();
    m_cellenum.SetSolid(solid);

    const Foam::scalarField& cv = m_mesh.V();

    scalar alpha = 0.0;
    while (!m_cellenum.Empty())
    {
        int icur = m_cellenum.GetCurCellInd();
        if (m_cellenum.GetCurCellType() == CellEnumerator::ALL_INSIDE)
            alpha = 1.0;
        else
            alpha = m_geotools.calcCellVolume(icur, solid, m_ON_TWOD)/cv[icur];

        m_ptr_As->operator[](icur) += alpha;
        if (m_ptr_As->operator[](icur) > 1.0)
            m_ptr_As->operator[](icur) = 1.0;

        m_cellenum.Next();
    }
}

void SolidCloud::writeVOF(const Foam::word& name)
{
    dimensionedScalar zero("zero", Foam::dimless, 0.0);
    volScalarField     tmp( IOobject("tmp", "0", m_mesh, IOobject::NO_READ,   IOobject::NO_WRITE  ), m_mesh, zero);
    surfaceScalarField phi( IOobject("phi", "0", m_mesh, IOobject::NO_READ,   IOobject::NO_WRITE  ), m_mesh, zero);
    volScalarField   field( IOobject( name, "0", m_mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), m_mesh);

    m_ptr_As = &tmp;
    for(Solid& solid : m_solids)
        solidFluidInteract(solid, 0.0);

    for(Solid& solid : m_planes)
        solidFluidInteract(solid, 0.0);

    tmp.correctBoundaryConditions();

    forAll(m_mesh.cells(), icell)
    {
        scalar alpha = tmp[icell];
        if(alpha < 0 || alpha > 1 + 1e-6)
            Foam::Info << "Wrong volume fraction " << alpha << " at cell " << icell << '\n';
        field[icell] = alpha;
    }
    field.correctBoundaryConditions();

    field.write();
    const Foam::scalarField& cell_vol = m_mesh.V();
    Foam::Info << "* Write smooth phase field [" << name << "] to ./0, total volume = "
               << Foam::gSum(tmp*cell_vol) << '\n';
    Foam::Info << "* Finished at " << sdfibm::GetTimeString() << '\n';

    logfile.close();
}

}
