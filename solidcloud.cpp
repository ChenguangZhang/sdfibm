#include <queue>
#include <cstdlib>
#include <iomanip>
#include <limits>

#include "Pstream.H"
#include "solidcloud.h"
#include "./libshape/vof.h"
#include "./libmotion/motionfactory.h"
#include "./libshape/shapefactory.h"

Foam::scalar calcLineFraction(const Foam::scalar& phia, const Foam::scalar& phib)
{
    if(phia > 0 && phib > 0)
        return 0;
    if(phia < 0 && phib < 0)
        return 1;
    if(phia > 0) // phib < 0
        return -phib/(phia-phib);
    else
        return -phia/(phib-phia);
}

// ctor
SolidCloud::SolidCloud()
{
    if(Foam::Pstream::master())
    {
        logfile.open("cloud.log");
        logfile << GenBanner("INIT SOLIDCLOUD");
    }
    m_solids.reserve(10);
    m_planes.reserve(10);

    //collision_pairs.reserve(100);
    //ptr_uniform_grid = nullptr;
    //collisionFunctionTableInit();

    m_ptr_Mesh = nullptr;
    m_ptr_Uf = nullptr;
    m_ptr_As = nullptr;
    m_ptr_Fs = nullptr;
    m_ptr_Ts = nullptr;
    m_ptr_ct = nullptr;

    // init output file, append if file exists, create if not
    statefile.open("cloud.out", std::fstream::out | std::fstream::app);
    statefile << std::scientific;

    dictionary root(Foam::IFstream("solidDict")());
    if(Foam::Pstream::master())
        logfile << "reading solids from: " << root.name() << '\n';

    // meta information
    const dictionary& meta = root.subDict("meta");
    m_ON_FLUID = Foam::readBool(meta.lookup("on_fluid"));
    m_ON_TWOD  = Foam::readBool(meta.lookup("on_twod"));
    m_gravity = meta.lookup("gravity");

    // output meta information
    if(Foam::Pstream::master())
    {
        logfile << GenBanner("SUMMARY");
        string dim = m_ON_TWOD ? "2D" : "3D";
        string type = m_ON_FLUID? "FSI": "DEM (fluid disabled)";
        logfile << dim << ' ' << type << " simulation with g = " << "("<<m_gravity[0]<<' '<< m_gravity[1]<<' '<< m_gravity[2]<<")\n";
        logfile << "Binary was compiled at " << __DATE__ << ' ' << __TIME__ << '\n';
        logfile << "Simulation starts   at " << GetTimeString() << '\n';
    }
    if (m_ON_TWOD)
        m_vofCalculator = cellFraction2D;
    else
        m_vofCalculator = cellFraction3D;


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
        m_radiusB = -1.0;
        for (int i=0; i < shapes.size(); ++i)
        {
            const dictionary& para = shapes.subDict(shapes.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));
            std::string name = Foam::word(para.lookup("name"));

            m_libshape[name] = ShapeFactory::create(type, para);
            if(m_libshape[name] == nullptr)
                throw("Unrecognized shape type " + type);
            m_radiusB = std::max(m_radiusB, m_libshape[name]->getRadiusB());

            if(Foam::Pstream::master())
                logfile << "[+] " << type << " as " << name << " (" << m_libshape[name]->description() << ")\n";
        }
        // find the maximum bounding raidus
        if(Foam::Pstream::master())
        {
            logfile << "maximum bounding raidus is " << m_radiusB << '\n';
        }

        // read and create motions
        if(Foam::Pstream::master())
        {
            logfile << GenBanner("CREATE: MOTIONS");
            logfile << "--> Available motion types:\n";
            MotionFactory::report(logfile);
            logfile << "--> Used motions:\n";
        }
        const dictionary& motions = root.subDict("motions");
        for (int i=0; i < motions.size(); ++i)
        {
            const dictionary& para = motions.subDict(motions.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));
            std::string name = Foam::word(para.lookup("name"));
            m_libmotion[name] = MotionFactory::create(type, para);
            if(m_libmotion[name] == nullptr)
               throw("Unrecognized motion type " + type);

            if(Foam::Pstream::master())
                logfile << "[+] " << type << " as " << name << " (" << m_libmotion[name]->description() << ")\n";
        }

        // read and create materials
        if(Foam::Pstream::master())
            logfile << GenBanner("CREATE: MATERIALS");
        const dictionary& materials = root.subDict("materials");
        for (int i=0; i < materials.size(); ++i)
        {
            const dictionary& para = materials.subDict(materials.toc()[i]);
            std::string type = Foam::word(para.lookup("type"));

            if(type == "General")
            {
                real rho = readScalar(para.lookup("rho"));
                std::string name = Foam::word(para.lookup("name"));
                m_libmat[name] = new IMaterial(rho);
            }
            else
            {
                throw("Unrecognizable material parameter!");
            }
       }
    }
    catch (const std::exception& e)
    {
        if(Foam::Pstream::master())
            logfile << "Error when creating shape/motion/material!" << e.what() << '\n';
    }

    // create solids
    if(Foam::Pstream::master())
        logfile << GenBanner("CREATE: SOLIDS");

    const dictionary &solids = root.subDict("solids");
    forAll(solids, i)
    {
        const dictionary &solid = solids.subDict(solids.toc()[i]);

        // read pos
        vector pos = solid.lookup("pos");
        if (m_ON_TWOD)
        {
            // sanity check for 2d simulation: the solid must have z = 0
            if(pos.z() !=  0)
            {
                Quit("Solid must has z=0 in 2D simulation, violated by solid number " + std::to_string(i));
            }
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
        if(mot_name!="free")
            s.setMotion(m_libmotion[mot_name]);
        s.setMaterialAndShape(m_libmat[mat_name], m_libshape[shp_name]);
        this->addSolid(s);

        if(Foam::Pstream::master())
            logfile << "Solid " << i << ":" << " motion = " << mot_name
                << ", material = " << mat_name << ", shape = " << shp_name << "\n";
    }
    if(Foam::Pstream::master())
    {
        logfile << "Totally [" << solids.size() << "] solids.\n";
        logfile << GenBanner("END OF INIT");
    }
}

SolidCloud::~SolidCloud()
{
    logfile << "Simulation finished at " << GetTimeString() << '\n';
    delete m_ms; m_ms = nullptr;
    delete m_ptr_ugrid, m_ptr_ugrid = nullptr;
    delete m_ptr_bbox, m_ptr_bbox = nullptr;
    if(isOnFluid())
    {
        delete m_ptr_cs;
        m_ptr_cs = nullptr;
    }
    statefile.close();
    logfile.close();
}

void SolidCloud::initialCorrect()
{
    this->interact(0, 1);
    m_ptr_As->write(); // write to the zero directory
    logfile << "Initial solid volume field written to 0 directory\n";
    for(Solid& solid : m_solids)
        solid.clearForceAndTorque();
}

void SolidCloud::fixInternal()
{
    // after updating solid velocity, fix fluid velocity in cells fully covered by solid
    const Foam::vectorField& cc = m_ptr_Mesh->C().internalField();
    const Foam::volScalarField & ct = (*m_ptr_ct);
    forAll(cc, icell)
    {
        if(ct[icell] >= 4) // if cell totally within a solid
        {
            label id = ct[icell] - 4; // id of the solid that contains this cell
            m_ptr_Uf->operator[](icell) =  m_solids[id].evalPointVelocity(cc[icell]);
        }
    }
    m_ptr_Uf->correctBoundaryConditions();
}

void SolidCloud::solidFluidInteract(Solid& solid, const real& dt)
{
    // interaction of a single solid with fluid
    const Foam::fvMesh&       mesh = *m_ptr_Mesh;
    const Foam::volVectorField& Uf = *m_ptr_Uf;
    const Foam::pointField &    pp = mesh.points();
    const Foam::vectorField&    cc = mesh.cellCentres();
    const Foam::scalarField&    cv = mesh.V();

    const Foam::vectorField& faceCentres = mesh.faceCentres();
    const Foam::vectorField& faceAreas   = mesh.faceAreas();
    // connectivity
    const Foam::labelListList& c2c = mesh.cellCells();
    const Foam::labelListList& c2p = mesh.cellPoints();

    Foam::volScalarField& ct = *m_ptr_ct; // cell type field

    Foam::volScalarField& As = *m_ptr_As;
    Foam::volVectorField& Fs = *m_ptr_Fs;
    Foam::volScalarField& Ts = *m_ptr_Ts;

    // handle host cell
    int hostid = m_ms->findCell(solid.getCenter());
    if(hostid < 0)
    {
        // fallback to more expensive option, scan cells of the current partition to find the hostid
        forAll(cc, icell)
        {
            if(solid.isInside(cc[icell]))
            {
                int insideCount = 0;
                // loop cell vertices
                forAll(c2p[icell], ivert)
                    insideCount += solid.isInside(pp[c2p[icell][ivert]]);
                if(insideCount == c2p[icell].size()) // equals #vertex
                {
                    hostid = icell;
                    break;
                }
            }
        }
    }

    // force and torque on the solid by the current partition
    vector force  = vector::zero;
    vector torque = vector::zero;

    if(hostid>0) // solid overlaps partition
    {
        label innerType = 4 + solid.getID();
        ct[hostid] = innerType;

        // use breadth first search (BFS) to visit cells intersected by the solid
        std::queue<int> q;
        q.push(hostid);
        std::vector<int> intersectCell;

        real dtINV = 1.0/dt;
        while(!q.empty())
        {
            // icur is the index of the current cell
            int icur = q.front();
            q.pop();

            // process current cell
            if(ct[icur] == innerType)
            {
                // cell is fully within solid
                real alpha = 1.0;
                vector us = solid.evalPointVelocity(cc[icur]);
                vector uf = Uf[icur];
                real   dV = cv[icur];

                vector localforce = alpha*dV*(uf - us)*dtINV;

                force += localforce;
                torque+= (cc[icur]-solid.getCenter()) ^ localforce;

                Fs[icur] = localforce/cv[icur]; // note division by cell volume
                As[icur] = alpha;
                Ts[icur] = alpha;
            }

            // loop neighbors of the current cell
            forAll(c2c[icur], itmp)
            {
                Foam::label inb = c2c[icur][itmp];
                if(ct[inb] == 0) // unvisited nb cell
                {
                    // detailed check of this nb cell
                    int insideCount = 0;
                    // loop cell vertices
                    forAll(c2p[inb], ivert)
                    {
                        insideCount += solid.isInside(pp[c2p[inb][ivert]]);
                    }

                    if(insideCount == 0)
                    {
                        // cell totally outside of solid, no need to proceed, but we still need to mark it as "visited"
                        ct[inb] = 1;
                        continue;
                    }

                    // at this point the cell is both unvisited & intersected for sure, so we add it to the queue
                    q.push(inb);

                    if(insideCount == c2p[inb].size())
                    {
                        ct[inb] = innerType;
                    }
                    else
                    {
                        intersectCell.push_back(inb);
                        // check cell center
                        if(solid.isInside(cc[inb]))
                            ct[inb] = 3; // cell intersected and has centroid inside solid
                        else
                            ct[inb] = 2; // cell intersected but has centroid outside of solid
                    }
                }
            }
        }

        for(label icur : intersectCell)
        {
            // build a hash map first to store the sdf of all vertices
            std::map<int, Foam::scalar> sdf_values; // vertex-id : value pairs
            const Foam::cell& cell = mesh.cells()[icur];
            const Foam::labelList& pointids = c2p[icur];
            forAll(pointids, ipoint)
            {
                sdf_values[pointids[ipoint]] = solid.getShape()->signedDistance(pp[pointids[ipoint]], solid.getCenter(), solid.getOrientation());
            }

            Foam::vector pcell;
            const Foam::vector& A = pp[pointids[0]];
            Foam::scalar phiA = sdf_values[pointids[0]];
            int i;
            Foam::vector B = Foam::vector::zero;
            Foam::scalar phiB = 0.0;
            for(i = 1; i < pointids.size(); ++i)
            {
                B = pp[pointids[i]];
                phiB = sdf_values[pointids[i]];
                if (phiA * phiB < 0)
                    break;
            }
            pcell = A - std::fabs(phiA)/(std::fabs(phiA)+std::fabs(phiB))*(A-B);
            if(m_ON_TWOD)
                pcell[2] = 0.0;

            Foam::scalar volume = 0.0;
            forAll(cell, iface)
            {
                Foam::label faceid = cell[iface];
                Foam::face myface = mesh.faces()[faceid];

                Foam::label nvertex = myface.size();
                std::vector<Foam::scalar> phis(nvertex);
                forAll(myface, ivertex)
                {
                    phis[ivertex] = sdf_values[myface[ivertex]];
                }

                Foam::scalar epsilon_f = 0.0;

                int sign_sum = 0;
                forAll(myface, ivertex)
                {
                    if(phis[ivertex] > 0)
                        sign_sum += 1;
                    else
                        sign_sum -= 1;
                }

                if(sign_sum == nvertex)
                    epsilon_f = 0.0;
                else if(sign_sum ==-nvertex)
                    epsilon_f = 1.0;
                else
                {
                    const Foam::vector& A = pp[myface[0]];
                    Foam::scalar phiA = sdf_values[myface[0]];

                    Foam::vector B = Foam::vector::zero;
                    Foam::scalar phiB = 0.0;
                    int i;
                    for(i = 1; i < nvertex; ++i)
                    {
                        B = pp[myface[i]];
                        phiB = sdf_values[myface[i]];
                        if (phiA * phiB < 0)
                            break;
                    }

                    Foam::vector pface = A - phiA/(phiA-phiB)*(A-B);

                    Foam::scalar area = 0.0;
                    for(int iseg = 0; iseg < nvertex; ++iseg)
                    {
                        const Foam::scalar& phiO = phis[iseg];
                        const Foam::scalar& phiA = phis[(iseg+1)%nvertex];
                        const Foam::vector& O = pp[myface[iseg]];
                        const Foam::vector& A = pp[myface[(iseg+1)%nvertex]];
                        area += std::fabs(0.5*Foam::mag((A-O) ^ (pface-O))) * calcLineFraction(phiO, phiA);
                    }
                    epsilon_f = area/Foam::mag(faceAreas[faceid]);
                }
                volume += (1.0/3.0)*epsilon_f*std::fabs((pcell - faceCentres[faceid]) & faceAreas[faceid]);
            }
            Foam::scalar alpha = volume/mesh.V()[icur];

            // below is identical for 2D and 3D
            vector us = solid.evalPointVelocity(cc[icur]);
            vector uf = Uf[icur];
            real   dV = cv[icur];

            vector localforce = alpha*dV*(uf - us)*dtINV;

            force += localforce;
            torque+= (cc[icur]-solid.getCenter()) ^ localforce;

            Fs[icur] = localforce/cv[icur];
            As[icur] = alpha;
            Ts[icur] = alpha;
        }

        force  *= m_rho;
        torque *= m_rho;
    }

    // a solid can cover several parallel partitions, with each contributing a "partial" force
    // they are added together
    if(Foam::Pstream::parRun())
    {
        // init *PerCPU to be local copies
        Foam::vector forcePerCPU ((force));
        Foam::vector torquePerCPU((torque));
        // parallel reduce
        Foam::reduce(forcePerCPU , Foam::sumOp<Foam::vector>());
        Foam::reduce(torquePerCPU, Foam::sumOp<Foam::vector>());
        // now *PerCPU contain the global sum
        force = forcePerCPU;
        torque= torquePerCPU;
    }
    solid.setFluidForceAndTorque(force, torque);
}

void SolidCloud::linkMesh(const Foam::fvMesh& mesh)
{
    m_ON_RESTART = false;
    m_ptr_Mesh= &(const_cast<Foam::fvMesh&>(mesh));
    // create bounding box
    Foam::vector low  = mesh.bounds().min();
    Foam::vector high = mesh.bounds().max();
    real a[3] = {low.x(), low.y(), low.z()};
    real b[3] = {high.x(),high.y(),high.z()};
    m_ptr_bbox = new BBox(a, b);
    m_ptr_bbox->report(logfile);
    m_ptr_ugrid = new UGrid(*m_ptr_bbox, 2*m_radiusB);

    if(mesh.time().value() > 0)
    {
        std::string filename;
        if(!Foam::Pstream::parRun())
            filename = "./" + mesh.time().timeName() + "/restart";
        else
            // when parallel, all init from the master
            filename = "./processor0/" + mesh.time().timeName() + "/restart";
        logfile << "Restart from " << filename << std::endl;
        loadRestart(filename);
        m_ON_RESTART = true;
    }
    else
    {
        std::cout << "New simulation, starting from 0\n";
    }

}

void SolidCloud::linkFluid(const Foam::volVectorField& U)
{
    const Foam::fvMesh& mesh = U.mesh();
    m_ms = new Foam::meshSearch(mesh);

    m_ptr_Uf  = &(const_cast<Foam::volVectorField&>(mesh.lookupObject<Foam::volVectorField>("U")));
    m_ptr_As  = &(const_cast<Foam::volScalarField&>(mesh.lookupObject<Foam::volScalarField>("As")));
    m_ptr_ct  = &(const_cast<Foam::volScalarField&>(mesh.lookupObject<Foam::volScalarField>("Ct")));
    m_ptr_Fs  = &(const_cast<Foam::volVectorField&>(mesh.lookupObject<Foam::volVectorField>("Fs")));
    m_ptr_Asf = &(const_cast<Foam::surfaceScalarField&>(mesh.lookupObject<Foam::surfaceScalarField>("Asf")));
    m_ptr_Ts  = &(const_cast<Foam::volScalarField&>(mesh.lookupObject<Foam::volScalarField>("Ts")));

    // created and destroyed by SolidCloud, not output
    m_ptr_cs = new Foam::scalarField(mesh.nCells(), 0.0);

    const Foam::dictionary& transportProperties = mesh.lookupObject<Foam::IOdictionary>("transportProperties");
    const Foam::dimensionedScalar&      density = transportProperties.lookup("rho");
    m_rho = density.value();

    // precalculate the mesh cell size field, must after creating solid cloud
    if(m_ON_TWOD)
    {
        forAll(mesh.V(), icell)
        {
            m_ptr_cs->operator[](icell) = std::sqrt(real(mesh.V()[icell]));
        }
    }
    else
    {
        forAll(mesh.V(), icell)
        {
            m_ptr_cs->operator[](icell) = std::cbrt(real(mesh.V()[icell]));
        }
    }

    // correct initial condtion to reduce the number of iteration at the 1st time step
    if(!m_ON_RESTART)
        initialCorrect();
}

void SolidCloud::interact(const real& time, const real& dt)
{
    // reset solid field, which are source terms for the fluid solver
    (*m_ptr_ct) = 0;
    (*m_ptr_As) = 0.0;
    (*m_ptr_Asf) = 0.0;
    (*m_ptr_Fs) = Foam::dimensionedVector("zero", Foam::dimAcceleration, Foam::vector::zero);

    (*m_ptr_Ts) = Foam::dimensionedScalar("zero", Foam::dimTemperature, 0.0);

    for(Solid& solid : m_solids)
        solidFluidInteract(solid, dt);

    for(Solid& solid : m_planes)
        solidFluidInteract(solid, dt);

    m_ptr_As->correctBoundaryConditions();
    m_ptr_Fs->correctBoundaryConditions();
    m_ptr_Ts->correctBoundaryConditions();
}

void SolidCloud::addMidEnvironment()
{
    // add env effect (like gravity) at time = t + dt/2 (only to solid, not to plane)
    for(Solid& solid : m_solids)
    {
        real rhos = solid.getMaterial()->getRho();
        vector gprime = ((rhos - m_rho)/rhos) * m_gravity;
        solid.addAcceleration(gprime);
    }
}

void SolidCloud::solidSolidInteract()
{
    m_ptr_ugrid->clear();
    for(const Solid& s : m_solids)
    {
        label id = s.getID();
        vector center = s.getCenter();
        m_ptr_ugrid->insert(center.x(), center.y(), center.z(), id);
    }
    std::vector<CollisionPair> pairs;
    m_ptr_ugrid->generateCollisionPairs(pairs);

    for(CollisionPair& pair : pairs)
    {
        // hardcoded parameters
        static real crange = 0.05;
        static real epsilon = 1.0e-3;
        static real msg = 2.710932;

        Solid& s1 = m_solids[pair.first];
        Solid& s2 = m_solids[pair.second];

        scalar c2c = Foam::mag(s1.getCenter() - s2.getCenter());
        scalar contact_depth = (crange+s1.getShape()->getRadiusB()+s2.getShape()->getRadiusB()-c2c)/crange;

        if(contact_depth < 0) continue;
        vector contact_normal = Foam::normalised(s1.getCenter() - s2.getCenter());
        vector force = msg/epsilon*contact_depth*contact_depth*contact_normal;
        vector torque = vector::zero;

        s1.addForceAndTorque( force, torque);
        s2.addForceAndTorque(-force,-torque);
    }
}

void SolidCloud::evolve(const real& time, const real& dt)
{
    static label N_SUBITER = 20;
    if (m_solids.size() == 1) N_SUBITER = 1;
    real dt_sub = dt / N_SUBITER;
    for(int i = 0; i < N_SUBITER; ++i)
    {
        // clean all forces
        for(Solid& solid : m_solids)
            solid.clearForceAndTorque();
        // NOW Fn = 0.0

        for(Solid& solid : m_solids)
        {
            solid.addMidFluidForceAndTorque();
        }
        // NOW Fn = F_f

        this->addMidEnvironment();
        // NOW Fn = F_f + m*g

        this->solidSolidInteract();
        // NOW Fn = F_f + F_c + m*g

        for(Solid& solid : m_solids)
        {
            solid.move(time, dt_sub);
        }
    }

    {
        // clean all forces
        for(Solid& solid : m_planes)
            solid.clearForceAndTorque();
        // NOW Fn = 0.0

        for(Solid& solid : m_planes)
        {
            solid.addMidFluidForceAndTorque();
        }
        // NOW Fn = F_f

        this->addMidEnvironment();
        // NOW Fn = F_f + m*g

        for(Solid& solid : m_planes)
        {
            solid.move(time, dt);
        }
    }

    for(Solid& solid : m_solids)
        solid.storeOldFluidForce();

    // store state
    if(Foam::Pstream::master())
        saveState(time);
}

void SolidCloud::checkAlpha() const
{
    forAll(m_ptr_Mesh->cells(), icell)
    {
        real alpha = m_ptr_As->operator[](icell);
        if(alpha < 0 || alpha > 1)
        {
           Foam::Info << "Unbounded cell volume fraction!\n";
           Foam::Info << "As[" << icell << "] = " << alpha << '\n';
           Quit("Quit\n");
        }
    }
}

real SolidCloud::totalSolidVolume() const
{
    const Foam::scalarField& cell_vol = m_ptr_Mesh->V();
    return Foam::gSum((*m_ptr_As)*cell_vol);
}

void SolidCloud::saveState(const real& time)
{
    vector v;
    for(Solid& solid : m_solids)
    {
        statefile << time << ' ';
        v = solid.getCenter();  statefile << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
        v = solid.getVelocity();statefile << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
        v = solid.getForce();   statefile << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';

        v = solid.getOrientation().eulerAngles(quaternion::XYZ);
                                statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getOmega();   statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getTorque();  statefile << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
    }

    for(Solid& solid : m_planes)
    {
        statefile << time << ' ';
        v = solid.getCenter();  statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getVelocity();statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getForce();   statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';

        v = solid.getOrientation().eulerAngles(quaternion::XYZ);
                                statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getOmega();   statefile << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getTorque();  statefile << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
    }
    statefile.flush();
}

void SolidCloud::saveRestart(const std::string& filename)
{
    int ndigits = std::numeric_limits<real>::max_digits10;
    std::ofstream fout(filename);
    if(!fout.good()){
        Quit("Can't open file " + filename + "\n!");
    }
    fout.precision(ndigits);
    vector v;
    for(Solid& solid : m_solids)
    {
        v = solid.getCenter();  fout << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
        v = solid.getVelocity();fout << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
        v = solid.getForce();   fout << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';

        v = solid.getOrientation().eulerAngles(quaternion::XYZ);
                                fout << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getOmega();   fout << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getTorque();  fout << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getFluidForce();  fout << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
        v = solid.getFluidTorque();  fout << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
    }

    fout.close();
}

void SolidCloud::loadRestart(const std::string& filename)
{
    std::ifstream fin(filename, std::ifstream::in);
    if(!fin.good()){
        Quit("Can't open file " + filename + "\n!");
    }

    real x, y, z;
    for(Solid& solid : m_solids)
    {
        fin >> x >> y >> z; solid.setCenter     (vector(x, y, z));
        fin >> x >> y >> z; solid.setVelocity   (vector(x, y, z));
        fin >> x >> y >> z; solid.setForce      (vector(x, y, z));
        fin >> x >> y >> z; solid.setOrientation(vector(x, y, z));
        fin >> x >> y >> z; solid.setOmega      (vector(x, y, z));
        fin >> x >> y >> z; solid.setTorque     (vector(x, y, z));
        fin >> x >> y >> z; solid.setFluidForceOld(vector(x, y, z));
        fin >> x >> y >> z; solid.setFluidTorqueOld(vector(x, y, z));
    }
    fin.close();
}

const Solid& SolidCloud::operator[](label i) const
{
    return m_solids[i];
}
