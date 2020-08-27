#ifndef MESHINFO_H
#define MESHINFO_H

#include "fvc.H"

namespace sdfibm {

class MeshInfo
{
protected:
    const Foam::fvMesh& m_mesh; // cell to cell connectivity
    const Foam::labelListList& m_c2c; // cell to cell connectivity
    const Foam::labelListList& m_c2p; // cell to point connectivity
    const Foam::pointField   & m_pp; // mesh points
    const Foam::vectorField  & m_cc; // mesh centers
    const Foam::scalarField  & m_cv; // cell volumes
    const Foam::vectorField  & m_fc; // face centers
    const Foam::vectorField  & m_fa; // face areas
public:
    MeshInfo(const Foam::fvMesh& mesh):
        m_mesh(mesh),
        m_c2c(mesh.cellCells()),
        m_c2p(mesh.cellPoints()),
        m_pp(mesh.points()),
        m_cc(mesh.cellCentres()),
        m_cv(mesh.V()),
        m_fc(mesh.faceCentres()),
        m_fa(mesh.faceAreas())
    {}
};

}
#endif // MESHINFO_H
