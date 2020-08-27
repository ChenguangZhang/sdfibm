#include "cellenumerator.h"
#include "solid.h"

namespace sdfibm {

void CellEnumerator::_next()
{
    int icur = m_queue.front();
    forAll(m_c2c[icur], _)
    {
        label inb = m_c2c[icur][_];
        if (m_ct[inb] == CELL_TYPE::UNVISITED) // unvisited
        {
            int n_v_inside = CountVertexInside(inb, *mp_solid);

            if (n_v_inside==0)
            {
                m_ct[inb] = CELL_TYPE::ALL_OUTSIDE;
                continue;
            }

            // at this point the cell is both unvisited and intersected
            m_queue.push(inb);
            if (n_v_inside == m_c2p[icur].size())
                m_ct[inb] = CELL_TYPE::ALL_INSIDE;
            else if (mp_solid->isInside(m_cc[inb]))
                m_ct[inb] = CELL_TYPE::CENTER_INSIDE;
            else
                m_ct[inb] = CELL_TYPE::CENTER_OUTSIDE;
        }
    }
}

/* count the number of vertices of a cell that falls within the solid
 * @param icell The index of the mesh cell considered
 * @param solid The solid considered, it checks each cell vertex for containment
 */
int CellEnumerator::CountVertexInside(int icell, const Solid& solid) const
{
    int n_v_inside = 0;
    forAll(m_c2p[icell], _)
    {
        int ivert = m_c2p[icell][_];
        const vector& p = m_pp[ivert];
        n_v_inside += solid.isInside(p);
    }
    return n_v_inside;
}

CellEnumerator::CellEnumerator(const Foam::fvMesh& mesh) : MeshInfo(mesh)
{
    m_ms = new Foam::meshSearch(mesh);
    m_ct.resize(m_cv.size(), CELL_TYPE::UNVISITED);
}

void CellEnumerator::SetSolid(const Solid& solid)
{
    std::fill(m_ct.begin(), m_ct.end(), UNVISITED);

    mp_solid = &solid;

    m_seed = m_ms->findNearestCell(mp_solid->getCenter());
    if (!m_mesh.pointInCell(mp_solid->getCenter(), m_seed))
    {
        int insideCount = 0;
        forAll(m_c2p[m_seed], ivert)
            insideCount += mp_solid->isInside(m_pp[m_c2p[m_seed][ivert]]);
        if (insideCount != m_c2p[m_seed].size())
            m_seed = -1;
    }

    // fallback to the most costly option: scan all cells of current partition to find seed
    if (m_seed == -1)
    {
        forAll(m_cc, icell)
        {
            if (mp_solid->isInside(m_cc[icell]))
            {
                int insideCount = 0;
                // loop cell vertices
                forAll(m_c2p[icell], ivert)
                    insideCount += mp_solid->isInside(m_pp[m_c2p[icell][ivert]]);
                if (insideCount == m_c2p[icell].size())
                {
                    m_seed = icell;
                    break;
                }
            }
        }
    }

    if (m_seed >= 0)
    {
        m_queue.push(m_seed);
        m_ct[m_seed] = CELL_TYPE::ALL_INSIDE;

        Foam::Info << "Solid " << solid.getID() <<  " has seed " << m_seed << '\n';
    }
}

CellEnumerator::~CellEnumerator()
{
    delete m_ms;
    m_ms = nullptr;
}

std::vector<CellEnumerator::CELL_TYPE> CellEnumerator::m_ct;

}
