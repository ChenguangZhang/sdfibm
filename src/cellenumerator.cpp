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

    int seed = m_ms->findNearestCell(mp_solid->getCenter());
    if (seed < 0)
    {
        forAll(m_cc, icell)
        {
            if (seed > 0) break;
            forAll(m_c2p[icell], ivert)
            {
                if (mp_solid->isInside(m_pp[m_c2p[icell][ivert]]))
                {
                    seed = icell;
                    break;
                }
            }
        }
    }

    if (seed >= 0)
    {
        m_queue.push(seed);

        if (CountVertexInside(seed, *mp_solid) == m_c2p[seed].size())
            m_ct[seed] = CELL_TYPE::ALL_INSIDE;
        else
            if (mp_solid->isInside(m_cc[seed]))
                m_ct[seed] = CELL_TYPE::CENTER_INSIDE;
            else
                m_ct[seed] = CELL_TYPE::CENTER_OUTSIDE;
    }
}

CellEnumerator::~CellEnumerator()
{
    delete m_ms;
    m_ms = nullptr;
}

std::vector<CellEnumerator::CELL_TYPE> CellEnumerator::m_ct;

}
