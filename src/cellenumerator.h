#ifndef CELLENUMERATOR_H
#define CELLENUMERATOR_H

#include <vector>
#include <queue>
#include "meshinfo.h"
#include "meshSearch.H"

// [08-27-2020] Thanks to Mr. Haorang Wang for the idea of using an enumerator class, and for
// the code from which the currently class is adapted.

namespace sdfibm {

class Solid;

class CellEnumerator : public MeshInfo
{
public:
    enum CELL_TYPE {UNVISITED, ALL_INSIDE, CENTER_INSIDE, CENTER_OUTSIDE, ALL_OUTSIDE};
private:
    static std::vector<CELL_TYPE> m_ct; // cell type, size = #mesh cell
    Foam::meshSearch* m_ms;

    std::queue<int> m_queue;
    const Solid* mp_solid;
    int m_seed;

    void _next();

    /* count the number of vertices of a cell that falls within the solid
     * @param icell The index of the mesh cell considered
     * @param solid The solid considered, its shape information is used to check each vertex of the cell
     */
    int CountVertexInside(int icell, const Solid& solid) const;

public:
    CellEnumerator(const Foam::fvMesh& mesh);
    ~CellEnumerator();
    void SetSolid(const Solid& solid);

    inline void Next()
    {
        _next();
        m_queue.pop();
    }


    inline int GetCurCellInd() const
    {
        return m_queue.front();
    }

    inline CELL_TYPE GetCurCellType() const
    {
        return m_ct[m_queue.front()];
    }

    inline bool Empty() const
    {
        return m_queue.empty();
    }
};

}
#endif // CELLENUMERATOR_H
