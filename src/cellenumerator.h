#ifndef CELLENUMERATOR_H
#define CELLENUMERATOR_H

#include <vector>
#include <queue>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include "meshinfo.h"
#include "meshSearch.H"
#include "types.h"

// [08-27-2020] Thanks to Mr. Haorang Wang for the idea of using an enumerator class, and for
// the code from which the currently class is adapted.

namespace sdfibm {

class CellEnumerator : public MeshInfo
{
public:
    enum CELL_TYPE {UNVISITED, ALL_INSIDE, CENTER_INSIDE, CENTER_OUTSIDE, ALL_OUTSIDE};
    using Predicate = std::function<bool(const vector&)>;
    struct IntersectionSet
    {
        std::set<size_t> all_inside_cells;
        std::set<size_t> center_inside_cells;
        std::set<size_t> center_outside_cells;
        std::set<size_t> all_outside_cells;
    };
private:
    std::vector<CELL_TYPE> m_ct; // cell type, size = #mesh cell
    Predicate pred_;

    std::queue<int> m_queue;

    void _next();

    int CountVertexInside(int icell, const Predicate& p) const;

public:
    CellEnumerator(const Foam::fvMesh& mesh, const Predicate& pred, int seed = -1); // TODO, use optional
    ~CellEnumerator() = default;

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
