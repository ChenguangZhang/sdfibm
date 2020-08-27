#ifndef GEOMETRICTOOLS_H
#define GEOMETRICTOOLS_H

#include <vector>
#include <unordered_map>
#include "meshinfo.h"
#include "types.h"

namespace sdfibm {

class Solid;

class GeometricTools : public MeshInfo
{
    typedef std::unordered_map<int, scalar> CacheMap; // vertexInd-phi
public:
    GeometricTools(const Foam::fvMesh& mesh) : MeshInfo(mesh) {}

    scalar UpdateCache(label vertexInd, const Solid& solid);

    vector calcApex(const Foam::labelList& vertexInds, CacheMap& phis) const;
    scalar calcFaceArea(const Foam::face& vertexInds, CacheMap& phis);
    scalar calcLineFraction(const scalar& phia, const scalar& phib) const;
    scalar calcFaceAreaFraction(const Foam::face& vertexInds, CacheMap& phis, label faceInd);
    scalar calcCellVolume(label cellInd, const Solid& solid, bool isTWOD);
    inline void clearCache() {phiCache.clear();};
private:
    CacheMap phiCache;
};

}
#endif // GEOMETRICTOOLS_H
