#include "geometrictools.h"
#include "solid.h"

namespace sdfibm {

scalar GeometricTools::UpdateCache(label vertexInd, const Solid& solid)
{
    if (phiCache.find(vertexInd) == phiCache.end())
        phiCache[vertexInd] = solid.phi(m_pp[vertexInd]);
    return phiCache[vertexInd];
}

scalar GeometricTools::calcLineFraction(const scalar& phia, const scalar& phib) const
{
    if (phia > 0 && phib > 0)
        return 0;
    if (phia <= 0 && phib <= 0)
        return 1;
    if (phia > 0) // phib < 0
        return -phib/(phia-phib);
    else
        return -phia/(phib-phia);
}

vector GeometricTools::calcApex(const Foam::labelList& vertexInds, CacheMap& phis) const
{
    label nvertex = vertexInds.size(); // #face vertex

    const vector& A = m_pp[vertexInds[0]];
    scalar     phiA = phis[vertexInds[0]];

    vector B = vector::zero;
    scalar phiB = 0.0;
    label i;
    for (i = 1; i < nvertex; ++i)
    {
        B    = m_pp[vertexInds[i]];
        phiB = phis[vertexInds[i]];

        if (phiA * phiB <= 0)
            break;
    }

    return A - std::abs(phiA)/(SMALL + std::fabs(phiA)+std::fabs(phiB))*(A-B);
}

scalar GeometricTools::calcCellVolume(label cellInd, const Solid& solid, bool isTWOD = false)
{
    const Foam::labelList& vertexInds = m_c2p[cellInd];
    // prepare cell's vertex phi
    forAll(vertexInds, ivertex)
    {
        UpdateCache(vertexInds[ivertex], solid);
    }

    vector apex = calcApex(vertexInds, phiCache);
    if (isTWOD)
        apex[2] = 0.0;

    scalar volume = 0.0;
    const Foam::cell& faceInds = m_mesh.cells()[cellInd]; // cell is a labelList of faces
    forAll(faceInds, iface)
    {
        // visit each cell face (face is labelList)
        label faceInd = faceInds[iface];
        Foam::face myface = m_mesh.faces()[faceInd];
        Foam::scalar eps_f = calcFaceAreaFraction(myface, phiCache, faceInd);

        volume += (1.0/3.0)*eps_f*std::fabs((apex - m_fc[faceInd]) & m_fa[faceInd]);
    }
    return volume;
}

scalar GeometricTools::calcFaceArea(const Foam::face& vertexInds, CacheMap& phis)
{
    label nvertex = vertexInds.size(); // #face vertex

    vector apex = calcApex(vertexInds, phis);

    std::vector<scalar> phiarr(nvertex);
    forAll(vertexInds, ivertex)
    {
        phiarr[ivertex] = phis[vertexInds[ivertex]];
    }

    scalar area = 0.0;
    for (int iseg = 0; iseg < nvertex; ++iseg)
    {
        const scalar& phiO = phiarr[iseg];
        const scalar& phiA = phiarr[(iseg+1)%nvertex];
        const vector& O = m_pp[vertexInds[iseg]];
        const vector& A = m_pp[vertexInds[(iseg+1)%nvertex]];
        area += std::fabs(0.5*Foam::mag((A-O) ^ (apex-O))) * calcLineFraction(phiO, phiA);
    }
    return area;
}

scalar GeometricTools::calcFaceAreaFraction(const Foam::face& vertexInds, CacheMap& phis, label faceInd)
{
    label nvertex = vertexInds.size();
    int sign_sum = 0;
    forAll(vertexInds, ivertex)
    {
        if (phis[vertexInds[ivertex]] > 0)
            ++sign_sum;
        else
            --sign_sum;
    }

    if (sign_sum == nvertex) // ALL_OUT
        return 0.0;
    if (sign_sum ==-nvertex) // ALL_IN
        return 1.0;

    return calcFaceArea(vertexInds, phis)/Foam::mag(m_fa[faceInd]);
}

}
