#include "ugrid.h"
#include <iostream>
namespace sdfibm{

UGrid::UGrid(BBox& bbox, double delta)
{
   m_bbox  = bbox;
   m_delta = delta;
   m_deltaINV = 1.0/m_delta;

   m_nx = std::ceil((bbox.high[0] - bbox.low[0])*m_deltaINV);
   m_ny = std::ceil((bbox.high[1] - bbox.low[1])*m_deltaINV);
   m_nz = std::ceil((bbox.high[2] - bbox.low[2])*m_deltaINV);

   m_nynz = m_ny*m_nz;
   m_num_cells = m_nx*m_ny*m_nz;

   std::vector<int> slot;
   slot.reserve(10); // by default each cell has room for 10 objects
   for(int i = 0; i < m_num_cells; ++i)
       m_map[i] = slot;
}

void UGrid::report(std::ostream& os, bool detail)
{
    os << "[UGrid]\n";
    m_bbox.report(os);
    os << "[UGrid] Cells: " << m_num_cells << std::endl;
    if (detail)
    {
        for(int i = 0; i < m_nx; ++i) {
            for(int j = 0; j < m_ny; ++j) {
                for(int k = 0; k < m_nz; ++k) {
                    std::vector<int>& lst = getObjectList(i, j, k);
                    if(!lst.empty())
                    {
                        char msg[80];
                        sprintf(msg, "\tC[%d,%d,%d](%ld):", i,j,k,lst.size());
                        os << msg;
                        for(int pid : lst)
                        {
                            os << pid << ' ';
                        }
                        os << std::endl;
                    }
        } } } // end of triple loop
    }
}

void UGrid::generateCollisionPairs(std::vector<CollisionPair>& collision_pairs)
{
    for(int i = 0; i < m_nx; ++i) {
        for(int j = 0; j < m_ny; ++j) {
            for(int k = 0; k < m_nz; ++k) {
               int myid = hash(i, j, k);
               if(m_map[myid].empty())
                   continue;
               // note: if use unsigned nbi, -1 is implicitly casted into 0
               for(int nbi = i-1; nbi <= i+1; ++nbi) {
                   for(int nbj = j-1; nbj <= j+1; ++nbj) {
                       for(int nbk = k-1; nbk <= k+1; ++nbk) {
                           if(nbi < 0 || nbi > m_nx-1) continue;
                           if(nbj < 0 || nbj > m_ny-1) continue;
                           if(nbk < 0 || nbk > m_nz-1) continue;

                           int nbid = hash(nbi, nbj, nbk);
                           if(m_map[nbid].empty()) continue;

                           std::vector<int>& lst = m_map[myid];
                           std::vector<int>& nblist = m_map[nbid];
                           for(int pi : lst)
                               for(int qi : nblist)
                                   if(pi < qi)
                                       collision_pairs.push_back(CollisionPair(pi, qi));
                } } } // end of inner triple loop
    } } } // end of outer triple loop
}

}
