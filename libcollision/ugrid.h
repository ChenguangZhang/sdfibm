#ifndef UGRID_HPP
#define UGRID_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>

#include "bbox.h"

typedef std::pair<int, int> CollisionPair;

class UGrid
{
private:
    BBox m_bbox;
    int m_nx, m_ny, m_nz, m_nynz, m_num_cells;
    double m_delta;
    double m_deltaINV;
    std::unordered_map<int, std::vector<int>> m_map;

public:
    UGrid(BBox& bbox, double delta);

    void generateCollisionPairs(std::vector<CollisionPair> &collision_pairs);
    void report();

    inline void insert(double x, double y, double z, int id)
    {
        m_map[hash(x, y, z)].push_back(id);
    }

    std::vector<int>& getObjectList(int i, int j, int k)
    {
        return m_map[hash(i,j,k)];
    }

    void clear()
    {
        for(auto& item : m_map)
            item.second.clear();
    }

private:
    inline int hash(int i, int j, int k)
    {
        return i*m_nynz + j*m_nz + k;
    }

    inline int hash(double x, double y, double z)
    {
        int i = std::floor((x-m_bbox.low[0])*m_deltaINV);
        int j = std::floor((y-m_bbox.low[1])*m_deltaINV);
        int k = std::floor((z-m_bbox.low[2])*m_deltaINV);
        return hash(i, j, k);
    }
};

#endif
