#ifndef BOUNDINGBOX_HPP
#define BOUNDINGBOX_HPP

#include <iostream>

namespace sdfibm {
struct BBox
{
    double low[3], high[3];

    BBox(){};
    BBox(const double* low_in, const double* high_in)
    {
        for(int i = 0; i < 3; ++i) {
            low[i] = low_in[i];
            high[i] = high_in[i];
        }
    }

    void report(std::ostream& os = std::cout) const
    {
        char msg[200];
        sprintf(msg, "[BBox] [%g,%g,%g]x[%g,%g,%g]\n", low[0], low[1], low[2], high[0], high[1], high[2]);
        os << msg;
    }
};

}
#endif // BOUNDINGBOX_HPP
