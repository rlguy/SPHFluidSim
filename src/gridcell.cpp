#include "gridcell.h"

GridCell::GridCell()
{
    i = 0; j = 0; k = 0;
}

void GridCell::initialize(int ii, int jj, int kk) {
    i = ii; j = jj; k = kk;
}

void GridCell::insertGridPoint(GridPoint *gp) {
    points.push_back(gp);
}

void GridCell::reset() {
    while (!points.empty())
    {
        points.pop_back();
    }
    i = 0; j = 0; k = 0;
}
