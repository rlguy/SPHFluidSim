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

void GridCell::removeGridPoint(GridPoint *gp) {

    for (int i=0; i<(int)points.size(); i++) {
        if (points[i]->id == gp->id) {
            points.erase(points.begin() + i);
            return;
        }
    }
}

bool GridCell::isEmpty() {
    return points.size() == 0;
}

void GridCell::reset() {
    points.clear();
    i = 0; j = 0; k = 0;
}
