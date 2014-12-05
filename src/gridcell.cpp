#include "gridcell.h"

GridCell::GridCell()
{
    i = 0; j = 0; k = 0;
}

void GridCell::initialize(int ii, int jj, int kk) {
    i = ii; j = jj; k = kk;
}

void GridCell::insertGridPoint(GridPoint *gp) {
    gp->i = i; gp->j = j; gp->k = k;
    gp->isInGridCell = true;
    points.push_back(gp);
}

void GridCell::removeGridPoint(GridPoint *gp) {

    for (int i=0; i<(int)points.size(); i++) {
        if (points[i]->id == gp->id) {
            gp->isInGridCell = false;
            points.erase(points.begin() + i);
            return;
        }
    }

    qDebug() << "Cannot find gridpoint" << gp->i << gp->j << gp->k << gp->id;
}

std::vector<GridPoint*> GridCell::getGridPoints() {
    return points;
}

bool GridCell::isEmpty() {
    return points.size() == 0;
}

void GridCell::reset() {
    for (uint i=0; i<points.size(); i++) {
        points[i]->isInGridCell = false;
    }
    points.clear();
    i = 0; j = 0; k = 0;
}
