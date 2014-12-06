#ifndef GRIDCELL_H
#define GRIDCELL_H

#include <vector>
#include <QDebug>
#include <cmath>
#include "glm/glm.hpp"
#include "gridpoint.h"

class GridCell
{
public:
    GridCell();
    void reset();
    void initialize(int i, int j, int k);
    void insertGridPoint(GridPoint *gp);
    void removeGridPoint(GridPoint *gp);
    bool isEmpty();
    std::vector<GridPoint*> getGridPoints();

    std::vector<GridCell*> neighbours;
    int i, j, k;

private:
    std::vector<GridPoint*> points;
};

#endif // GRIDCELL_H
