#ifndef GRIDCELL_H
#define GRIDCELL_H

#include <vector>
#include "glm/glm.hpp"

struct GridPoint;

class GridCell
{
public:
    GridCell();
    void reset();
    void initialize(int i, int j, int k);
    void insertGridPoint(GridPoint *gp);

    int i, j, k;

private:
    std::vector<GridPoint*> points;
};

#endif // GRIDCELL_H
