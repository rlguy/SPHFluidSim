#ifndef SPATIALGRID_H
#define SPATIALGRID_H

#include <vector>
#include <unordered_map>
#include <GL/glu.h>
#include <QDebug>
#include <cmath>
#include "glm/glm.hpp"
#include "cellhash.h"
#include "gridcell.h"
#include "utils.h"

struct GridPoint {
    glm::vec3 position;
    int id;
};

class SpatialGrid
{
public:
    SpatialGrid();
    SpatialGrid(double cell_size);
    int insertPoint(glm::vec3 point);
    void movePoint(int id, glm::vec3 position);
    void draw();

private:
    int generateUniqueGridPointID();
    void initializeFreeCells();
    void insertGridPointIntoGrid(GridPoint *p);
    void positionToIJK(glm::vec3 p, int *i, int *j, int *k);
    glm::vec3 IJKToPosition(int i, int j, int k);
    GridCell* getNewGridCell();

    double size;
    int currentGridPointID;
    std::vector<GridPoint*> points;
    std::unordered_map<int,GridPoint*> gridPointsByID;
    std::vector<GridCell*> freeCells;
    int numInitialFreeCells;
    CellHash cellHashTable;
};

#endif // SPATIALGRID_H













