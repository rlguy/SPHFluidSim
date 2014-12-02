#include "spatialgrid.h"

SpatialGrid::SpatialGrid() {

}

SpatialGrid::SpatialGrid(double cell_size)
{
    size = cell_size;
    currentGridPointID = 0;
    numInitialFreeCells = 1000;

    initializeFreeCells();
}

void SpatialGrid::initializeFreeCells() {
    for (int i=0; i<numInitialFreeCells; i++) {
        GridCell *cell = new GridCell();
        freeCells.push_back(cell);
    }
}

int SpatialGrid::insertPoint(glm::vec3 p) {
    GridPoint *point = new GridPoint();
    point->position = p;
    point->id = generateUniqueGridPointID();
    points.push_back(point);

    std::pair<int,GridPoint*> pair(point->id, point);
    gridPointsByID.insert(pair);

    insertGridPointIntoGrid(point);
    return point->id;
}

void SpatialGrid::positionToIJK(glm::vec3 p, int *i, int *j, int *k) {
    double inv = 1 / size;
    *i = ceil(p.x*inv);
    *j = ceil(p.y*inv);
    *k = ceil(p.z*inv);

    double eps = 0.0000001;
    if (fmod(p.x, size) > eps) { *i = *i - 1; }
    if (fmod(p.y, size) > eps) { *j = *j - 1; }
    if (fmod(p.z, size) > eps) { *k = *k - 1; }
}

glm::vec3 SpatialGrid::IJKToPosition(int i, int j, int k) {
    return glm::vec3(i*size, j*size, k*size);
}

GridCell* SpatialGrid::getNewGridCell() {
    if (freeCells.size() == 0) {
        int n = 200;
        for (int i=0; i<n; i++) {
            freeCells.push_back(new GridCell());
        }
    }

    GridCell *cell = freeCells.back();
    freeCells.pop_back();
    return cell;
}

void SpatialGrid::insertGridPointIntoGrid(GridPoint *p) {
    int i, j, k;
    positionToIJK(p->position, &i, &j, &k);

    bool isCellInTable = cellHashTable.isGridCellInHash(i, j, k);
    if (isCellInTable) {
        GridCell *cell = cellHashTable.getGridCell(i, j, k);
        cell->insertGridPoint(p);
    } else {
        GridCell *cell = getNewGridCell();
        cell->initialize(i, j, k);
        cell->insertGridPoint(p);
        cellHashTable.insertGridCell(i, j, k, cell);
    }

}

int SpatialGrid::generateUniqueGridPointID() {
    int id = currentGridPointID;
    currentGridPointID++;
    return id;
}

void SpatialGrid::movePoint(int id, glm::vec3 newPos) {
    if (gridPointsByID.find(id) == gridPointsByID.end()) {
        return;
    }

    GridPoint *point = gridPointsByID[id];
    point->position = newPos;
}

void SpatialGrid::draw() {

    glColor3f(1.0, 0.0, 0.0);
    glPointSize(8.0);
    glBegin(GL_POINTS);
    for (int i=0; i<(int)points.size(); i++) {
      glm::vec3 p = points[i]->position;
      glVertex3f(p.x, p.y, p.z);
      //qDebug() << p.x << p.y << p.z << points[i]->id << points.size();
    }
    glEnd();

    std::vector<GridCell*> cells;
    cellHashTable.getGridCells(&cells);

    for (int i=0; i<(int)cells.size(); i++) {
        GridCell *c = cells[i];
        glm::vec3 pos = IJKToPosition(c->i, c->j, c->k);
        pos = pos + glm::vec3(0.5*size, 0.5*size, 0.5*size);
        utils::drawWireframeCube(pos, size);
    }

}
