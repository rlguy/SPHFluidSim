#include "spatialgrid.h"

SpatialGrid::SpatialGrid() {

}

SpatialGrid::SpatialGrid(double cell_size)
{
    size = cell_size;
    currentGridPointID = 0;
    numInitialFreeCells = 10000;

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
    point->isMarkedForRemoval = false;
    points.push_back(point);

    std::pair<int,GridPoint*> pair(point->id, point);
    gridPointsByID.insert(pair);

    insertGridPointIntoGrid(point);
    return point->id;
}

void SpatialGrid::positionToIJK(glm::vec3 p, int *i, int *j, int *k) {
    double inv = 1 / size;
    *i = ceil(p.x*inv)-1;
    *j = ceil(p.y*inv)-1;
    *k = ceil(p.z*inv)-1;

    double eps = 0.00000000001;
    if (fabs(fmod(p.x, size)) < eps) {
        *i = *i + 1;
    }
    if (fabs(fmod(p.y, size)) < eps) {
        *j = *j + 1;
    }
    if (fabs(fmod(p.z, size)) < eps) {
        *k = *k + 1;
    }
}

glm::vec3 SpatialGrid::IJKToPosition(int i, int j, int k) {
    return glm::vec3(i*size, j*size, k*size);
}

GridCell* SpatialGrid::getNewGridCell(int i, int j, int k) {
    if (freeCells.size() == 0) {
        int n = 200;
        for (int i=0; i<n; i++) {
            freeCells.push_back(new GridCell());
        }
    }

    GridCell *cell = freeCells.back();
    freeCells.pop_back();

    cell->initialize(i, j, k);
    return cell;
}

void SpatialGrid::insertGridPointIntoGrid(GridPoint *p) {
    int i, j, k;
    positionToIJK(p->position, &i, &j, &k);

    bool isCellInTable = false;
    GridCell *cell = cellHashTable.findGridCell(i, j, k, &isCellInTable);
    if (isCellInTable) {
        cell->insertGridPoint(p);
    } else {
        cell = getNewGridCell(i, j ,k);
        cell->insertGridPoint(p);
        cellHashTable.insertGridCell(cell);
    }

    // offset used for updating position
    updateGridPointCellOffset(p, i, j, k);
}

int SpatialGrid::generateUniqueGridPointID() {
    int id = currentGridPointID;
    currentGridPointID++;
    return id;
}

void SpatialGrid::updateGridPointCellOffset(GridPoint *gp, int i, int j, int k) {
    glm::vec3 cp = IJKToPosition(i, j, k);
    gp->tx = gp->position.x - cp.x;
    gp->ty = gp->position.y - cp.y;
    gp->tz = gp->position.z - cp.z;
}

void SpatialGrid::movePoint(int id, glm::vec3 newPos) {
    if (gridPointsByID.find(id) == gridPointsByID.end()) {
        return;
    }

    GridPoint *point = gridPointsByID[id];
    int i = point->i;
    int j = point->j;
    int k = point->k;

    glm::vec3 trans = newPos - point->position;
    point->tx += trans.x;
    point->ty += trans.y;
    point->tz += trans.z;
    point->position = newPos;

    // point has moved to new cell
    if (point->tx >= size || point->ty >= size || point->tz >= size ||
              point->tx < 0 || point->ty < 0 || point->tz < 0) {
        int nexti, nextj, nextk;
        positionToIJK(point->position, &nexti, &nextj, &nextk);

        // remove grid point from old cell
        GridCell *oldCell = cellHashTable.getGridCell(i, j, k);
        oldCell->removeGridPoint(point);

        // remove cell from hash if empty
        if (oldCell->isEmpty()) {
            cellHashTable.removeGridCell(oldCell);
            oldCell->reset();
            freeCells.push_back(oldCell);
        }

        // insert into new cell
        bool isCellInTable = false;
        GridCell *cell = cellHashTable.findGridCell(nexti, nextj, nextk, &isCellInTable);
        if (isCellInTable) {
            cell->insertGridPoint(point);
        } else {
            GridCell *cell = getNewGridCell(nexti, nextj, nextk);
            cell->insertGridPoint(point);
            cellHashTable.insertGridCell(cell);
        }

        updateGridPointCellOffset(point, nexti, nextj, nextk);
    }
}

void SpatialGrid::removePoint(int id) {
    if (gridPointsByID.find(id) == gridPointsByID.end()) {
        return;
    }
    isCellRemoved = true;

    GridPoint *point = gridPointsByID[id];
    gridPointsByID.erase(id);

    int i = point->i;
    int j = point->j;
    int k = point->k;

    point->isMarkedForRemoval = true;
    bool isCellInTable = false;
    GridCell *cell = cellHashTable.findGridCell(i, j, k, &isCellInTable);

    if (!isCellInTable) {
        return;
    }

    cell->removeGridPoint(point);
    if (cell->isEmpty()) {
        cellHashTable.removeGridCell(cell);
        cell->reset();
        freeCells.push_back(cell);
    }
}

std::vector<glm::vec3> SpatialGrid::getObjectsInRadiusOfPoint(int ref, double r) {
    std::vector<glm::vec3> objects;

    if (gridPointsByID.find(ref) == gridPointsByID.end()) {
        return objects;
    }

    GridPoint *p = gridPointsByID[ref];
    double tx = p->tx;
    double ty = p->ty;
    double tz = p->tz;
    int i, j, k;
    positionToIJK(p->position, &i, &j, &k);
    double inv = 1/size;
    double rsq = r*r;

    int imin = i - fmax(0, ceil((r-tx)*inv));
    int jmin = j - fmax(0, ceil((r-ty)*inv));
    int kmin = k - fmax(0, ceil((r-tz)*inv));
    int imax = i + fmax(0, ceil((r-size+tx)*inv));
    int jmax = j + fmax(0, ceil((r-size+ty)*inv));
    int kmax = k + fmax(0, ceil((r-size+tz)*inv));

    GridCell *cell;
    GridPoint *gp;
    glm::vec3 v;
    std::vector<GridPoint*> points;
    for (int ii=imin; ii<=imax; ii++) {
      for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=kmin; kk<=kmax; kk++) {

            bool isInHash = false;
            cell = cellHashTable.findGridCell(ii, jj, kk, &isInHash);
            if (isInHash) {
                points = cell->getGridPoints();
                for (int idx=0; idx<(int)points.size(); idx++) {
                    gp = points[idx];
                    if (gp->id != ref) {
                        v = p->position - gp->position;
                        if (glm::dot(v, v) < rsq) {
                            objects.push_back(gp->position);
                        }
                    }
                }
            }

        }
      }
    }

    return objects;
}

std::vector<int> SpatialGrid::fastIDNeighbourSearch(int ref, double r, GridPoint *p) {
    std::vector<int> objects;

    bool isInHash = false;
    GridCell *cell = cellHashTable.findGridCell(p->i, p->j, p->k, &isInHash);
    if (!isInHash) {
        return objects;
    }

    std::vector<GridPoint*> points = cell->getGridPoints();
    GridPoint *gp;
    glm::vec3 v;
    double rsq = r*r;
    for (uint i=0; i<points.size(); i++) {
        gp = points[i];
        if (gp->id != ref) {
            v = p->position - gp->position;
            if (glm::dot(v, v) < rsq) {
                objects.push_back(gp->id);
            }
        }
    }

    std::vector<GridCell*> neighbours = cell->neighbours;
    for (uint i=0; i<neighbours.size(); i++) {
        points = neighbours[i]->getGridPoints();
        for (uint j=0; j<points.size(); j++) {
            gp = points[j];
            v = p->position - gp->position;
            if (glm::dot(v, v) < rsq) {
                objects.push_back(gp->id);
            }
        }
    }

    return objects;
}

std::vector<int> SpatialGrid::getIDsInRadiusOfPoint(int ref, double r) {

    if (gridPointsByID.find(ref) == gridPointsByID.end()) {
        std::vector<int> objects;
        return objects;
    }

    GridPoint *p = gridPointsByID[ref];
    double tx = p->tx;
    double ty = p->ty;
    double tz = p->tz;
    int i, j, k;
    positionToIJK(p->position, &i, &j, &k);
    double inv = 1/size;
    double rsq = r*r;

    int imin = i - fmax(0, ceil((r-tx)*inv));
    int jmin = j - fmax(0, ceil((r-ty)*inv));
    int kmin = k - fmax(0, ceil((r-tz)*inv));
    int imax = i + fmax(0, ceil((r-size+tx)*inv));
    int jmax = j + fmax(0, ceil((r-size+ty)*inv));
    int kmax = k + fmax(0, ceil((r-size+tz)*inv));

    if (imax - imin <= 3 and imax - imin >= 1) {
        return fastIDNeighbourSearch(ref, r, p);
    }

    std::vector<int> objects;
    GridPoint *gp;
    GridCell *cell;
    glm::vec3 v;
    std::vector<GridPoint*> points;
    for (int ii=imin; ii<=imax; ii++) {
      for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=kmin; kk<=kmax; kk++) {

            bool isInHash = false;
            cell = cellHashTable.findGridCell(ii, jj, kk, &isInHash);
            if (isInHash) {
                points = cell->getGridPoints();
                for (int idx=0; idx<(int)points.size(); idx++) {
                    gp = points[idx];
                    if (gp->id != ref) {
                        v = p->position - gp->position;
                        if (glm::dot(v, v) < rsq) {
                            objects.push_back(gp->id);
                        }
                    }
                }
            }

        }
      }
    }

    return objects;
}

void SpatialGrid::removeGridPointsMarkedForRemoval() {
    if (points.size() == 0 || !isCellRemoved) {
        return;
    }

    GridPoint *p;
    for (int i=(int)points.size()-1; i>=0; i--) {
        if (points[i]->isMarkedForRemoval) {
            p = points[i];
            points.erase(points.begin() + i);

            delete p;
        }
    }

    isCellRemoved = false;
}

void SpatialGrid::update() {
    removeGridPointsMarkedForRemoval();

    // update each cell's cell neighbours
    std::vector<GridCell*> cells;
    cellHashTable.getGridCells(&cells);

    GridCell* cell;
    GridCell* gc;
    for (uint idx=0; idx<cells.size(); idx++) {
        cell = cells[idx];
        cell->neighbours.clear();

        int ii = cell->i;
        int jj = cell->j;
        int kk = cell->k;

        for (int k=kk-1; k<=kk+1; k++) {
            for (int j=jj-1; j<=jj+1; j++) {
                for (int i=ii-1; i<=ii+1; i++) {
                    if (!(i==ii && j==jj && k==kk)) {
                        bool isInTable = false;
                        gc = cellHashTable.findGridCell(i, j, k, &isInTable);
                        if (isInTable) {
                            cell->neighbours.push_back(gc);
                        }
                    }
                }
            }
        }

    }
}

void SpatialGrid::draw() {
    if (points.size() == 0) { return; }

    glColor3f(1.0, 0.4, 0.0);
    glPointSize(6.0);
    glBegin(GL_POINTS);
    for (int i=0; i<(int)points.size(); i++) {
      glm::vec3 p = points[i]->position;
      glVertex3f(p.x, p.y, p.z);
    }
    glEnd();

    std::vector<GridCell*> cells;
    cellHashTable.getGridCells(&cells);

    glLineWidth(1.0);
    glColor4f(0.0, 0.0, 1.0, 0.4);
    for (int i=0; i<(int)cells.size(); i++) {
        GridCell *c = cells[i];
        glm::vec3 pos = IJKToPosition(c->i, c->j, c->k);
        pos = pos + glm::vec3(0.5*size, 0.5*size, 0.5*size);
        utils::drawWireframeCube(pos, size);
    }
}










