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
    int i, j, k;
    positionToIJK(point->position, &i, &j, &k);

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

        // verify that point has moved to new cell
        if (!(i==nexti && j == nextj && k == nextk)) {
            // remove grid point from old cell
            GridCell *oldCell = cellHashTable.getGridCell(i, j, k);
            oldCell->removeGridPoint(point);

            // remove cell from hash if empty
            if (oldCell->isEmpty()) {
                cellHashTable.removeGridCell(i, j, k);
                oldCell->reset();
                freeCells.push_back(oldCell);
            }

            // insert into new cell
            bool isCellInTable = cellHashTable.isGridCellInHash(nexti, nextj, nextk);
            if (isCellInTable) {
                GridCell *cell = cellHashTable.getGridCell(nexti, nextj, nextk);
                cell->insertGridPoint(point);
            } else {
                GridCell *cell = getNewGridCell();
                cell->initialize(nexti, nextj, nextk);
                cell->insertGridPoint(point);
                cellHashTable.insertGridCell(nexti, nextj, nextk, cell);
            }

            updateGridPointCellOffset(point, nexti, nextj, nextk);
        } else {
            //debug
            /*
            qDebug() << point->position.x <<
                        point->position.y <<
                        point->position.z <<
                        point->tx <<
                        point->ty <<
                        point->tz << i << j << k << nexti << nextj << nextk;
            positionToIJK(point->position, &nexti, &nextj, &nextk);
            qDebug() << nexti << nextj << nextk <<
                        newPos.x << newPos.y << newPos.z <<
                        trans.x << trans.y << trans.z;
            qDebug() << cellHashTable.isGridCellInHash(i, j, k);
            */

            // Current work around for a bug where point may not belong to the
            // grid cell at it's current position
            if (!cellHashTable.isGridCellInHash(i, j, k)) {
                GridCell *cell = getNewGridCell();
                cell->initialize(i, j, k);
                cell->insertGridPoint(point);
                cellHashTable.insertGridCell(nexti, nextj, nextk, cell);
                updateGridPointCellOffset(point, nexti, nextj, nextk);
            }
        }
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

    for (int ii=imin; ii<=imax; ii++) {
      for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=kmin; kk<=kmax; kk++) {

            if (cellHashTable.isGridCellInHash(ii, jj, kk)) {
                GridCell *cell = cellHashTable.getGridCell(ii, jj, kk);
                for (int idx=0; idx<(int)cell->points.size(); idx++) {
                    GridPoint *gp = cell->points[idx];
                    if (gp->id != ref) {
                        glm::vec3 v = p->position - gp->position;
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

std::vector<int> SpatialGrid::getIDsInRadiusOfPoint(int ref, double r) {
    std::vector<int> objects;

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

    GridPoint *gp;
    GridCell *cell;
    glm::vec3 v;
    for (int ii=imin; ii<=imax; ii++) {
      for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=kmin; kk<=kmax; kk++) {

            if (cellHashTable.isGridCellInHash(ii, jj, kk)) {
                cell = cellHashTable.getGridCell(ii, jj, kk);
                for (int idx=0; idx<(int)cell->points.size(); idx++) {
                    gp = cell->points[idx];
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

    /*
    timespec ts_beg, ts_end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_beg);

    float avg = 0;
    double r = 0.5;
    for (int i=0; i<(int)points.size(); i++) {
      GridPoint *p = points[i];
      std::vector<glm::vec3> objects = getObjectsInRadiusOfPoint(p->id, r);
      avg += objects.size();
    }
    avg /= points.size();

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);

    qDebug() << "average neighbours" << avg << (ts_end.tv_sec - ts_beg.tv_sec) +
                                (ts_end.tv_nsec - ts_beg.tv_nsec) / 1e9 << " sec";
    */
}










