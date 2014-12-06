#include "cellhash.h"

CellHash::CellHash()
{
    maxNumHashValues = 10000;
}

inline long CellHash::computeHash(int i, int j, int k) {
    return (abs(541*(long)i + 79*(long)j + 31*(long)k) % maxNumHashValues);
}

void CellHash::insertGridCell(GridCell *cell) {
    long h = computeHash(cell->i, cell->j, cell->k);

    if (cellMap.find(h) == cellMap.end()) {
        std::vector<GridCell*> newChain;
        std::pair<long,std::vector<GridCell*>> pair(h, newChain);
        cellMap.insert(pair);
    }

    cellMap[h].push_back(cell);
}

void CellHash::removeGridCell(GridCell *cell) {
    int i = cell->i;
    int j = cell->j;
    int k = cell->k;
    long h = computeHash(i, j, k);

    if (cellMap.find(h) == cellMap.end()) {
        qDebug() << "cant find cell" << i << j << k << h;
        return;
    }

    // remove from hash chain
    bool isRemoved = false;
    std::vector<GridCell*> chain = cellMap[h];
    for (int idx=0; idx<(int)chain.size(); idx++) {
        GridCell *c = (cellMap[h])[idx];
        if (c->i == i && c->j == j && c->k == k) {
            cellMap[h].erase(cellMap[h].begin() + idx);
            isRemoved = true;
            break;
        }
    }

    if (!isRemoved) {
        qDebug() << "Could not find/remove gridcell" << i << j << k;
    }

    // remove chain from map if empty
    if (chain.size() == 0) {
        cellMap.erase(h);
    }

}

GridCell* CellHash::getGridCell(int i, int j, int k) {
    long h = computeHash(i, j, k);

    GridCell *c;
    std::vector<GridCell*> chain = cellMap[h];
    for (int idx=0; idx<(int)chain.size(); idx++) {
        c = chain[idx];
        if (c->i == i && c->j == j and c->k == k) {
            return c;
        }
    }

    return c;
}

GridCell* CellHash::findGridCell(int i, int j, int k, bool *isGridCellFound) {
    long h = computeHash(i, j, k);

    GridCell *c;
    std::vector<GridCell*> chain = cellMap[h];
    for (int idx=0; idx<(int)chain.size(); idx++) {
        c = chain[idx];
        if (c->i == i && c->j == j and c->k == k) {
            *isGridCellFound = true;
            return c;
        }
    }

    *isGridCellFound = false;
    return c;
}

bool CellHash::isGridCellInHash(int i, int j, int k) {
    long h = computeHash(i, j, k);

    if (cellMap.find(h) == cellMap.end()) {
        return false;
    }

    GridCell *c;
    std::vector<GridCell*> chain = cellMap[h];
    for (int idx=0; idx<(int)chain.size(); idx++) {
        c = chain[idx];
        if (c->i == i && c->j == j and c->k == k) {
            return true;
        }
    }

    return false;
}



void CellHash::getGridCells(std::vector<GridCell*> *cells) {
    for (std::pair<int, std::vector<GridCell*>> pair: cellMap) {
        for (int i=0; i < (int)pair.second.size(); i++) {
            cells->push_back(pair.second[i]);
        }
    }

}










