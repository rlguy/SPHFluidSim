#include "cellhash.h"

CellHash::CellHash()
{
    maxNumHashValues = 10000;
}

int CellHash::computeHash(int i, int j, int k) {
    return (541*i + 79*j + 31*k) % maxNumHashValues;
}

void CellHash::insertGridCell(int i, int j, int k, GridCell *cell) {
    int h = computeHash(i, j, k);

    if (cellMap.find(h) == cellMap.end()) {
        std::vector<GridCell*> newChain;
        std::pair<int,std::vector<GridCell*>> pair(h, newChain);
        cellMap.insert(pair);
    }

    cellMap[h].push_back(cell);
}

GridCell* CellHash::getGridCell(int i, int j, int k) {
    int h = computeHash(i, j, k);

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

bool CellHash::isGridCellInHash(int i, int j, int k) {
    int h = computeHash(i, j, k);

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

bool CellHash::isGridCellInHash(int i, int j, int k, GridCell *cell) {
    int h = computeHash(i, j, k);

    if (cellMap.find(h) == cellMap.end()) {
        return false;
    }

    GridCell *c;
    std::vector<GridCell*> chain = cellMap[h];
    for (int idx=0; idx<(int)chain.size(); idx++) {
        c = chain[idx];
        if (c->i == i && c->j == j && c->k == k) {
            cell = c;
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










