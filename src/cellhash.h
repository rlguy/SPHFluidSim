#ifndef CELLHASH_H
#define CELLHASH_H

#include <QDebug>
#include <unordered_map>
#include <vector>
#include "gridcell.h"

class CellHash
{
public:
    CellHash();
    bool isGridCellInHash(int i, int j, int k);
    bool isGridCellInHash(int i, int j, int k, GridCell *cell);
    void insertGridCell(int i, int j, int k, GridCell *cell);
    GridCell* getGridCell(int i, int j, int k);
    void getGridCells(std::vector<GridCell*> *cells);

private:
    int computeHash(int i, int j, int k);

    int maxNumHashValues;
    std::unordered_map<int, std::vector<GridCell*>> cellMap;
};

#endif // CELLHASH_H
