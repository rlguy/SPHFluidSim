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
    void insertGridCell(int i, int j, int k, GridCell *cell);
    void removeGridCell(int i, int j, int k);
    GridCell* getGridCell(int i, int j, int k);
    void getGridCells(std::vector<GridCell*> *cells);

private:
    inline long computeHash(int i, int j, int k);

    long maxNumHashValues;
    std::unordered_map<long, std::vector<GridCell*>> cellMap;
};

#endif // CELLHASH_H
