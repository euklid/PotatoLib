/* 
 * File:   tree2d.cpp
 * Author: euklid
 * 
 * Created on June 13, 2015, 2:55 PM
 */

#include "tree2d.h"
#include <queue>

Tree2D::Tree2D() : Tree(2) {
}

Tree2D::~Tree2D() {
}

void Tree2D::build_tree(int max_elements) {
    std::queue<Cell*> undivided_cells;
    undivided_cells.push(m_root);
    while (!undivided_cells.empty()) 
    {
        Cell* cur_cell = undivided_cells.front();
        undivided_cells.pop();
        std::vector<Cell*> new_cells = cur_cell->divide();
        for (int i = 0; i<new_cells.size(); i++)
        {
            undivided_cells.push(new_cells[i]);
        }
    }
}
