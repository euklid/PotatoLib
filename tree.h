#ifndef TREE_H
#define TREE_H

#include "cell.h"

class Tree
{
public:
    Tree();
    Tree(unsigned int dim);
    virtual void set_root(Cell* root);
    virtual void build_tree(int max_elements) = 0;

    class Tree_Iterator {
    public:
        Tree_Iterator();
        Cell* next();
    };

private:
    unsigned int m_dim;
    Cell* m_root;
    std::vector<Cell*> m_cells;
};


#endif
