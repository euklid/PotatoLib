#ifndef TREE_H
#define TREE_H

#include "cell.h"

class Tree
{
public:
    Tree(unsigned int dim);
    virtual void set_root(Cell* root);
    virtual void build_tree(int max_elements) = 0;
    virtual Cell* get_root();
    virtual const Cell* get_root() const;

    class Tree_Iterator {
    public:
        Tree_Iterator(Tree* tree);
        virtual Cell* next() = 0;
    protected:
        Tree* m_tree;
    };

protected:
    unsigned int m_dim;
    Cell* m_root;
    std::vector<Cell*> m_cells;
};


#endif
