#ifndef TREE_H
#define TREE_H

#include "cell.h"

class Tree_Iterator;

class Tree
{
public:
    Tree(unsigned int dim);
    virtual void set_root(Cell* root);
    virtual void build_tree(int max_elements) = 0;
    virtual Cell* get_root();
    virtual const Cell* get_root() const;
    virtual const std::vector<Cell*> & get_leaves() const;
    virtual std::vector<Cell*> & get_leaves();
    virtual const std::vector<Cell*> & get_cells() const;
    virtual ~Tree();
    virtual Tree_Iterator*  upward_iterator() = 0;
    virtual Tree_Iterator* downward_iterator() = 0;
protected:
    virtual Tree_Iterator* bfs_iterator() = 0;
    unsigned int m_dim;
    Cell* m_root;
    std::vector<Cell*> m_cells;
    std::vector<Cell*> m_leaves;
};

class Tree_Iterator
{
public:
    virtual Cell* next() = 0;
    virtual bool has_next() = 0;
};

#endif
