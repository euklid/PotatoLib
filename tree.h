#ifndef TREE_H
#define TREE_H

#include "cell.h"

class Tree_Iterator;

class Tree
{
public:
    
    static Point get_cell_grid_pos(Point const & cell_center,
                                   unsigned int lvl,
                                   Point const & root_cell_center,
                                   double root_cell_size);
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
    
    //for tree code based representation
    std::vector<std::vector<unsigned int> > m_lvl_ids;
    unsigned int m_max_level;
};

class Tree_Iterator
{
public:
    virtual Cell* next() = 0;
    virtual bool has_next() = 0;
};

#endif
