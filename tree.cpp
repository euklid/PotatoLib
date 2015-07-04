#include "tree.h"

void Tree::set_root(Cell *root)
{
    assert(m_dim == root->get_dimension());
    m_root = root;
}

Tree::Tree(unsigned int dim):
    m_dim(dim)
{}

Tree::~Tree() {}

Cell* Tree::get_root()
{
    return m_root;
}

const Cell* Tree::get_root() const
{
    return m_root;
}

const std::vector<Cell*> & Tree::get_leaves() const
{
    return m_leaves;
}

std::vector<Cell*> Tree::get_leaves()
{
    return m_leaves;
}

const std::vector<Cell*> & Tree::get_cells() const
{
    return m_cells;
}