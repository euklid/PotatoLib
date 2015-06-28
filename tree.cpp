#include "tree.h"

void Tree::set_root(Cell *root)
{
    assert(m_dim == root->get_dimension());
    m_root = root;
}

Tree::Tree(unsigned int dim):
    m_dim(dim)
{}

Cell* Tree::get_root() {
    return m_root;
}

const Cell* Tree::get_root() const {
    return m_root;
}