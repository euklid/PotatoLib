#include "tree.h"

void Tree::set_root(Cell *root)
{
    assert(m_dim == root->get_dimension());
}
