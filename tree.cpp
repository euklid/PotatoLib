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

std::vector<Cell *> & Tree::get_leaves()
{
    return m_leaves;
}

const std::vector<Cell*> & Tree::get_cells() const
{
    return m_cells;
}

Point Tree::get_cell_grid_pos(Point const & cell_center,
                       unsigned int lvl,
                       Point const & root_cell_center,
                       double root_cell_size)
{
#if DEBUG
    assert(cell_center.get_dimension() == root_cell_center.get_dimension());
#endif
    Point center_shift(root_cell_center);
    
    for (int i = 0; i<center_shift.get_dimension(); i++) {
        center_shift[i] = root_cell_size/2 - root_cell_center[i];
    }
    
    Point shiftet_cell_center = cell_center + center_shift;
    double cell_size = root_cell_size/(1<<lvl);
    Point grid_coords(cell_center.get_dimension());
    
    for (int i = 0; i<root_cell_center.get_dimension(); i++) {
        grid_coords[i] = (int)(shiftet_cell_center[i]/cell_size);
    }
    return grid_coords;
}