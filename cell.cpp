#include "cell.h"

Cell::Cell(unsigned int dimension, unsigned int terms) :
    m_dim(dimension),
    m_terms(terms),
    has_moment(false)
{
    m_center = Point(dimension);
}

void Cell::set_elements(std::vector<Element*> const & elements)
{
    m_elements = std::vector<Element*>(elements);
}

Cell* const Cell::get_father() const
{
    return m_father;
}

std::vector<Cell*> const & Cell::get_children()
{
    return m_children;
}

unsigned int Cell::get_id()
{
    return m_id;
}

unsigned int Cell::get_dimension()
{
    return m_dim;
}

Point const & Cell::get_center()
{
    return m_center;
}

void Cell::set_center(Point const & center)
{
    m_center = center;
}

const double Cell::get_moment()
{
    if(!has_moment)
    {
        calc_moment();
    }

    return m_moment;
}

bool Cell::is_leaf() const {
    return m_children.empty();
}
