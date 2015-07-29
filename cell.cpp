#include "cell.h"

Cell::Cell(unsigned int dimension, unsigned int terms, double size, Point const & center) :
    m_dim(dimension),
    m_terms(terms),
    m_size(size),
    m_center(center),
    m_has_moment(false),
    m_has_grid_pos(false)
{
    assert(m_size > 0);
    m_grid_pos = Point(dimension);
}

Cell::~Cell() {}

void Cell::set_elements(std::vector<Element*> const & elements)
{
    m_elements = elements;
}

std::vector<Element*> & Cell::get_elements()
{
    return m_elements;
}

const std::vector<Element*> & Cell::get_elements() const
{
    return m_elements;
}

unsigned int Cell::number_of_elements() const
{
    return m_elements.size();
}

Cell * const Cell::get_father() const
{
    return m_father;
}

void Cell::set_father(Cell * father)
{
    m_father = father;
}

std::vector<Cell*> const & Cell::get_children() const
{
    return m_children;
}

unsigned int Cell::get_id() const
{
    return m_id;
}

void Cell::set_id(unsigned int index)
{
    m_id = index;
}

unsigned int Cell::get_dimension() const
{
    return m_dim;
}

Point const & Cell::get_center() const
{
    return m_center;
}

void Cell::set_center(Point const & center)
{
    m_center = center;
}

const std::vector<double> & Cell::get_moments() const
{
    return m_moments;
}

const std::vector<complex_t>& Cell::get_moments_cmp() const
{
    return m_moments_cmp;
}

bool Cell::is_leaf() const
{
    return m_children.empty();
}

void Cell::set_level(unsigned int lvl)
{
    m_level = lvl;
}

unsigned int Cell::get_level() const
{
    return m_level;
}

const double Cell::get_size() const
{
    return m_size;
}

void Cell::set_moments(std::vector<double> const & moments)
{
    m_moments = moments;
}

void Cell::set_moments_cmp(std::vector<complex_t> const & moments_cmp)
{
    m_moments_cmp = moments_cmp;
}

void Cell::add_to_interaction_list(Cell *cell)
{
    m_interaction_list.insert(std::make_pair(cell->get_id(),cell));
}

void Cell::add_to_direct_list(Cell *cell)
{
    m_direct_list.insert(std::make_pair(cell->get_id(),cell));
}

std::vector<Cell *> Cell::get_interaction_list() const
{
    std::vector<Cell*> linearized_interaction_list;
    std::map<unsigned int, Cell*>::const_iterator it = m_interaction_list.begin();
    for(; it != m_interaction_list.end();  ++it)
    {
        linearized_interaction_list.push_back(it->second);
    }
    return linearized_interaction_list;
}

std::vector<unsigned int> Cell::get_interaction_list_ids() const
{
    std::vector<unsigned int> linearized_interaction_list_ids;
    std::map<unsigned int, Cell*>::const_iterator it = m_interaction_list.begin();
    for(; it != m_interaction_list.end();  ++it)
    {
        linearized_interaction_list_ids.push_back(it->first);
    }
    return linearized_interaction_list_ids;
}

bool Cell::has_level_grid_position() const
{
    return m_has_grid_pos;
}

const Point &Cell::get_level_grid_position() const
{
    return m_grid_pos;
}

void Cell::set_level_grid_position(const Point &grid_pos)
{
    assert(grid_pos.get_dimension() == m_dim);
    m_grid_pos = grid_pos;
}

const std::vector<double> & Cell::get_local_exps() const
{
    return m_local_exps;
}

const std::vector<complex_t> & Cell::get_local_exps_cmp() const
{
    return m_local_exps_cmp;
}

void Cell::set_local_exps(std::vector<double> const & local_exps)
{
    m_local_exps = local_exps;
}

void Cell::set_local_exps_cmp(std::vector<complex_t> const & local_exps_cmp)
{
    m_local_exps_cmp = local_exps_cmp;
}
