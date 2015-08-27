#include "cell.h"
#include <sstream>

Cell::Cell(unsigned int dimension, double size, Point const & center) :
    m_dim(dimension),
    m_center(center),
    m_has_grid_pos(false),
    m_size(size),
    m_leaf_block_start_pos(-1),
    m_leaf_number(-1)
{
    assert(m_size > 0);
    m_grid_pos = Point(dimension);
}

Cell::~Cell() {}

const std::vector<Element *> &Cell::get_source_elements() const
{
    return m_src_elements;
}

void Cell::set_source_elements(const std::vector<Element *> &el)
{
#if DEBUG
    for(int i = 0; i< el.size(); i++)
    {
        assert(el[i]->get_type() & Element::SOURCE);
    }
#endif
    m_src_elements = el;
}

const std::vector<Element *> &Cell::get_target_elements() const
{
    return m_target_elements;
}

void Cell::set_target_elements(const std::vector<Element *> &el)
{
#if DEBUG
    for(int i = 0; i< el.size(); i++)
    {
        assert(el[i]->get_type() & Element::TARGET);
    }
#endif
    m_target_elements = el;
}

const std::vector<Element*> Cell::get_elements() const
{
    std::vector<Element*> elements(m_src_elements);
    elements.insert(elements.end(),m_target_elements.begin(),m_target_elements.end());
    return elements;
}

unsigned int Cell::number_of_elements() const
{
    return m_src_elements.size() + m_target_elements.size();
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

void Cell::init_empty_lists(unsigned int num_lists)
{
    m_lists = std::vector<std::map<unsigned int, Cell*> >(num_lists,std::map<unsigned int, Cell*>());
}

void Cell::add_to_list(Cell* cell, unsigned int list_number)
{
    assert(m_lists.size() > list_number);
    m_lists[list_number].insert(std::make_pair(cell->get_id(),cell));
}

std::vector<Cell*> Cell::get_list(unsigned int list_number) const
{
    assert(list_number < m_lists.size());
    std::vector<Cell*> linearized_list;
    std::map<unsigned int, Cell*> const & sel_list = m_lists.at(list_number);
    std::map<unsigned int, Cell*>::const_iterator it = sel_list.begin();
    for(; it != sel_list.end();  ++it)
    {
        linearized_list.push_back(it->second);
    }
    return linearized_list;
}

std::vector<unsigned int> Cell::get_list_ids(unsigned int list_number) const
{
    assert(list_number < m_lists.size());
    std::vector<unsigned int> linearized_list_ids;
    std::map<unsigned int, Cell*> const & sel_list = m_lists.at(list_number);
    std::map<unsigned int, Cell*>::const_iterator it = sel_list.begin();
    for(; it != sel_list.end();  ++it)
    {
        linearized_list_ids.push_back(it->first);
    }
    return linearized_list_ids;
}


//void Cell::add_to_interaction_list(Cell *cell)
//{
//    m_interaction_list.insert(std::make_pair(cell->get_id(),cell));
//}
//
//void Cell::add_to_direct_list(Cell *cell)
//{
//    m_direct_list.insert(std::make_pair(cell->get_id(),cell));
//}
//
//std::vector<Cell *> Cell::get_interaction_list() const
//{
//    std::vector<Cell*> linearized_interaction_list;
//    std::map<unsigned int, Cell*>::const_iterator it = m_interaction_list.begin();
//    for(; it != m_interaction_list.end();  ++it)
//    {
//        linearized_interaction_list.push_back(it->second);
//    }
//    return linearized_interaction_list;
//}
//
//std::vector<unsigned int> Cell::get_interaction_list_ids() const
//{
//    std::vector<unsigned int> linearized_interaction_list_ids;
//    std::map<unsigned int, Cell*>::const_iterator it = m_interaction_list.begin();
//    for(; it != m_interaction_list.end();  ++it)
//    {
//        linearized_interaction_list_ids.push_back(it->first);
//    }
//    return linearized_interaction_list_ids;
//}
//
//std::vector<Cell*> Cell::get_direct_list() const
//{
//    std::vector<Cell*> linearized_direct_list;
//    std::map<unsigned int, Cell*>::const_iterator it = m_direct_list.begin();
//    for(; it != m_direct_list.end();  ++it)
//    {
//        linearized_direct_list.push_back(it->second);
//    }
//    return linearized_direct_list;
//}
//std::vector<unsigned int> Cell::get_direct_list_ids() const
//{
//    std::vector<unsigned int> linearized_direct_list_ids;
//    std::map<unsigned int, Cell*>::const_iterator it = m_direct_list.begin();
//    for(; it != m_direct_list.end();  ++it)
//    {
//        linearized_direct_list_ids.push_back(it->first);
//    }
//    return linearized_direct_list_ids;
//}
//

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
    m_has_grid_pos = true;
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

void Cell::set_leaf_block_start_pos(int leaf_block_start_pos)
{
    m_leaf_block_start_pos = leaf_block_start_pos;
}

int Cell::get_leaf_block_start_pos() const
{
    return m_leaf_block_start_pos;
}

 int Cell::get_leaf_number() const
 {
     return m_leaf_number;
 }
 
 void Cell::set_leaf_number(int leaf_number)
 {
     m_leaf_number = leaf_number;
 }
 
std::string Cell::debug_info() const
{
    std::stringstream sstr;
    sstr << "Cell(id:"  <<m_id<<", lvl:"<<m_level << 
            ", pos:" << m_grid_pos << ((is_leaf())?", lf":", fa") << ")";
    return sstr.str();
}
