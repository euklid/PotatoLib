#ifndef CELL_H
#define CELL_H

#include <vector>
#include "element.h"
#include "point.h"

class Cell
{
public:
    Cell();
    Cell(unsigned int dimension, unsigned int terms, double size, Point const & center);
    virtual void calc_moment() = 0;
    virtual const double get_moment();
    virtual std::vector<Cell*> divide() = 0;
    virtual void set_elements(std::vector<Element*> const & elements);
    virtual Cell* const get_father() const;
    virtual std::vector<Cell*> const & get_children() ;
    virtual unsigned int get_id() ;
    virtual unsigned int get_dimension() ;
    virtual Point const & get_center() ;
    virtual void set_center(Point const & center);
    virtual bool is_leaf() const;

protected:
    std::vector<Cell*> m_children;
    Cell* m_father;
    std::vector<Cell*> m_interaction_list;
    std::vector<Element*> m_elements;
    unsigned int m_id;
    unsigned int m_dim;
    unsigned int m_terms;
    Point m_center;
    double m_moment;
    bool has_moment;
    double m_size;
};

#endif
