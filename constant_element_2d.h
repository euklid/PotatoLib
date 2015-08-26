/* 
 * File:   constant_element_2d.h
 * Author: euklid
 *
 * Created on August 24, 2015, 6:16 PM
 */

#ifndef CONSTANT_ELEMENT_2D_H
#define	CONSTANT_ELEMENT_2D_H

#include "element.h"
#include "point.h"

class ConstEl2D : public Element
{
public:

    ConstEl2D(Point const & start_node, Point const & end_node, unsigned int id, int type = SOURCE) :
    Element(2, id, type),
            m_start_node(start_node),
            m_end_node(end_node)
    {
        assert(start_node.get_dimension() == 2 && end_node.get_dimension() == 2);
        m_pos = (start_node + end_node)*0.5;
        m_len = Point::dist(start_node,end_node);
        assert(m_len > 0);
        m_el_norm = end_node - start_node;
        double tmp = m_el_norm[0];
        m_el_norm[0] = m_el_norm[1]/m_len;
        m_el_norm[1] = -tmp/m_len;
    }
    
    Point const & get_start_node() const
    {
        return m_start_node;
    }

    Point const & get_end_node() const
    {
        return m_end_node;
    }
    
    virtual const Point& get_position() const
    {
        return m_pos;
    }
    
    virtual double get_length() const
    {
        return m_len;
    }
    
    virtual Point const & get_el_norm() const
    {
        return m_el_norm;
    }
    
private:
    const Point m_start_node;
    const Point m_end_node;
    Point m_pos;
    double m_len;
    Point m_el_norm;

};

#endif	/* CONSTANT_ELEMENT_2D_H */

