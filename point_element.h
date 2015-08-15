#ifndef POINT_ELEMENT_H
#define POINT_ELEMENT_H

#include "element.h"
#include "point.h"

class PointElement : public Element
{
public:
    
    PointElement(Point const & position,
                 unsigned int id,
                 int type = SOURCE) :
    Element(position.get_dimension(),id,type),
    m_pos(position)
    {}

    
    virtual const Point & get_position() const
    {
        return m_pos;
    }
    
private:
    Point m_pos;
};

#endif