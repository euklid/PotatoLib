#ifndef ELEMENT_H
#define ELEMENT_H

#include "point.h"

class Element
{
public:
    Element(unsigned int dimension) :
        m_dim(dimension)
    {}
    
    virtual double get_value() const
    {
        return m_value;
    }

    virtual Point get_position() const = 0;
    
    virtual bool operator==(Element const & el) = 0;

protected:
    double m_value;
    unsigned int m_dim;
};

#endif
