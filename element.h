#ifndef ELEMENT_H
#define ELEMENT_H

#include "point.h"

class Element
{
public:
    Element(unsigned int dimension) :
        m_dim(dimension)
    {}
    virtual double get_value()
    {
        return m_value;
    }

    virtual Point get_position() = 0;

protected:
    double m_value;
    unsigned int m_dim;
};

#endif
