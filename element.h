#ifndef ELEMENT_H
#define ELEMENT_H

#include "point.h"
#include "complex_t.h"

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
    
    virtual complex_t get_cmp_value() const
    {
        return m_cmp_value;
    }

    virtual Point get_position() const = 0;
    
    virtual bool operator==(Element const & el) const = 0 ;

protected:
    double m_value;
    unsigned int m_dim;
    complex_t m_cmp_value;
};

#endif
