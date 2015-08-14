#ifndef ELEMENT_H
#define ELEMENT_H

#include "point.h"

class Element
{
public:

    typedef enum
    {
        SOURCE = 1 << 0,
        TARGET = 1 << 1
    } element_type;

    Element(unsigned int dimension, unsigned int id, int type = SOURCE) :
        m_dim(dimension),
        m_type(type),
        m_id(id)
    {}
    
    virtual double get_value() const
    {
        return m_value;
    }
    
    virtual element_type get_type() const
    {
        return (element_type)m_type;
    }

    virtual Point get_position() const = 0;
    
    virtual bool operator==(Element const & el) const = 0 ;

    virtual void set_value(double val)
    {
        m_value = val;
    }

    virtual void set_target_value(double val)
    {
        m_target_value = val;
    }

    virtual double get_target_value() const
    {
        return m_target_value;
    }

protected:
    double m_value;
    double m_target_value;
    unsigned int m_dim;
    unsigned int m_id;
    int m_type;
};

#endif
