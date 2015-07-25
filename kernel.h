#ifndef KERNEL_H
#define KERNEL_H

#include "element.h"
#include "complex_t.h"
#include "cell.h"

class Kernel
{
public:
    Kernel(){}
    virtual double direct(Element const & el1, Element const & el2) = 0;
    virtual complex_t direct_cmp(Element const & el1, Element const & el2) = 0;
    virtual void M2M(Cell const & first, Cell & second) = 0;
    virtual void M2L(Cell const & first, Cell & second) = 0;
    virtual void L2L(Cell const & first, Cell & second) = 0;

};

#endif
