#ifndef KERNEL_H
#define KERNEL_H

#include "element.h"
#include "complex_t.h"

class Kernel
{
public:
    Kernel(){}
    virtual double direct(Element const & el1, Element const & el2) = 0;
};

#endif
