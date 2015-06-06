#ifndef KERNEL_H
#define KERNEL_H

#include "element.h"

class Kernel
{
public:
    Kernel(){}
    double kernel(Element const & el1, Element const & el2) = 0;
};

#endif
