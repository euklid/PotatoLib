#ifndef KERNEL_H
#define KERNEL_H

#include "element.h"
#include "complex_t.h"
#include "cell.h"
#include <vector>

class Kernel
{
public:
    Kernel(){}
    virtual double direct(Element const & el1, Element const & el2) const = 0;
    virtual complex_t direct_cmp(Element const & el1, Element const & el2) const  = 0;
    virtual std::vector<complex_t> calc_moments_cmp(Cell const & cell, int num_moments) const = 0;
    virtual std::vector<double> calc_moments(Cell const & cell, int num_moments) const = 0;
    virtual void M2M(Cell const & first, Cell & second) const = 0;
    virtual void M2L(Cell const & first, Cell & second) const = 0;
    virtual void L2L(Cell const & first, Cell & second) const = 0;

};

#endif
