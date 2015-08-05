#ifndef LAPLACE2D_H
#define LAPLACE2D_H

#include "kernel.h"

/**
 *  2-dimensional Laplace kernel for particle elements
 */
class Laplace2DKernel : public Kernel {

public:
    Laplace2DKernel(){}
    virtual double direct(Element const & el1, Element const & el2) const;
    virtual complex_t direct_cmp(Element const & el1, Element const & el2) const;
    virtual std::vector<complex_t> calc_moments_cmp(Cell const & cell, int num_moments) const;
    virtual std::vector<double> calc_moments(Cell const & cell, int num_moments) const;
    virtual void M2M(Cell const & first, Cell & second) const;
    virtual void M2L(Cell const & first, Cell & second) const;
    virtual void L2L(Cell const & first, Cell & second) const;
};

#endif