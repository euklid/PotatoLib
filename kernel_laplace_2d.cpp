#include "kernel_laplace_2d.h"

double Laplace2DKernel::direct(Element const & el1, Element const & el2) const
{
    return direct_cmp(el1, el2).real;
}

complex_t Laplace2DKernel::direct_cmp(Element const & el1, Element const & el2) const
{
    // on same object no calculation
    if (&el1 == &el2)
    {
        return 0;
    }
    
    
    
    
}

std::vector<complex_t> Laplace2DKernel::calc_moments_cmp(Cell const & cell, int num_moments) const
{
    return std::vector<complex_t>();
}

std::vector<double> Laplace2DKernel::calc_moments(Cell const & cell, int num_moments) const
{
    return std::vector<double>();
}

void Laplace2DKernel::M2M(Cell const & first, Cell & second) const
{
    
}

void Laplace2DKernel::M2L(Cell const & first, Cell & second) const
{
    
}

void Laplace2DKernel::L2L(Cell const & first, Cell & second) const
{
    
}