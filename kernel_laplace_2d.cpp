#include "kernel_laplace_2d.h"

double Laplace2DKernel::direct(Element const & el1, Element const & el2) const
{
    return direct_cmp(el1, el2).real;
}

complex_t Laplace2DKernel::direct_cmp(Element const & src, Element const & pt) const
{
    // on identical points no calculation
    if (src == pt)
    {
        return 0;
    }
    
    return -complex_t::log(pt.get_cmp_value());
}

std::vector<complex_t> Laplace2DKernel::calc_moments_cmp(const std::vector<Element *> &elements,
                                                         int num_moments) const
{
    
}

std::vector<double> Laplace2DKernel::calc_moments(const std::vector<Element *> &elements,
                                                  int num_moments) const
{
    return std::vector<double>();
}

void Laplace2DKernel::L2L(std::vector<double> const & loc_in,
                          Point const & loc_in_center,
                          std::vector<double> & loc_out,
                          Point const & loc_out_center) const
{}

void Laplace2DKernel::L2L_cmp(std::vector<complex_t> const & loc_in,
                              complex_t const & loc_in_center,
                              std::vector<complex_t> & loc_out,
                              complex_t const & loc_out_center) const
{
    
}

void Laplace2DKernel::M2M(std::vector<double> const & mom_in,
                          Point const & mom_in_center,
                          std::vector<double> & mom_out,
                          Point const & mom_out_center) const
{}

void Laplace2DKernel::M2M_cmp(std::vector<complex_t> const & mom_in,
                              complex_t const & mom_in_center,
                              std::vector<complex_t> & mom_out,
                              complex_t const & mom_out_center) const
{
    
}

void Laplace2DKernel::M2L_cmp(std::vector<complex_t> const & moments,
                              complex_t const & mom_center,
                              std::vector<complex_t> & loc_exp,
                              complex_t const & loc_center) const
{
    
}

void Laplace2DKernel::M2L(std::vector<double> const & moments,
                          Point const & mom_center,
                          std::vector<double> & loc_exp,
                          Point const & loc_center) const
{}


