#include "kernel_laplace_2d.h"

double Laplace2DKernel::direct(Element const & el1, Element const & el2) const
{
    return direct_cmp(el1, el2).real;
}

complex_t Laplace2DKernel::direct_cmp(const Element &el1, const Element &pt) const
{
    // on identical points no calculation
    if (el1 == pt)
    {
        return 0;
    }
    
    return -el1.get_value()*std::log(Point::dist(el1.get_position(),pt.get_position()));
}

std::vector<complex_t> Laplace2DKernel::calc_moments_cmp(const std::vector<Element *> &elements,
                                                         const complex_t &mom_center,
                                                         int num_moments) const
{
    // using recursive multiplication get all exponentiations for a fixed element to compute
    // moments using M_k(z_c) = \sum_{i=1}^m q(z_i)(z_i-z_c)^k/k!
    std::vector<complex_t> moments(num_moments,0);
    for(int i = 0; i<elements.size(); i++)
    {
        complex_t contrib = elements[i]->get_value();
        const complex_t fac = complex_t(elements[i]->get_position())-mom_center;
        for(int j = 0; j<num_moments; j++)
        {
            moments[j]+= contrib;
            contrib*=fac/(j+1);
        }
    }
    return moments;
}

std::vector<double> Laplace2DKernel::calc_moments(const std::vector<Element *> &elements,
                                                  Point const & mom_center,
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
    loc_out.resize(loc_in.size(),0);
    // compute (z_l' - z_l)^k/k! for all k from 0 to mom_in.size()-1 to reuse them
    std::vector<complex_t> factors(loc_in.size(),0);
    const complex_t diff = loc_out_center-loc_in_center;
    factors[0] = complex_t(1);
    for(int i = 1; i<loc_in.size(); i++)
    {
        factors[i] = factors[i-1]*(diff/i);
    }

    for(int i = 0; i<loc_in.size(); i++)
    {
        for(int j = 0; j<loc_in.size()-i; j++)
        {
            loc_out[i] += factors[j]*loc_in[i+j];
        }
    }
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
    mom_out.resize(mom_in.size(),0);
    // compute (z_c - z_c')^k/k! for all k from 0 to mom_in.size()-1 to reuse them
    std::vector<complex_t> factors(mom_in.size(),0);
    const complex_t diff = mom_in_center-mom_out_center;
    factors[0] = complex_t(1);
    for(int i = 1; i<mom_in.size(); i++)
    {
        factors[i] = factors[i-1]*(diff/i);
    }

    for(int i = 0; i< mom_in.size(); i++)
    {
        for(int j = 0; j<=i; j++)
        {
            mom_out[i] += factors[i-j]*mom_in[j];
        }
    }
}

complex_t Laplace2DKernel::L2element_cmp(const std::vector<complex_t> &local_in, const complex_t &local_center, Element const &el) const
{
    complex_t dist = (complex_t)el.get_position() - local_center;
    unsigned int num_local_exp = local_in.size();
    complex_t fac;
    complex_t res = local_in[0];
    for(int i = 1; i<num_local_exp; i++)
    {
        fac*=dist/i;
        res+=fac*local_in[i];
    }
    return res;
}

double Laplace2DKernel::L2element(const std::vector<double> &local_in, Point loc_in_center, const Element &el) const
{
    return 0;
}

void Laplace2DKernel::M2L_cmp(std::vector<complex_t> const & moments,
                              complex_t const & mom_center,
                              std::vector<complex_t> & loc_exp,
                              complex_t const & loc_center) const
{
    loc_exp.resize(moments.size(),0);
    // compute (k-1)!/(z_l-z_c)^k for k>0 and -log(z_l-z_c) for k = 0 for reuse
    assert(mom_center != loc_center);

    const complex_t diff = loc_center-mom_center;
    std::vector<complex_t> factors(2*moments.size(), 0);
    factors[0] = -complex_t::log(diff);
    factors[1] = complex_t(1)/diff;
    for(int i=2; i<2*moments.size(); i++)
    {
        factors[i] = (factors[i-1]*(i-1))/diff;
    }

    int sign = 1;
    for(int i = 0; i<loc_exp.size(); i++)
    {
        for(int j = 0; j<moments.size(); j++)
        {
            loc_exp[i] += factors[i+j]*moments[i];
        }
        loc_exp[i] *= sign;
        sign *= -1;
    }
}

void Laplace2DKernel::M2L(std::vector<double> const & moments,
                          Point const & mom_center,
                          std::vector<double> & loc_exp,
                          Point const & loc_center) const
{}
