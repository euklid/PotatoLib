#include "kernel_laplace_point_2d.h"
#ifdef DEBUG
#include <iostream>
#endif

double Laplace2DKernel::direct(Element const & target, Element const & src) const
{
    assert(src.get_type() & Element::SOURCE);
    assert(target.get_type() & Element::TARGET);
    // on identical points no calculation
    
    if(((PointElement&)src).get_position() == ((PointElement&)target).get_position())
    {
        return 0;
    }
    return -std::log(Point::dist(src.get_position(),target.get_position()));
}

std::vector<complex_t> Laplace2DKernel::calc_moments_cmp(const std::vector<Element *> &elements,
                                                         const complex_t &mom_center,
                                                         int num_moments) const
{
    // using recursive multiplication get all exponentiations for a fixed element to compute
    // moments using M_k(z_c) = \sum_{i=1}^m q(z_i)(z_i-z_c)^k/k!
    std::vector<complex_t> moments(num_moments,0);
    unsigned int num_el = elements.size();
    for(unsigned int i = 0; i<num_el; i++)
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
    assert(loc_in.size() > 0);
    assert(loc_in.size() == loc_out.size());
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
    assert(mom_in.size());
    if(mom_out.empty())
    {
        mom_out.resize(mom_in.size(),0);
    }
    // compute (z_c - z_c')^k/k! for all k from 0 to mom_in.size()-1 to reuse them
    std::vector<complex_t> factors(mom_in.size(),0);
    const complex_t diff = mom_in_center-mom_out_center;
    factors[0] = complex_t(1,0);
    unsigned int num_mom = mom_in.size();
    for(unsigned int i = 1; i<num_mom; i++)
    {
        factors[i] = factors[i-1]*(diff/i);
    }

    for(unsigned int i = 0; i< num_mom; i++)
    {
        for(unsigned int j = 0; j<=i; j++)
        {
            mom_out[i] += factors[i-j]*mom_in[j];
        }
    }
}

complex_t Laplace2DKernel::L2element_cmp(const std::vector<complex_t> &local_in, const complex_t &local_center, Element const &el) const
{
    complex_t dist = (complex_t)el.get_position() - local_center;
    unsigned int num_local_exp = local_in.size();
    complex_t fac(1,0);
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
    // compute (k-1)!/(z_l-z_c)^k for k>0 and -log(z_l-z_c) for k = 0 for reuse
    assert(mom_center != loc_center);

    const complex_t diff = loc_center-mom_center;
    std::vector<complex_t> factors(moments.size()+loc_exp.size(), 0);
    factors[0] = -complex_t::log(diff);
    factors[1] = complex_t(1)/diff;
    for(int i=2; i<factors.size(); i++)
    {
        factors[i] = (factors[i-1]*(i-1))/diff;
    }

    double sign = 1;
    for(int i = 0; i<loc_exp.size(); i++)
    {
        for(int j = 0; j<moments.size(); j++)
        {
            loc_exp[i] += factors[i+j]*moments[j];
        }
        loc_exp[i] *= sign;
        sign *= -1;
    }
}

void Laplace2DKernel::M2L(std::vector<double> const & moments,
                          Point const & mom_center,
                          std::vector<double> & loc_exp,
                          Point const & loc_center) const
{
}

std::vector<double> Laplace2DKernel::calc_local_exp(const std::vector<Element*>& elements, 
                                                    const Point& loc_center, 
                                                    int num_loc_exps) const
{
    return std::vector<double>();
}

std::vector<complex_t> Laplace2DKernel::calc_local_exp_cmp(const std::vector<Element*>& elements, 
                                                           const complex_t& loc_center, 
                                                           int num_loc_exps) const
{
    // using recursive multiplication get all exponentiations for a fixed element to compute
    // local expansions using L_k(z_l) = \sum_{i=1}^m q(z_i)/(z_i-z_l)^k*(k-1)! for k >=1
    // L_0(z_l) = \sum_{i=1}^m q(z_i)*-log(z_i-z_l)
    std::vector<complex_t> loc_exps(num_loc_exps,0);
    unsigned int num_el = elements.size();
    for(unsigned int i = 0; i<num_el; i++)
    {
        complex_t contrib = elements[i]->get_value();
        const complex_t fac = complex_t(elements[i]->get_position())-loc_center;
        loc_exps[0] += -contrib*(complex_t::log(fac));
        loc_exps[1] += (complex_t(1)*contrib)/fac;
        for(int j = 2; j<num_loc_exps; j++)
        {
            contrib*=complex_t(j-1)/fac;
            loc_exps[j]+= contrib;
        }
    }
    return loc_exps;
}

complex_t Laplace2DKernel::M2element_cmp(const std::vector<complex_t>& moments_in, const complex_t& moment_center, const Element& el) const
{
    const complex_t dist = complex_t(el.get_position()) - moment_center;
    unsigned int num_moments = moments_in.size();
    complex_t fac = complex_t(1,0)/dist;
    complex_t res = -complex_t::log(dist)*moments_in[0] + moments_in[1]*fac;
    for(int i = 2; i<num_moments; i++)
    {
        fac*=complex_t(i-1)/dist;
        res+=fac*moments_in[i];
    }
    return res;
}

double Laplace2DKernel::M2element(const std::vector<double>& moments_in, const Point& moment_center, const Element& el) const
{
    return 0;
}

