/* 
 * File:   kernel_laplace_constant_element_2d.cpp
 * Author: euklid
 * 
 * Created on August 24, 2015, 6:16 PM
 */

#include "kernel_laplace_constant_element_2d.h"
#include "constant_element_2d.h"
#ifdef DEBUG
#include <iostream>
#endif

double KernLapConstEl2D::direct(Element const & target, Element const & src) const
{
    assert(target.get_type() & Element::TARGET && src.get_type() & Element::SOURCE);
    ConstEl2D const & s = static_cast<const ConstEl2D&>(src);
    Point const & end_node = s.get_end_node();
    Point const & start_node = s.get_start_node();
    Point const & t_pos = target.get_position();
    Point r1_vec = start_node - t_pos;
    Point r2_vec = end_node - t_pos;
    
    double el_length = s.get_length();
    double r1;
    if(t_pos == s.get_position())
    {
        r1 = el_length/2;
        return r1*(1-std::log(r1))*M_1_PI;
    }
    
    r1 = Point::dist(t_pos,start_node);
    double r2 = Point::dist(t_pos,end_node);
    Point const & norm = s.get_el_norm();
    double d = std::fabs(r1_vec[0]*norm[0] + r1_vec[1]*norm[1]);
    double t1 = -r1_vec[0]*norm[1]+r1_vec[1]*norm[0];
    double t2 = -r2_vec[0]*norm[1]+r2_vec[1]*norm[0];
    double theta = std::atan2(d*el_length,d*d+t1*t2);
    
    return M_1_PI/2*(-theta*d + el_length+t1*std::log(r1) - t2*std::log(r2));
}

std::vector<complex_t> KernLapConstEl2D::calc_moments_cmp(std::vector<Element*> const & elements,
                                                          complex_t const & mom_center,
                                                          int num_moments) const
{
    // using recursive multiplication get all exponentiations for a fixed element to compute
    // moments using M_k(z_c) = \sum_{i=1}^n q(e_i)*conj(tangent_normed_i)*{(end_i-center)^(k+1)/(k+1)! - (start_i-center)^(k+1)/(k+1)!}
    std::vector<complex_t> moments(num_moments,0);
    unsigned int num_el = elements.size();
    for(unsigned int i = 0; i<num_el; i++)
    {
        ConstEl2D* el = dynamic_cast<ConstEl2D*>(elements[i]);
        complex_t tangent_norm_conj(el->get_end_node()-el->get_start_node());
        tangent_norm_conj = tangent_norm_conj/ tangent_norm_conj.abs();
        tangent_norm_conj = tangent_norm_conj.conj();
        
        const complex_t fac_tangent_norm = tangent_norm_conj*el->get_value();
        
        complex_t fac_end(el->get_end_node());
        fac_end = fac_end - mom_center;
        
        complex_t fac_start(el->get_start_node());
        fac_start = fac_start - mom_center;
        
        complex_t contrib_p = fac_tangent_norm;
        complex_t contrib_n = fac_tangent_norm;
        for(int j = 0; j<num_moments; j++)
        {
            contrib_p *= fac_end/(j+1);
            contrib_n *= fac_start/(j+1);
            moments[j]+= contrib_p - contrib_n;
        }
    }
    return moments;
}
std::vector<double> KernLapConstEl2D::calc_moments(std::vector<Element*> const & elements,
                                                   Point const & mom_center,
                                                   int num_moments) const
{
    return std::vector<double>();
}

complex_t KernLapConstEl2D::L2element_cmp(const std::vector<complex_t>& local_in,
                                          const complex_t& local_center,
                                          const Element& el) const
{
    return  Laplace2DKernel::L2element_cmp(local_in, local_center, el)*(M_1_PI/2);
}
