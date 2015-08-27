/* 
 * File:   kernel_laplace_element_2d.h
 * Author: euklid
 *
 * Created on August 24, 2015, 6:16 PM
 */

#ifndef KERNEL_LAPLACE_CONSTANT_ELEMENT_2D_H
#define	KERNEL_LAPLACE_CONSTANT_ELEMENT_2D_H

#include "kernel_laplace_point_2d.h"


class KernLapConstEl2D : public Laplace2DKernel
{
public:
    KernLapConstEl2D()
    {}
    virtual double direct(Element const & target, Element const & src) const ;
    virtual std::vector<complex_t> calc_moments_cmp(std::vector<Element*> const & elements,
                                                    complex_t const & mom_center,
                                                    int num_moments) const;
    virtual std::vector<double> calc_moments(std::vector<Element*> const & elements,
                                             Point const & mom_center,
                                             int num_moments) const;

    virtual complex_t L2element_cmp(const std::vector<complex_t>& local_in, 
                                    const complex_t& local_center, 
                                    const Element& el) const;
    virtual std::vector<complex_t> calc_local_exp_cmp(std::vector<Element*> const & elements,
                                                      complex_t const & loc_center,
                                                      int num_loc_exps) const;
    virtual std::vector<double> calc_local_exp(std::vector<Element*> const & elements,
                                               Point const & loc_center,
                                               int num_loc_exps) const;
    virtual complex_t M2element_cmp(std::vector<complex_t> const & moments_in,
                                    complex_t const & moment_center,
                                    Element const & el) const;
    virtual double M2element(std::vector<double> const &  moments_in,
                             Point const & moment_center,
                             Element const & el) const;


};

#endif	/* KERNEL_LAPLACE_ELEMENT_2D_H */

