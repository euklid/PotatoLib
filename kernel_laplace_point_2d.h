#ifndef LAPLACE2D_H
#define LAPLACE2D_H

#include "kernel.h"
#include "point_element.h"

/**
 *  2-dimensional Laplace kernel for particle elements
 */
class Laplace2DKernel : public Kernel {

public:
    Laplace2DKernel(){}
    virtual double direct(Element const & target, Element const & src) const;
    virtual std::vector<complex_t> calc_moments_cmp(std::vector<Element*> const & elements,
                                                    complex_t const & mom_center,
                                                    int num_moments) const;
    virtual std::vector<double> calc_moments(std::vector<Element*> const & elements,
                                             const Point &mom_center,
                                             int num_moments) const;
    virtual void M2L_cmp(std::vector<complex_t> const & moments,
                         complex_t const & mom_center,
                         std::vector<complex_t> & loc_exp,
                         complex_t const & loc_center) const;
    virtual void M2L(std::vector<double> const & moments,
                     Point const & mom_center,
                     std::vector<double> & loc_exp,
                     Point const & loc_center) const;
    virtual void L2L(std::vector<double> const & loc_in,
                     Point const & loc_in_center,
                     std::vector<double> & loc_out,
                     Point const & loc_out_center) const;
    virtual void L2L_cmp(std::vector<complex_t> const & loc_in,
                         complex_t const & loc_in_center,
                         std::vector<complex_t> & loc_out,
                         complex_t const & loc_out_center) const;
    virtual void M2M(std::vector<double> const & mom_in,
                     Point const & mom_in_center,
                     std::vector<double> & mom_out,
                     Point const & mom_out_center) const;
    virtual void M2M_cmp(std::vector<complex_t> const & mom_in,
                         complex_t const & mom_in_center,
                         std::vector<complex_t> & mom_out,
                         complex_t const & mom_out_center) const;
    virtual complex_t L2element_cmp(std::vector<complex_t> const & local_in,
                               complex_t const & local_center,
                               Element const & el) const ;
    virtual double L2element(std::vector<double> const & local_in,
                           Point loc_in_center,
                           Element const & el) const;
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

#endif
