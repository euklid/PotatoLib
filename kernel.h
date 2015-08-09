#ifndef KERNEL_H
#define KERNEL_H

#include "element.h"
#include "complex_t.h"
#include <vector>

class Kernel
{
public:
    Kernel(){}
    virtual double direct(Element const & el1, Element const & el2) const = 0;
    virtual complex_t direct_cmp(Element const & el1, Element const & el2) const  = 0;
    virtual std::vector<complex_t> calc_moments_cmp(std::vector<Element*> const & elements,
                                                    int num_moments) const = 0;
    virtual std::vector<double> calc_moments(std::vector<Element*> const & elements,
                                             int num_moments) const = 0;
    virtual void M2L_cmp(std::vector<complex_t> const & moments,
                         complex_t const & mom_center,
                         std::vector<complex_t> & loc_exp,
                         complex_t const & loc_center) const = 0;
    virtual void M2L(std::vector<double> const & moments,
                     Point const & mom_center,
                     std::vector<double> & loc_exp,
                     Point const & loc_center) const = 0;
    virtual void L2L(std::vector<double> const & loc_in,
                     Point const & loc_in_center,
                     std::vector<double> & loc_out,
                     Point const & loc_out_center) const = 0;
    virtual void L2L_cmp(std::vector<complex_t> const & loc_in,
                     complex_t const & loc_in_center,
                     std::vector<complex_t> & loc_out,
                     complex_t const & loc_out_center) const = 0;
private:
    virtual void M2M(std::vector<double> const & mom_in,
                     Point const & mom_in_center,
                     std::vector<double> & mom_out,
                     Point const & mom_out_center) const = 0;
    virtual void M2M_cmp(std::vector<complex_t> const & mom_in,
                     complex_t const & mom_in_center,
                     std::vector<complex_t> & mom_out,
                     complex_t const & mom_out_center) const = 0;

};

#endif
