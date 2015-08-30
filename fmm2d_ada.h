/* 
 * File:   fmm2d_ada.h
 * Author: euklid
 *
 * Created on August 27, 2015, 6:39 PM
 */

#ifndef FMM2D_ADA_H
#define	FMM2D_ADA_H
#include "fmm2d.h"


class FMM2D_ADA : public FMM2D
{
public:
    FMM2D_ADA(std::vector<Element*> const & src_elements,
              std::vector<Element*> const & tgt_elements,
              unsigned int exp_terms,
              unsigned int loc_terms,
              unsigned int max_cell_elements,
              bool src_eq_tgt = false);

protected:
    virtual void build_tree();
    virtual void downward_pass();
private:
    virtual void moment_to_element_downward_pass(Cell * const target);
    virtual void add_locals_list4(Cell * const target);
};

#endif	/* FMM2DADA_H */

