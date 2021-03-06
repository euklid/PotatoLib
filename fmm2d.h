/* 
 * File:   FMM2D.h
 * Author: euklid
 *
 * Created on June 7, 2015, 3:51 PM
 */

#ifndef FMM2D_H
#define	FMM2D_H

#include "fmm.h"
#include "tree2d.h"

class FMM2D : public FMM
{
public:
    FMM2D(std::vector<Element*> const & src_elements,
          std::vector<Element*> const & tgt_elements,
          unsigned int exp_terms,
          unsigned int loc_terms,
          unsigned int max_cell_elements,
          bool src_eq_tgt = false);
    virtual void calculate(bool precond = false);
    virtual void recalculate();

protected:
    virtual void build_tree();
    virtual void upward_pass();
    virtual void downward_pass();
    virtual void evaluate();
    virtual void reset();
    virtual void m2l_downward_pass(Cell* cell);
    virtual void l2l_downward_pass(Cell* cell);
    virtual void direct_downward_pass(Cell* target);
    virtual void evaluate_far_interactions(Cell* cell);
    virtual void init_precond();

    template<class T>
    static void add_moments(std::vector<T> const & summand,
                            std::vector<T> & moments)
    {
#if DEBUG
        assert(summand.size() == moments.size());
#endif
        unsigned int num_moments = moments.size();
        for (unsigned int i = 0; i < num_moments; i++)
        {
            moments[i] += summand[i];
        }
    }

};

#endif	/* FMM2D_H */

