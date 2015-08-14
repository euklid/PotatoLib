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
	FMM2D(std::vector<Element*> const & elements, 
         unsigned int terms, 
         unsigned int max_cell_elements);
	virtual std::vector<double> calculate();
private:
	virtual void build_tree();
	virtual void upward_pass();
	virtual void downward_pass();
    virtual void evaluate();
    virtual void m2l_downward_pass(Cell* cell);
    virtual void l2l_downward_pass(Cell* cell);
    virtual void direct_downward_pass(Cell* target);
    virtual void evaluate_far_interactions(Cell* cell);

};

#endif	/* FMM2D_H */

