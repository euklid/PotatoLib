/* 
 * File:   FMM2D.h
 * Author: euklid
 *
 * Created on June 7, 2015, 3:51 PM
 */

#ifndef FMM2D_H
#define	FMM2D_H

#include "fmm.h"

class FMM2D : public FMM
{
public:
	FMM2D(std::vector<Element*> const & elements, unsigned int terms);
	virtual std::vector<double> calculate();
private:
	virtual void build_tree();
	virtual void upward_pass();
	virtual void downward_pass();

};

#endif	/* FMM2D_H */

