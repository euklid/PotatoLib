#ifndef FMM_H
#define FMM_H

#include "element.h"
#include "kernel.h"
#include "tree.h"
#include <map>

class FMM
{
public:
	
    FMM(std::vector<Element*> const & elements, 
        unsigned int terms, 
        unsigned int max_cell_elements);
    virtual void set_kernel(Kernel const & kernel);
    virtual void calculate() = 0;
    
    /**
     *  using the already constructed tree and given elements (that may have
     *  different source values) reevaluate the target values
     *  (source value, target value) of element
     */
    virtual void recalculate() = 0;

protected:
    virtual void build_tree() = 0;
    virtual void upward_pass() = 0;
    virtual void downward_pass() = 0;
    virtual void evaluate() = 0;

    std::vector<Element*> m_elements;
    const Kernel* m_kernel;
    Tree* m_tree;
    unsigned int m_terms;
    unsigned int m_max_cell_elements;
};

#endif
