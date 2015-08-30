#ifndef FMM_H
#define FMM_H

#include "element.h"
#include "kernel.h"
#include "tree.h"
#include <vector>
#include <armadillo>

class FMM
{
public:
	
    FMM(std::vector<Element*> const & src_elements,
        std::vector<Element*> const & tgt_elements,
        unsigned int exp_terms,
        unsigned int loc_terms,
        unsigned int max_cell_elements,
        bool src_eq_tgt = false);
    
    virtual void set_kernel(Kernel const & kernel);
    
    /**
     * calculates the target values of target elements and if wanted a
     * preconditioner for GMRES.
     * @param precond If true, calculate a preconditioning matrix.
     * When computing the preconditioner it is assumed that each element has 
     * source and target values (e.g. a composition of an element that has a
     * node on it that is the target for potential calculations
     * ----*----   ------- is source of element, * is target node on it)
     */
    virtual void calculate(bool precond = false) = 0;
    virtual std::vector<Element*> const & get_src_elements() const;
    virtual std::vector<Element*> const & get_tgt_elements() const;
    
    /**
     * If preconditioner should've been calculated, return it as array of
     * block matrices it consists of 
     * It is suggested to use the get_permuted_elements() to find out how
     * the entries in the preconditioner relate to the original order of 
     * elements
     * @return preconditioner sparse matrix for preconditioning GMRES
     */
    virtual std::vector<arma::mat> const & get_precond() const;
    /**
     *  using the already constructed tree and given elements (that may have
     *  different source values) reevaluate the target values
     *  (source value, target value) of element
     */
    virtual void recalculate() = 0;
    
    /**
     * Assuming that each source element has a target node on it like described
     * in calculate(bool), return the permuted elements how they
     * occur in the preconditioning matrix.
     * @return permutation of elements
     */
    virtual std::vector<Element*> const & get_permuted_elements();
    
     /**
     * Assuming that each source element has a target node on it like described
     * in calculate(bool), return the permutation of the elements how they
     * occur in the preconditioning matrix.
     * @return permutation of elements
     */
    virtual std::vector<unsigned int> const & get_element_permutation();
    
    virtual std::vector<unsigned int> const & get_prec_block_starts() const;

protected:
    virtual void build_tree() = 0;
    virtual void upward_pass() = 0;
    virtual void downward_pass() = 0;
    virtual void evaluate() = 0;
    virtual void reset() = 0;
    
    /**
     * get smallest cube containing all elements
     * @param elements all n-dimensional elements
     * @return pair of cube length and center coordinate
     */
    static std::pair<double, Point> get_bounding_cube(std::vector<Element*> const & elements);
    
    std::vector<Element*> const & m_src_elements;
    std::vector<Element*> const & m_tgt_elements;
    std::vector<Element*> m_elements;
    std::vector<Element*> m_perm_el;
    std::vector<unsigned int> m_permutation;
    const Kernel* m_kernel;
    Tree* m_tree;
    unsigned int m_exp_terms;
    unsigned int m_loc_terms;
    unsigned int m_max_cell_elements;
    bool m_make_prec;
    std::vector<arma::mat> m_precond;
    std::vector<unsigned int> m_prec_block_starts;
};

#endif
