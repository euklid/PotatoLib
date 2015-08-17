#include "fmm.h"

FMM::FMM(std::vector<Element*> const & src_elements,
         std::vector<Element*> const & tgt_elements,
         unsigned int exp_terms,
         unsigned int loc_terms,
         unsigned int max_cell_elements) :
    m_src_elements(src_elements),
    m_tgt_elements(tgt_elements),
    m_tree(NULL),
    m_exp_terms(exp_terms),
    m_loc_terms(loc_terms),
    m_max_cell_elements(max_cell_elements)
{
    m_elements.insert(m_elements.begin(), m_src_elements.begin(), m_src_elements.end());
    m_elements.insert(m_elements.begin(), m_tgt_elements.begin(), m_tgt_elements.end());
}


std::vector<Element*> const & FMM::get_src_elements() const
{
    return m_src_elements;
    
}

std::vector<Element*> const & FMM::get_tgt_elements() const
{
    return m_tgt_elements;
}

void FMM::set_kernel(Kernel const & kernel)
{
    m_kernel = &kernel;
}