#include "fmm.h"

FMM::FMM(std::vector<Element*> const & src_elements,
         std::vector<Element*> const & tgt_elements,
         unsigned int exp_terms,
         unsigned int loc_terms,
         unsigned int max_cell_elements,
         bool src_eq_tgt) :
    m_src_elements(src_elements),
    m_tgt_elements(tgt_elements),
    m_tree(0),
    m_exp_terms(exp_terms),
    m_loc_terms(loc_terms),
    m_max_cell_elements(max_cell_elements),
    m_make_prec(false)
{
    // we assume that source and target elements have ascending order by id
    for(int i = 0; i<src_elements.size(); i++)
    {
        assert(src_elements[i]->get_id() == i);
    }
    if(!src_eq_tgt)
    {
        for(int i = 0; i<tgt_elements.size(); i++)
        {
            assert(tgt_elements[i]->get_id() == i);
        }
    } else {
        // source, that is target, is counted twice --> double criterion for 
        // subdivision of cell
        m_max_cell_elements =max_cell_elements << 1;
    }
    
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

const std::vector<arma::mat> &FMM::get_precond() const
{
   if(!m_precond.size())
   {
       std::cerr << "No preconditioner has been computed" << std::endl;
       return m_precond;
   }
   assert(m_precond.size() == m_tree->get_leaves().size());
   return m_precond;
}

std::vector<Element*> const & FMM::get_permuted_elements()
{
    if (!m_perm_el.empty())
    {
        return m_perm_el;
    }

    std::vector<Cell*> const & leaves = m_tree->get_leaves();
    std::vector<Cell*>::const_iterator it = leaves.begin();
    m_perm_el.reserve(m_src_elements.size());
    for (; it != leaves.end(); ++it)
    {
        std::vector<Element*> const & cell_el = (*it)->get_source_elements();
        m_perm_el.insert(m_perm_el.end(), cell_el.begin(), cell_el.end());
    }
    
    m_permutation.reserve(m_perm_el.size());
    for(int i = 0; i<m_permutation.size(); i++)
    {
        unsigned int id = m_perm_el[i]->get_id();
        assert(id < m_perm_el.size());
        m_permutation[i] = id;
    }
    assert(m_perm_el.size() == m_src_elements.size());
    return m_perm_el;
}
    
const std::vector<unsigned int>& FMM::get_element_permutation()
{
    if(!m_permutation.empty())
    {
        return m_permutation;
    }
    
    get_permuted_elements();
    return m_permutation;
}

const std::vector<unsigned int> & FMM::get_prec_block_starts() const
{
    return m_prec_block_starts;
}
