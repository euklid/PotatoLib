#include "fmm.h"

FMM::FMM(std::vector<Element*> const & elements, unsigned int terms, unsigned int max_cell_elements) :
    m_terms(terms), m_max_cell_elements(max_cell_elements)
{
    m_elements = std::vector<Element*>(elements);
}

void FMM::set_kernel(Kernel const & kernel)
{
    m_kernel = &kernel;
}

