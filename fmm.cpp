#include "fmm.h"

FMM::FMM(std::vector<Element*> const & elements, unsigned int terms) :
    m_terms(terms)
{
    m_elements = std::vector<Element*>(elements);
}

void FMM::set_kernel(Kernel const & kernel)
{
    m_kern = &kernel;
}

