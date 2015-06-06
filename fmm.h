#ifndef FMM_H
#define FMM_h
#include "element.h"

class FMM
{
public:
    FMM(std::vector<Element> const & elements ,unsigned int terms);

private:
    std::vector<Element> const & m_elements;
};

#endif
