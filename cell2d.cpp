/* 
 * File:   cell2d.cpp
 * Author: euklid
 * 
 * Created on June 9, 2015, 4:35 PM
 */

#include "cell2d.h"
#include <cmath>


Cell2D::Cell2D(unsigned int terms, double size, Point const & center)
    : Cell(2, terms, size, center) {}

std::vector<Cell*> & Cell2D::divide()
{
    double half_size = m_size/2;
    double quarter_size = m_size/4;
    // don't want to assign an element to two or more cells


    // TODO: instead of searching the cell an element belongs to, use
    //       multiple sorts after their coordinates (smaller or larger than
    //       center coordinate) and use this information to directly assign
    //       the element to a cell. --> source: FMBEM book implementation details

    std::vector<bool> copied_src_elements(m_src_elements.size(),false);
    std::vector<bool> copied_target_elements(m_target_elements.size(), false);
    for (int i = 1; i>-2; i=i-2)
    {
        for (int j = -1; j<2; j=j+2)
        {
            //generate new cell
            Point new_center(2);
            new_center[0] = m_center[0] + i * quarter_size;
            new_center[1] = m_center[1] + j * quarter_size;
            Cell2D* new_cell = new Cell2D(m_terms,half_size, new_center);
            new_cell->set_father(this);
            //set its elements
            std::vector<Element*> new_cell_src_elements;
            std::vector<Element*> new_cell_tgt_elements;
            for (int k = 0; k<m_src_elements.size(); k++)
            {
                if (!copied_src_elements[k])
                {
                    if (new_cell->contains_point(m_src_elements[k]->get_position()))
                    {
                        copied_src_elements[k] = true;
                        new_cell_src_elements.push_back(m_src_elements[k]);
                    }
                }
            }
            for (int k = 0; k<m_target_elements.size(); k++)
            {
                if (!copied_target_elements[k])
                {
                    if(new_cell->contains_point(m_target_elements[k]->get_position()))
                    {
                        copied_target_elements[k] = true;
                        new_cell_tgt_elements.push_back(m_target_elements[k]);
                    }
                }
            }

            if(new_cell->number_of_elements() == 0)
            {
                delete new_cell;
            }
            else
            {
                m_children.push_back(new_cell);
            }
        }
    }
    //set children
    return m_children;
}

bool Cell2D::contains_point(Point const &pt) const
{
    if (pt.get_dimension() != m_dim)
    {
        return false;
    }
    
    Point diff = m_center - pt;
    int half_size = m_size/2;
    for (int i = 0; i<m_dim; i++)
    {
        if (std::abs(diff[i]) >  half_size)
        {
            return false;
        }
    }
    
    return true;
}
