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

std::vector<Cell*> Cell2D::divide()
{
    double half_size = m_size/2;
    double quarter_size = m_size/4;
    std::vector<Cell*> new_cells;
    // don't want to assign an element to two or more cells


    // TODO: instead of searching the cell an element belongs to, use
    //       multiple sorts after their coordinates (smaller or larger than
    //       center coordinate) and use this information to directly assign
    //       the element to a cell. --> source: FMBEM book implementation details


    std::vector<bool> copied_elements(m_elements.size(),false);
    for (int i = 1; i>-2; i=i-2)
    {
        for (int j = -1; j<2; j=j+2)
        {
            //generate new cell
            Point new_center(2);
            new_center[0] = m_center[0] + i * quarter_size;
            new_center[1] = m_center[1] + j * quarter_size;
            Cell2D* new_cell = new Cell2D(m_terms,half_size, new_center);
            new_cells.push_back(new_cell);
            new_cell->set_father(this);
            //set its elements
            std::vector<Element*> new_cell_elements;
            for (int k = 0; k<m_elements.size(); k++)
            {
                if (copied_elements[k])
                {
                    continue;
                }
                if (new_cell->contains_point(m_elements[k]->get_position()))
                {
                    copied_elements[k] = true;
                    new_cell_elements.push_back(m_elements[k]);
                }
            }
        }
    }
    //set children
    m_children = new_cells;
    return new_cells;
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
