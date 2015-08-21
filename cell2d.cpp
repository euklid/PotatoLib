/* 
 * File:   cell2d.cpp
 * Author: euklid
 * 
 * Created on June 9, 2015, 4:35 PM
 */

#include "cell2d.h"
#include <cmath>


Cell2D::Cell2D( double size, Point const & center)
    : Cell(2, size, center) {}

std::vector<Cell*> & Cell2D::divide()
{
    double half_size = m_size/2;
    double quarter_size = m_size/4;
    // don't want to assign an element to two or more cells


    //       instead of searching the cell an element belongs to, use
    //       multiple sorts after their coordinates (smaller or larger than
    //       center coordinate) and use this information to directly assign
    //       the element to a cell.

    //std::vector<bool> copied_src_elements(m_src_elements.size(),false);
    //std::vector<bool> copied_target_elements(m_target_elements.size(), false);
    unsigned int num_src_elements = m_src_elements.size();
    unsigned int num_tgt_elements = m_target_elements.size();
    std::vector<Cell*> new_cells;
    std::vector<std::vector<Element*> > new_cells_src_elements;
    std::vector<std::vector<Element*> > new_cells_tgt_elements;
    int max_childs = 1 << m_dim;

    //create buckets == cells
    for(int i = 0; i<max_childs; i++)
    {
        Point new_center(m_dim);
        for(int j = 0; j<m_dim; j++)
        {
            int dir = (((i & (1<<j))>>j)<<1)-1; //1 if (j+1)-lowest bit is 1, otherwise -1
            new_center[j] = m_center[j] + dir * quarter_size;
        }
        Cell2D* new_cell = new Cell2D(half_size, new_center);
        new_cell->set_father(this);
        new_cells.push_back(new_cell);
        new_cells_src_elements.push_back(std::vector<Element*>());
        new_cells_tgt_elements.push_back(std::vector<Element*>());
    }
    //identify buckets for each element
    int cell_bucket ;
    for(int j = 0; j<num_src_elements; j++)
    {
        cell_bucket = 0;
        for(int dim = 0; dim<m_dim; dim++)
        {
            cell_bucket+= (m_src_elements[j]->get_position()[dim] > m_center[dim]) << dim;
        }
        new_cells_src_elements[cell_bucket].push_back(m_src_elements[j]);
    }
    for(int j = 0; j<num_tgt_elements; j++)
    {
        cell_bucket = 0;
        for(int dim = 0; dim<m_dim; dim++)
        {
            cell_bucket+=(m_target_elements[j]->get_position()[dim] > m_center[dim]) << dim;
        }
        new_cells_tgt_elements[cell_bucket].push_back(m_target_elements[j]);
    }

    for(int i = 0; i<max_childs; i++)
    {
        new_cells[i]->set_target_elements(new_cells_tgt_elements[i]);
        new_cells[i]->set_source_elements(new_cells_src_elements[i]);
        if(new_cells[i]->number_of_elements())
        {
            m_children.push_back(new_cells[i]);
        }
    }
    return m_children;
    /*for (int i = 1; i>-2; i=i-2)
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
    */
}

bool Cell2D::contains_point(Point const &pt) const
{
    if (pt.get_dimension() != m_dim)
    {
        return false;
    }
    
    Point diff = m_center - pt;
    double half_size = m_size/2;
    for (int i = 0; i<m_dim; i++)
    {
        double abs_diff = diff[i];
        abs_diff = (abs_diff>0)?abs_diff:-abs_diff;
        if (std::abs(diff[i]) >  half_size)
        {
            return false;
        }
    }
    
    return true;
}
