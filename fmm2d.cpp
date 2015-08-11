/* 
 * File:   FMM2D.cpp
 * Author: euklid
 * 
 * Created on June 7, 2015, 3:51 PM
 */

#include "fmm2d.h"
#include "cell2d.h"

std::pair<double,Point> get_bounding_cube(const std::vector<Element*>& elements);

FMM2D::FMM2D(std::vector<Element*> const & elements, 
        unsigned int terms, 
        unsigned int max_cell_elements) : 
	FMM(elements, terms, max_cell_elements) 
{
}

std::vector<double> FMM2D::calculate() 
{
    build_tree();
    upward_pass();
    downward_pass();
}

void FMM2D::build_tree() 
{
    m_tree = new Tree2D;
    std::pair<double,Point> bounding_box = get_bounding_cube(m_elements);
    Cell2D* root = new Cell2D(m_terms,bounding_box.first, bounding_box.second);
    root->set_father(0);
    root->set_elements(m_elements);
    m_tree->set_root(root);
    m_tree->build_tree(m_max_cell_elements);
}

void FMM2D::upward_pass() 
{
    Tree_Iterator *it = m_tree->upward_iterator();
    while (it->has_next())
    {
        Cell *cur_cell = it->next();

        // leaf node has to calculate its moments
        if (cur_cell->is_leaf())
        {
            std::vector<complex_t> moments = m_kernel->calc_moments_cmp(cur_cell->get_elements(),
                                                                        cur_cell->get_center(),
                                                                        m_terms);
            cur_cell->set_moments_cmp(moments);
        }
        else
        {
            // sum up translated moments from children
            std::vector<Cell*> const & children = cur_cell->get_children();
            Cell* cur_child = children[0];
            std::vector<complex_t> moments;
            m_kernel->M2M_cmp(cur_child->get_moments_cmp(),
                              cur_child->get_center(),
                              moments,
                              cur_cell->get_center());
            unsigned int num_children = children.size();
            for(unsigned int i = 1; i<num_children; i++)
            {
                cur_child = children[i];
                std::vector<complex_t> mom_summand;
                m_kernel->M2M_cmp(cur_child->get_moments_cmp(),
                                  cur_child->get_center(),
                                  mom_summand,
                                  cur_cell->get_center());
                assert(mom_summand.size() == m_terms);
                for(int j = 0; j<m_terms; j++)
                {
                    moments[j] += mom_summand[j];
                }
            }
            cur_cell->set_moments_cmp(moments);
        }
    }
}

void FMM2D::downward_pass() 
{
	
}

/**
 * get smallest cube containing all elements
 * @param elements all n-dimensional elements
 * @return pair of cube length and center coordinate
 */
std::pair<double, Point> get_bounding_cube(std::vector<Element*> const & elements) 
{
    std::vector<Element*>::const_iterator it = elements.begin();
    Point min, max, center;
    min = elements.front()->get_position();
    max = min;
    for(;it!=elements.end(); ++it) 
    {
        Element* cur_el = *it;
        for(int i = 0; i<min.get_dimension(); i++)
        {
            Point cur_el_pos = cur_el->get_position();
            if (cur_el_pos[i] < min[i])
            {
                min[i] = cur_el_pos[i];
            }
            if (cur_el_pos[i] > max[i])
            {
                max[i] = cur_el_pos[i];
            }
        }
    }

    double max_dist = 0;
    for (int i = 0; i< min.get_dimension(); i++) 
    {
        double dist = max[i]-min[i];
        center[i] = (max[i]+min[i])/2;
        if (dist > max_dist) 
        {
            max_dist = dist;
        }
    }
    
    //extend cube a little bit in all directions to be sure to contain all elements
    max_dist = 1.02*max_dist;
    return std::make_pair(max_dist,center);
}
