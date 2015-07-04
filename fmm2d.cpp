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
        if (dist > max_dist) 
        {
            max_dist = dist;
        }
    }
    
    center = min;
    double half_dist = max_dist/2;
    for (int i = 0; i< min.get_dimension(); i++) 
    {
        center[i] += half_dist;
    }
    
    return std::make_pair(max_dist,center);
}