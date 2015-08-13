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
            std::vector<complex_t> moments = m_kernel->calc_moments_cmp(cur_cell->get_source_elements(),
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
template<class T>
void add_moments_cmp(std::vector<T> const & summand, std::vector<T> & moments)
{
    assert(summand.size() == moments.size());
    unsigned int num_moments = moments.size();
    for(unsigned int i = 0; i<num_moments; i++)
    {
        moments[i] += summand[i];
    }
}

void FMM2D::m2l_downward_pass(Cell* cur_cell)
{
    std::vector<complex_t> shifted_moments, local_exps = cur_cell->get_local_exps_cmp();
    std::vector<Cell*> const & interaction_list = cur_cell->get_interaction_list();
    unsigned int num_interaction_list = interaction_list.size();
    for(unsigned int i = 0; i<num_interaction_list; i++)
    {
        m_kernel->M2L_cmp(interaction_list[i]->get_moments_cmp(),
                          interaction_list[i]->get_center(),
                          shifted_moments,
                          cur_cell->get_center());
        add_moments_cmp(shifted_moments,local_exps);
    }
    cur_cell->set_local_exps_cmp(local_exps);
}

void FMM2D::l2l_downward_pass(Cell *cell)
{
    std::vector<complex_t> local_exp = cell->get_local_exps_cmp();
    std::vector<complex_t> shifted_local_exp(m_terms,0);
    Cell* father = cell->get_father();
    std::vector<complex_t> const & father_local_exp = father->get_local_exps_cmp();
    m_kernel->L2L_cmp(father_local_exp,father->get_center(),shifted_local_exp,cell->get_center());
    add_moments_cmp(shifted_local_exp,local_exp);
    cell->set_local_exps_cmp(local_exp);
}

void FMM2D::direct_downward_pass(Cell *target, Cell *source)
{
    std::vector<Element*> & target_elements = target->get_target_elements();
    std::vector<Element*> const & source_elements = source->get_source_elements();
    unsigned int num_tgt_el = target_elements.size();
    unsigned int num_src_el = source_elements.size();
    for(unsigned int i = 0; i < num_tgt_el; i++)
    {
        double contrib = 0;
        Element * t = target_elements[i];
        for(unsigned int j = 0; j < num_src_el; j++)
        {
            Element * s = source_elements[j];
            contrib += m_kernel->direct(*t,*s);
        }
        t->set_target_value(t->get_target_value()+contrib);
    }
}

void FMM2D::evaluate_far_interactions(Cell* cell)
{
    std::vector<Element*> & target_elements = cell->get_target_elements();
    complex_t cell_center = cell->get_center();
    std::vector<complex_t> const & local_exps = cell->get_local_exps_cmp();
    complex_t contrib;
    for(int i = 0; i< target_elements.size(); i++)
    {
        contrib = m_kernel->L2element_cmp(local_exps,cell_center,*(target_elements[i]));
        target_elements[i]->set_target_value(target_elements[i]->get_target_value() + contrib.real);
    }
}

void FMM2D::downward_pass() 
{
    Tree_Iterator* it = m_tree->downward_iterator();
    // level 2 only M2L and direct, no L2L
    Cell* cur_cell;
    while(it->has_next())
    {
        cur_cell = it->next();
        if(cur_cell->get_level() > 2) break;
        //M2L
        m2l_downward_pass(cur_cell);
        // FIXME direct evaluation

    }

    //reached lvl 3 or higher --> use L2L and M2L
    do
    {
        m2l_downward_pass(cur_cell);
        l2l_downward_pass(cur_cell);
        // FIXME direct evaluation

        // FIXME local expansions to cell if cell is leaf

        cur_cell = it->next();
    } while (it->has_next());
}

void FMM2D::evaluate()
{
    // in Liu's version we do the direct evaluation while we are
    // moving the tree downward because of the way the interaction
    // and direct lists are built
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
