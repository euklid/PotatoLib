/* 
 * File:   fmm2d_ada.cpp
 * Author: euklid
 * 
 * Created on August 27, 2015, 6:39 PM
 */

#include "fmm2d_ada.h"
#include "cell2d.h"
#include "tree2d_ada.h"

FMM2D_ADA::FMM2D_ADA(const std::vector<Element*>& src_elements, 
                     const std::vector<Element*>& tgt_elements, 
                     unsigned int exp_terms, 
                     unsigned int loc_terms, 
                     unsigned int max_cell_elements, 
                     bool src_eq_tgt) 
: FMM2D(src_elements,
        tgt_elements,
        exp_terms,
        loc_terms,
        max_cell_elements,
        src_eq_tgt)                     
{}

void FMM2D_ADA::build_tree()
{
    if(m_tree)
    {
        delete m_tree;
    }
    m_tree = new Tree2D_Ada;
    std::pair<double,Point> bounding_box = get_bounding_cube(m_elements);
    Cell2D* root = new Cell2D(bounding_box.first, bounding_box.second);
    root->set_father(0);
    root->set_source_elements(m_src_elements);
    root->set_target_elements(m_tgt_elements);
    m_tree->set_root(root);
    m_tree->build_tree(m_max_cell_elements,1);
}

void FMM2D_ADA::moment_to_element_downward_pass(Cell* const target)
{
    std::vector<Cell*> list3 = target->get_list(2);
    unsigned int num_list_elements = list3.size();
    if(!num_list_elements) return;
    std::vector<Element*> const & tgt_elements = target->get_target_elements();
    unsigned int num_tgt_elements = tgt_elements.size();
    for(unsigned int i = 0; i<num_tgt_elements; i++)
    {
        Element* t = tgt_elements[i];
        double contrib = 0;
        for(unsigned int j = 0; j< num_list_elements; j++)
        {
            contrib += m_kernel->M2element_cmp(list3[j]->get_moments_cmp(),
                                               list3[j]->get_center(),
                                               *t).real;
        }
        t->set_target_value(t->get_target_value() + contrib);
    }
}

void FMM2D_ADA::add_locals_list4(Cell * const target)
{
    std::vector<Cell*> list4 = target->get_list(3);
    unsigned int num_list_elements = list4.size();
    if(!num_list_elements) return;
    std::vector<complex_t> loc_exps = target->get_local_exps_cmp();
    if(!loc_exps.size())
    {
        loc_exps.resize(m_loc_terms,0);
    }
    
    for(unsigned int i = 0; i<num_list_elements; i++)
    {
        std::vector<Element*> const & other_cell_el = list4[i]->get_source_elements();
        std::vector<complex_t> summand_loc_exps =
                m_kernel->calc_local_exp_cmp(other_cell_el,
                                             target->get_center(),
                                             m_loc_terms);
        add_moments(summand_loc_exps,loc_exps);
    }
    
    target->set_local_exps_cmp(loc_exps);
}

void FMM2D_ADA::downward_pass()
{
    Tree_Iterator *it = m_tree->downward_iterator();
    
    Cell* cur_cell;
    while(it->has_next())
    {
        cur_cell = it->next();
        if(cur_cell->get_level() > 1) break;
        if(cur_cell->get_target_elements().empty()) continue;
#ifdef DEBUG
    std::cout << "Downward pass: " << cur_cell->debug_info() << std::endl;
#endif
        // move moments to local expansion for list 2 (interaction list)
        m2l_downward_pass(cur_cell);
    
        // L2L only from level 2 on!
        
        if (cur_cell->is_leaf())
        {
            // go through neighbor list 1
            direct_downward_pass(cur_cell);
            // moment to element, list 3
            moment_to_element_downward_pass(cur_cell);
            // create local expansions about center of cur_cell for list 4 elements
            add_locals_list4(cur_cell);
            // evaluate local expansion in center of current cell for each 
            // element
            evaluate_far_interactions(cur_cell);
        }
    }
    // there are no level 2 cells
    if(cur_cell->get_level() < 2) 
    {
        delete it;
        return;
    }
    
    //reached lvl 3 or higher --> use additionally L2L
    while(1)
    {
        if(cur_cell->get_target_elements().empty())
        {
            // skip this and go to next if available
            if(it->has_next())
            {
                cur_cell = it->next();
                continue;
            }
            else //last element, leave loop
            {
                break;
            }
        }
#ifdef DEBUG
    std::cout << "Downward pass: " << cur_cell->debug_info() << std::endl;
#endif
        // move moments to local expansion for list 2 (interaction list)
        m2l_downward_pass(cur_cell);
        
        // L2L
        l2l_downward_pass(cur_cell);
        
        if (cur_cell->is_leaf())
        {
            // go through neighbor list 1
            direct_downward_pass(cur_cell);
            // moment to element, list 3
            moment_to_element_downward_pass(cur_cell);
            // create local expansions about center of cur_cell for list 4 elements
            add_locals_list4(cur_cell);
            
            // evaluate local expansion in center of current cell for each 
            // element
            evaluate_far_interactions(cur_cell);
        }
        if(it->has_next())
        {
            cur_cell = it->next();
        } 
        else //last element, leave loop
        {
            break;
        }
    }
    delete it;        
}