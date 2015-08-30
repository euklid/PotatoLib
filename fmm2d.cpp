/* 
 * File:   FMM2D.cpp
 * Author: euklid
 * 
 * Created on June 7, 2015, 3:51 PM
 */

#include "fmm2d.h"
#include "cell2d.h"

std::pair<double,Point> get_bounding_cube(const std::vector<Element*>& elements);

FMM2D::FMM2D(std::vector<Element*> const & src_elements,
             std::vector<Element*> const & tgt_elements,
             unsigned int exp_terms,
             unsigned int loc_terms,
             unsigned int max_cell_elements,
             bool src_eq_tgt) :
	FMM(src_elements, 
         tgt_elements, 
         exp_terms, 
         loc_terms, 
         max_cell_elements, 
         src_eq_tgt)
{
}

void FMM2D::calculate(bool precond)
{
    reset();
    m_make_prec = precond;
    build_tree();
    upward_pass();
    downward_pass();
    m_make_prec = false;
}

void FMM2D::recalculate()
{
    if (!m_tree) {
        build_tree();
    }
    reset();
    upward_pass();
    downward_pass();
}

void FMM2D::build_tree() 
{
    if(m_tree)
    {
        delete m_tree;
    }
    m_tree = new Tree2D;
    std::pair<double,Point> bounding_box = get_bounding_cube(m_elements);
    Cell2D* root = new Cell2D(bounding_box.first, bounding_box.second);
    root->set_father(0);
    root->set_source_elements(m_src_elements);
    root->set_target_elements(m_tgt_elements);
    m_tree->set_root(root);
    m_tree->build_tree(m_max_cell_elements);
}

void FMM2D::upward_pass() 
{
    Tree_Iterator *it = m_tree->upward_iterator();
    while (it->has_next())
    {
        Cell *cur_cell = it->next();
#ifdef DEBUG
        std::cout << "current upward pass cell is " << cur_cell->debug_info() << std::endl;
#endif

        // leaf node has to calculate its moments
        if(cur_cell->get_source_elements().empty()) 
        {
            cur_cell->set_moments_cmp(std::vector<complex_t>(m_exp_terms,0));
            continue;
        }
        if (cur_cell->is_leaf())
        {
            std::vector<complex_t> moments = m_kernel->calc_moments_cmp(cur_cell->get_source_elements(),
                                                                        cur_cell->get_center(),
                                                                        m_exp_terms);
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
                std::vector<complex_t> mom_summand(m_exp_terms,0);
                m_kernel->M2M_cmp(cur_child->get_moments_cmp(),
                                  cur_child->get_center(),
                                  mom_summand,
                                  cur_cell->get_center());
                assert(mom_summand.size() == m_exp_terms);
                for(int j = 0; j<m_exp_terms; j++)
                {
                    moments[j] += mom_summand[j];
                }
            }
            cur_cell->set_moments_cmp(moments);
        }
    }
    delete it;
}

template<class T>
void FMM2D::add_moments(std::vector<T> const & summand, std::vector<T> & moments)
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
    std::vector<complex_t> local_exps = cur_cell->get_local_exps_cmp();
    if(!local_exps.size())
    {
        local_exps.resize(m_loc_terms,0);
    }
    std::vector<Cell*> const & interaction_list = cur_cell->get_list(1);
    unsigned int num_interaction_list = interaction_list.size();
    for(unsigned int i = 0; i<num_interaction_list; i++)
    {
        std::vector<complex_t> shifted_moments(m_loc_terms,0);
        if(interaction_list[i]->get_source_elements().empty())
        {
            continue;
        }
        m_kernel->M2L_cmp(interaction_list[i]->get_moments_cmp(),
                          interaction_list[i]->get_center(),
                          shifted_moments,
                          cur_cell->get_center());
        add_moments(shifted_moments,local_exps);
    }
    /*printf("M2L for cell: %d at level %d with posx %f and posy %f\n",
           cur_cell->get_id(), 
           cur_cell->get_level(),
           cur_cell->get_level_grid_position()[0], 
           cur_cell->get_level_grid_position()[1]);
    fflush(stdin);
    */
    cur_cell->set_local_exps_cmp(local_exps);
}

void FMM2D::l2l_downward_pass(Cell *cell)
{
    std::vector<complex_t> local_exp = cell->get_local_exps_cmp();
    if(local_exp.empty())
    {
        local_exp.resize(m_loc_terms,0);
    }
    std::vector<complex_t> shifted_local_exp(m_loc_terms,0);
    Cell* father = cell->get_father();
    std::vector<complex_t> const & father_local_exp = father->get_local_exps_cmp();
    m_kernel->L2L_cmp(father_local_exp,
                      father->get_center(),
                      shifted_local_exp,
                      cell->get_center());
    add_moments(shifted_local_exp,local_exp);
    cell->set_local_exps_cmp(local_exp);
}

/**
 * @brief FMM2D::direct_downward_pass evaluates all direct interactions of cell target
 *        with its interaction list
 * @param target target cell
 */
void FMM2D::direct_downward_pass(Cell *target)
{
    std::vector<Cell*> direct_neighbours = target->get_list(0);
    unsigned int num_dir_neighbours = direct_neighbours.size();
    if (num_dir_neighbours == 0) return;
    std::vector<Element*> const & target_elements = target->get_target_elements();
    std::vector<Element*> source_elements;

    unsigned int in_leaf_start = 0;
    //aggregate all source elements
    for(unsigned int i = 0; i<num_dir_neighbours; i++)
    {
        if (direct_neighbours[i] == target)
        {
            in_leaf_start = source_elements.size();
        }
        source_elements.insert(source_elements.end(),
                               direct_neighbours[i]->get_source_elements().begin(),
                               direct_neighbours[i]->get_source_elements().end());
    }
    unsigned int num_tgt_el = target_elements.size();
    unsigned int num_src_el = source_elements.size();
    arma::mat prec_block(num_tgt_el,num_tgt_el);
    bool make_prec = m_make_prec && target->is_leaf();
    for(unsigned int i = 0; i < num_tgt_el; i++)
    {
        double contrib = 0;
        double coeff;
        Element * t = target_elements[i];
        for(unsigned int j = 0; j < num_src_el; j++)
        {
            Element * s = source_elements[j];
            coeff = m_kernel->direct(*t,*s);
            if (make_prec)
            {
                // IMPORTANT:
                // we know that the order within the preconditioner is the order
                // of the elements within the leaf. 
                if(j>= in_leaf_start && j<in_leaf_start + num_tgt_el)
                {
                    prec_block(i,j-in_leaf_start) = coeff;
                }
            }
            contrib += s->get_value()*coeff;
        }
        t->set_target_value(t->get_target_value()+contrib);
    }
    if(make_prec)
    {
        // set block matrix in preconditioner        
        m_precond[target->get_leaf_number()] = prec_block;
    }
}

void FMM2D::evaluate_far_interactions(Cell* cell)
{
    std::vector<Element*> const & target_elements = cell->get_target_elements();
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
    
    if(m_make_prec)
    {
        if(!m_precond.size())
        { 
            //blocks within block matrix are ordered in order of leaf array
            std::vector<Cell*> const & leaves = m_tree->get_leaves();
            unsigned int num_leaves = leaves.size();
            m_precond.resize(num_leaves);
            m_prec_block_starts.resize(num_leaves,0);
            for(int i = 0; i<num_leaves; i++)
            {
                m_prec_block_starts[i] = leaves[i]->get_leaf_block_start_pos();
            }
        }
    }
    
    // level 2 only M2L and direct, no L2L
    Cell* cur_cell;
    
    while(it->has_next())
    {
        cur_cell = it->next();
        if(cur_cell->get_level() > 2) break;
        if(cur_cell->get_target_elements().empty())
        {
            continue;
        }
#ifdef DEBUG
    std::cout << "Downward pass: " << cur_cell->debug_info() << std::endl;
#endif
        
        //M2L
        m2l_downward_pass(cur_cell);
        direct_downward_pass(cur_cell);
        // leaf at lvl 2
        if(cur_cell->is_leaf())
        {
            evaluate_far_interactions(cur_cell);
        }
    }

    // there are no level 3 cells
    if(cur_cell->get_level() < 3) 
    {
        delete it;
        return;
    }
    
    //reached lvl 3 or higher --> use L2L and M2L
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
        m2l_downward_pass(cur_cell);
        l2l_downward_pass(cur_cell);

        // because of Liu's direct lists we need to do this for each cell
        direct_downward_pass(cur_cell);

        if(cur_cell->is_leaf())
        {
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

void FMM2D::evaluate()
{
    // in Liu's version we do the direct evaluation while we are
    // moving the tree downward because of the way the interaction
    // and direct lists are built
}

void FMM2D::reset()
{
    //set target values to 0
    unsigned int num_tgt_el = m_tgt_elements.size();
    for(unsigned int i = 0; i<num_tgt_el; i++)
    {
        m_tgt_elements[i]->set_target_value(0);
    }
    
    //set computed moments and local expansions to 0 for each cell
    if(!m_tree) return;
    std::vector<Cell*> const & cells = m_tree->get_cells();
    unsigned int num_cells = cells.size();
    for(unsigned int i = 0; i<num_cells; i++)
    {
        Cell* cur_cell = cells[i];
        cur_cell->set_local_exps(std::vector<double>());
        cur_cell->set_local_exps_cmp(std::vector<complex_t>());
        cur_cell->set_moments(std::vector<double>());
        cur_cell->set_moments_cmp(std::vector<complex_t>());
    }
}