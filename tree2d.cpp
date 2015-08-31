/* 
 * File:   tree2d.cpp
 * Author: euklid
 * 
 * Created on June 13, 2015, 2:55 PM
 */

#include "tree2d.h"
#include <queue>
#include <iostream>
#include <cmath>

Tree2D::Tree2D() : Tree(2), m_built(false)
{}

Tree2D::~Tree2D()
{
    for (int i = 0; i<m_cells.size(); i++)
    {
        delete m_cells[i];
    }
    m_cells.clear();
}

void Tree2D::build_tree(int max_elements, int min_level)
{
    if (m_built)
    {
        std::cout << "Tree has already been built." << std::endl;
        return;
    }
    else
    {
        m_built = true;
    }
    
    m_min_level = min_level;
    //creating new cells by dividing them according to division criterion
    generate_cells(max_elements, min_level);
    //generating interaction lists
    generate_interaction_lists();
}

void Tree2D::generate_cells(int max_elements,int min_level)
{
    int index = 0;
    std::queue<Cell*> undivided_cells;
    undivided_cells.push(m_root);
    m_root->set_level(0);
    int leaf_elements = 0;
    
    m_max_level = 0;
    m_lvl_ids.push_back(std::vector<unsigned int>());
    
    while (!undivided_cells.empty()) 
    {
        Cell* cur_cell = undivided_cells.front();
        undivided_cells.pop();
        
        if(cur_cell->number_of_elements() == 0) continue;

        m_cells.push_back(cur_cell);
        cur_cell->set_id(index);
        if (cur_cell->get_level() > m_max_level) {
#if DEBUG
            assert(cur_cell->get_level() == m_max_level+1);
#endif
            m_max_level++;
            m_lvl_ids.push_back(std::vector<unsigned int>());
        }
        
        m_lvl_ids[cur_cell->get_level()].push_back(cur_cell->get_id());
#if DEBUG
        assert(cur_cell->get_id() == m_cells.size()-1);
        assert(m_cells.size()-1 == index);
#endif
        index++;
        
	if(cur_cell->number_of_elements() > max_elements || cur_cell->get_level()<min_level)
        {
            std::vector<Cell*> & new_cells = cur_cell->divide();
            for (int i = 0; i<new_cells.size(); i++)
            {
                new_cells[i]->set_level(cur_cell->get_level()+1);
                undivided_cells.push(new_cells[i]);
            }
        }
        else
        {
            cur_cell->set_leaf_block_start_pos(leaf_elements);
            cur_cell->set_leaf_number(m_leaves.size());
            leaf_elements += cur_cell->get_source_elements().size();
            m_leaves.push_back(cur_cell);
        }
    }
}

Point Tree2D::lazy_get_set_cell_grid_pos(Cell* cell,
                                        unsigned int lvl,
                                        Cell const * const root)
{
    Point cell_grid_pos(cell->get_dimension());
    if (cell->has_level_grid_position())
    {
        cell_grid_pos = cell->get_level_grid_position();
    }
    else
    {
        cell_grid_pos = Tree::get_cell_grid_pos(cell->get_center(),
                                                lvl,
                                                root->get_center(),
                                                root->get_size());
        cell->set_level_grid_position(cell_grid_pos);
    }
    return cell_grid_pos;
}



void Tree2D::init_first_level_lists()
{
    if(m_min_level == 0)
    {
        m_root->init_empty_lists(3);
        m_root->add_to_list(m_root,3);
        return;
    }
    std::vector<unsigned int> & lvl_cells = m_lvl_ids[m_min_level];
    std::vector<unsigned int> & lvl_father_cells = m_lvl_ids[m_min_level - 1];
    std::vector<unsigned int >::const_iterator it = lvl_cells.begin();
    std::vector<unsigned int >::const_iterator other_it;
    for (; it != lvl_cells.end(); ++it)
    {
        Cell* cur_cell = m_cells[*it];

        // list 1: direct neighbors,
        // list 2: interaction list,
        // list 3: colleagues on same level
        cur_cell->init_empty_lists(3);
        assert(cur_cell->get_id() == *it);

        Cell* cur_cell_father = cur_cell->get_father();

        Point cur_cell_grid_pos = lazy_get_set_cell_grid_pos(cur_cell, m_min_level, m_root);
        Point cur_cell_father_grid_pos = lazy_get_set_cell_grid_pos(cur_cell_father, m_min_level - 1, m_root);
        for (other_it = lvl_father_cells.begin(); other_it != lvl_father_cells.end(); ++other_it)
        {
            Cell* other_cell_father = m_cells[*other_it];
            Point other_cell_father_grid_pos = lazy_get_set_cell_grid_pos(other_cell_father,
                                                                          m_min_level - 1,
                                                                          m_root);
            if (Point::max_norm_dist(cur_cell_father_grid_pos, other_cell_father_grid_pos) > 1)
            {
                continue;
            }
            //otherwise we have neighboring father cells (or the own father)
            //and want to retrieve the other cell's children

            if (other_cell_father->is_leaf())
            {
                continue;
            }
            std::vector<Cell*> children = other_cell_father->get_children();
            std::vector<Cell*>::const_iterator children_it = children.begin();

            for (; children_it != children.end(); ++children_it)
            {
                Point children_cell_grid_pos = lazy_get_set_cell_grid_pos((*children_it),
                                                                          m_min_level,
                                                                          m_root);
                if (Point::max_norm_dist(children_cell_grid_pos, cur_cell_grid_pos) > 1)
                {
                    cur_cell->add_to_list(*children_it, 1);
                }
                else 
                {
                    if ((*children_it)->is_leaf() || cur_cell->is_leaf() || m_min_level == m_max_level)
                    {
                        cur_cell->add_to_list(*children_it, 0);
                    }
                    cur_cell->add_to_list(*children_it,2);
                }
            }
        }
#ifdef DEBUG
        std::cout << "Cell " << cur_cell->get_id() << " at level " << cur_cell->get_level() << " at " << cur_cell->get_level_grid_position() << " interaction list is" << std::endl;
        std::vector<unsigned int> interaction_list = cur_cell->get_list_ids(1);
        for (unsigned int i = 0; i < interaction_list.size(); i++)
        {
            std::cout << interaction_list[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Cell " << cur_cell->get_id() << " direct list is" << std::endl;
        std::vector<unsigned int> direct_list = cur_cell->get_list_ids(0);
        for (unsigned int i = 0; i < direct_list.size(); i++)
        {
            std::cout << direct_list[i] << " ";
        }
        std::cout << "\n" << std::endl;
#endif
    }
}

void Tree2D::make_other_levels_lists()
{
    for (unsigned int i = m_min_level+1; i<=m_max_level; i++)
    {
        std::vector<unsigned int> & lvl_cells = m_lvl_ids[i];
        std::vector<unsigned int >::const_iterator it = lvl_cells.begin();
        std::vector<Cell*>::const_iterator other_it;
        for (; it!=lvl_cells.end(); ++it)
        {
            Cell* cur_cell = m_cells[*it];
            
            // list 1: direct neighbors,
            // list 2: interaction list,
            // list 3: colleagues on same level
            cur_cell->init_empty_lists(3);
            assert(cur_cell->get_id() == *it);

            Cell* cur_cell_father = cur_cell->get_father();

            Point cur_cell_grid_pos = lazy_get_set_cell_grid_pos(cur_cell,i,m_root);
            std::vector<Cell*> father_coll = cur_cell->get_father()->get_list(2);
            for (other_it = father_coll.begin(); other_it != father_coll.end(); ++other_it)
            {
                Cell* other_cell_father = *other_it;
                                
                if (other_cell_father->is_leaf())
                {
                    continue;
                }
                std::vector<Cell*> const & children = other_cell_father->get_children();
                std::vector<Cell*>::const_iterator children_it = children.begin();

                for(;children_it != children.end(); ++children_it)
                {
                    Point children_cell_grid_pos = lazy_get_set_cell_grid_pos((*children_it),
                                                                     i,
                                                                     m_root);
                    if (Point::max_norm_dist(children_cell_grid_pos,cur_cell_grid_pos) > 1)
                    {
                        cur_cell->add_to_list(*children_it,1);
                    }
                    else 
                    {
                        if ((*children_it)->is_leaf() || cur_cell->is_leaf() || i == m_max_level)
                        {
                            cur_cell->add_to_list(*children_it,0);
                        }
                        cur_cell->add_to_list(*children_it,2);
                    }  
                }
            }
#ifdef DEBUG
            std::cout << "Cell " << cur_cell->get_id() << " at level "<< cur_cell->get_level() <<" at " << cur_cell->get_level_grid_position() << " interaction list is" << std::endl;
            std::vector<unsigned int> interaction_list = cur_cell->get_list_ids(1);
            for(unsigned int i = 0; i< interaction_list.size(); i++)
            {
                std::cout << interaction_list[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "Cell " << cur_cell->get_id() << " direct list is" << std::endl;
            std::vector<unsigned int> direct_list = cur_cell->get_list_ids(0);
            for(unsigned int i = 0; i< direct_list.size(); i++)
            {
                std::cout << direct_list[i] << " ";
            }
            std::cout << "\n" << std::endl;
#endif
        }
    }
}

void Tree2D::generate_interaction_lists()
{   
    init_first_level_lists();
    make_other_levels_lists();
}

Tree_Iterator *Tree2D::bfs_iterator()
{
    return new Tree2D_BFS_Iterator(this);
}

Tree_Iterator * Tree2D::upward_iterator()
{
    return new Tree2D_Upward_Code_Iterator(m_lvl_ids, m_cells, m_min_level);
}

Tree_Iterator * Tree2D::downward_iterator()
{
    return new Tree2D_Downward_Code_Iterator(m_lvl_ids, m_cells, m_min_level);
}

Tree2D_BFS_Iterator::Tree2D_BFS_Iterator(Tree2D* tree)
{
    m_last = NULL;
    m_tree = tree;
    m_cell_queue.push(tree->get_root());
}

Cell* Tree2D_BFS_Iterator::next()
{
    m_last = m_cell_queue.front();
    m_cell_queue.pop();
    std::vector<Cell*> children = m_last->get_children();
    for (int i = 0; i<children.size(); i++)
    {
        m_cell_queue.push(children.at(i));
    }
    return m_last;
}

bool Tree2D_BFS_Iterator::has_next()
{
    return !m_cell_queue.empty();
}

Tree2D_Upward_Code_Iterator::Tree2D_Upward_Code_Iterator(
        std::vector<std::vector<unsigned int> > const & lvl_ids,
        std::vector<Cell*> const & cells,
        int min_level)
    :  m_lvl_ids(lvl_ids), m_cells(cells), m_min_level(min_level)
{
    m_outer_size = m_lvl_ids.size();
    m_outer_idx = m_outer_size-1;
    if(m_outer_idx > -1)
    {
        m_inner_size = m_lvl_ids[m_outer_idx].size();
    }
    m_inner_idx = -1;
}

Cell* Tree2D_Upward_Code_Iterator::next()
{
    ++m_inner_idx;
    if (m_inner_idx == m_inner_size)
    {
        m_inner_idx = 0;
        m_outer_idx--;
        m_inner_size = m_lvl_ids[m_outer_idx].size();
    }
#if DEBUG
    assert(m_cells.size() > m_lvl_ids[m_outer_idx][m_inner_idx]);
#endif
    return m_cells[m_lvl_ids[m_outer_idx][m_inner_idx]];
}

bool Tree2D_Upward_Code_Iterator::has_next()
{
    if (m_outer_idx > m_min_level)
    {
        return true;
    }
    // we are at highest level
    return (m_inner_idx+1!=m_inner_size);
}

Tree2D_Downward_Code_Iterator::Tree2D_Downward_Code_Iterator(
        const std::vector<std::vector<unsigned int> > &lvl_ids,
        const std::vector<Cell *> &cells,
        int min_level)
    :  m_lvl_ids(lvl_ids), m_cells(cells), m_min_level(min_level)
{
    m_outer_size = m_lvl_ids.size();
    m_outer_idx = min_level;
    if(m_outer_idx < m_outer_size)
    {
        m_inner_size = m_lvl_ids[m_outer_idx].size();
    }
    m_inner_idx = -1;
}

Cell *Tree2D_Downward_Code_Iterator::next()
{
    ++m_inner_idx;
    if (m_inner_idx == m_inner_size)
    {
        m_inner_idx = 0;
        m_outer_idx++;
        m_inner_size = m_lvl_ids[m_outer_idx].size();
    }
#if DEBUG
    assert(m_cells.size() > m_lvl_ids[m_outer_idx][m_inner_idx]);
#endif
    return m_cells[m_lvl_ids[m_outer_idx][m_inner_idx]];
}

bool Tree2D_Downward_Code_Iterator::has_next()
{
    if (m_outer_idx < m_outer_size-1)
    {
        return true;
    }
    // we are at highest level
    return (m_inner_idx+1!=m_inner_size);
}
