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

void Tree2D::build_tree(int max_elements)
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
    
    //creating new cells by dividing them according to division criterion
    generate_cells(max_elements);
    //generating interaction lists
    generate_interaction_lists();
}

void Tree2D::generate_cells(int max_elements)
{
    int index = 0;
    std::queue<Cell*> undivided_cells;
    undivided_cells.push(m_root);
    m_root->set_level(0);
    
    m_max_level = 0;
    m_lvl_ids.push_back(std::vector<unsigned int>());
    
    while (!undivided_cells.empty()) 
    {
        Cell* cur_cell = undivided_cells.front();
        undivided_cells.pop();
        
        m_cells.push_back(cur_cell);
        cur_cell->set_id(index);
        if (cur_cell->get_level() > m_max_level) {
            assert(cur_cell->get_level() == m_max_level+1);
            m_max_level++;
            m_lvl_ids.push_back(std::vector<unsigned int>());
        }
        
        m_lvl_ids[cur_cell->get_level()].push_back(cur_cell->get_id());
        
        assert(cur_cell->get_id() == m_cells.size()-1);
        index++;
        
		if(cur_cell->number_of_elements() > max_elements)
        {
            std::vector<Cell*> new_cells = cur_cell->divide();
            for (int i = 0; i<new_cells.size(); i++)
            {
                new_cells[i]->set_level(cur_cell->get_level()+1);
                undivided_cells.push(new_cells[i]);
            }
        }
        else
        {
            m_leaves.push_back(cur_cell);
        }
    }
}

void Tree2D::generate_interaction_lists()
{
    // first only like the Liu book
    for (unsigned int i = 0; i<=m_max_level; i++) {
        std::vector<unsigned int> & lvl_cells = m_lvl_ids[i];
        std::vector<unsigned int >::iterator it = lvl_cells.begin();
        std::vector<unsigned int >::const_iterator other_it;
        for (; it!=lvl_cells.end(); ++it) {
            Cell* cur_cell = m_cells[*it];
            assert(cur_cell->get_id() == *it);
            Point cur_cell_grid_pos = get_cell_grid_pos(cur_cell->get_center(),
                                                        i,
                                                        m_root->get_center(),
                                                        m_root->get_size());
            for (other_it = lvl_cells.begin(); other_it != lvl_cells.end(); ++other_it) {
                Cell* other_cell = m_cells[*other_it];
                Point other_cell_grid_pos = get_cell_grid_pos(other_cell->get_center(),
                                                              i,
                                                              m_root->get_center(),
                                                              m_root->get_size());
                if (Point::max_norm_dist(cur_cell_grid_pos, other_cell_grid_pos) != 1) {
                    continue;
                }
                
                //otherwise we have a neighbor
            }
        }
    }
}

Tree_Iterator *Tree2D::bfs_iterator()
{
    return new Tree2D_BFS_Iterator(this);
}

Tree_Iterator * Tree2D::upward_iterator()
{
    return new Tree2D_Upward_Iterator(this);
}

Tree_Iterator * Tree2D::downward_iterator()
{
    return new Tree2D_BFS_Iterator(this);
}

Tree2D_BFS_Iterator::Tree2D_BFS_Iterator(Tree2D* tree)
{
    m_last = NULL;
    m_tree = tree;
    m_cell_queue.push(tree->get_root());
}

Cell* Tree2D_BFS_Iterator::next()
{
    if (!has_next())
    {
        return NULL;
    }
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

Tree2D_Upward_Iterator::Tree2D_Upward_Iterator(Tree2D* tree)
{
    m_last = NULL;
    m_tree = tree;
    const std::vector<Cell*> leaf_cells = tree->get_leaves();
    for (int i = 0; i<leaf_cells.size(); i++)
    {
        m_leaf_queue.push(leaf_cells.at(i));
    }
    m_used_ids = std::vector<bool>(tree->get_cells().size(),false);
}

Cell* Tree2D_Upward_Iterator::next()
{
    m_last = m_leaf_queue.front();
    m_leaf_queue.pop();
    m_used_ids[(m_last->get_id())] = true;
    if (!m_used_ids.at(m_last->get_father()->get_id())) {
        Cell * const father = m_last->get_father();
        m_leaf_queue.push(father);
    }
    return m_last;
}

bool Tree2D_Upward_Iterator::has_next()
{
    return !m_leaf_queue.empty();
}
