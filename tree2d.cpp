/* 
 * File:   tree2d.cpp
 * Author: euklid
 * 
 * Created on June 13, 2015, 2:55 PM
 */

#include "tree2d.h"
#include <queue>

Tree2D::Tree2D() : Tree(2), m_built(false)
{}

Tree2D::~Tree2D()
{
    for (int i = 0; i<m_cells.size(); i++) {
        delete m_cells[i];
    }
    m_cells.clear();
}

void Tree2D::build_tree(int max_elements)
{
    int index = 0;
    if (m_built)
    {
        return;
    }
    else
    {
        m_built = true;
    }
    
    std::queue<Cell*> undivided_cells;
    undivided_cells.push(m_root);
    m_root->set_level(1);
    while (!undivided_cells.empty()) 
    {
        Cell* cur_cell = undivided_cells.front();
        undivided_cells.pop();
        
        m_cells.push_back(cur_cell);
        cur_cell->set_id(index);
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

Tree_Iterator Tree2D::bfs_iterator()
{
    return Tree2D_BFS_Iterator(this);
}

Tree_Iterator Tree2D::upward_iterator()
{
    return Tree2D_Upward_Iterator(this);
}

Tree_Iterator Tree2D::downward_iterator()
{
    return Tree2D_BFS_Iterator(this);
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
    m_used_ids.at(m_last->get_id());
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