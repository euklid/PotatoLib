/* 
 * File:   tree2d_ada.cpp
 * Author: euklid
 * 
 * Created on August 28, 2015, 12:57 AM
 */

#include "tree2d_ada.h"

#if DEBUG
#include <iostream>
#endif

Tree2D_Ada::Tree2D_Ada()
: Tree2D()
{}

Tree2D_Ada::~Tree2D_Ada()
 {
    for (int i = 0; i<m_cells.size(); i++)
    {
        delete m_cells[i];
    }
    m_cells.clear();
}

void Tree2D_Ada::generate_interaction_lists()
{
    // initialize lists for all cells
    unsigned int num_cells = m_cells.size();
    for(unsigned int i = 0; i<num_cells; i++)
    {
        m_cells[i]->init_empty_lists(4);
    }
    
    //iterate through the tree to set all lists
    for (unsigned int i = m_min_level; i <= m_max_level; i++)
    {
        std::vector<unsigned int> & lvl_cells = m_lvl_ids[i];
        std::vector<unsigned int> & lvl_father_cells = m_lvl_ids[i - 1];
        std::vector<unsigned int >::const_iterator it = lvl_cells.begin();
        std::vector<unsigned int >::const_iterator other_it;
        for (; it != lvl_cells.end(); ++it)
        {
            Cell* cur_cell = m_cells[*it];
            assert(cur_cell->get_id() == *it);

            Cell* cur_cell_father = cur_cell->get_father();

            Point cur_cell_grid_pos = lazy_get_set_cell_grid_pos(cur_cell, i, m_root);
            Point cur_cell_father_grid_pos = lazy_get_set_cell_grid_pos(cur_cell_father, i - 1, m_root);
            for (other_it = lvl_father_cells.begin(); other_it != lvl_father_cells.end(); ++other_it)
            {
                Cell* other_cell_father = m_cells[*other_it];
                Point other_cell_father_grid_pos = lazy_get_set_cell_grid_pos(other_cell_father,
                                                                              i - 1,
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
                                                                              i,
                                                                              m_root);
                    if (Point::max_norm_dist(children_cell_grid_pos, cur_cell_grid_pos) > 1)
                    {
                        // well seperated on same level as cur_cell --> list 2
                        cur_cell->add_to_list(*children_it, 1);
                    }
                    else if (cur_cell->is_leaf())
                    {
                        // the "neighbor" is cell itself --> add to list1
                        if(cur_cell == *children_it)
                        {
                            cur_cell->add_to_list(cur_cell,0);
                        }
                        else
                        {
                            // have a look at neighbors of b (and their 
                            // descendents to find out who belongs to list 1 
                            // (direct neighbor) or list 3
                            generate_lists134(cur_cell, *children_it);
                        }
                    }
                }
            }
            
#ifdef DEBUG
            std::cout << "Cell " << cur_cell->debug_info() << " lists are"  << std::endl;
            for(unsigned int list_nr = 0; list_nr < 4; list_nr++)
            {
                std::vector<unsigned int> list_ids = cur_cell->get_list_ids(list_nr);
                std::cout << "list " << list_nr << " : " << std::endl;
                for(unsigned int i = 0; i< list_ids.size(); i++)
                {
                    std::cout << list_ids[i] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "\n" << std::endl;
#endif
        }
    }
}

void Tree2D_Ada::generate_lists134(Cell * const cell, Cell * const neighbor)
{
    if(neighbor->is_leaf())
    {
        // don't go any further, just add to direct list
        cell->add_to_list(neighbor,0);
        neighbor->add_to_list(cell,0);
        return;
    }
    
    // the following two vectors determine the minimal and maximal positions
    // in every direction in the level grid of the cell in relative depth to 
    // the cells depth(== level)
    std::vector<Point> cell_min_pos;
    std::vector<Point> cell_max_pos;
    cell_min_pos.push_back(cell->get_level_grid_position());
    cell_max_pos.push_back(cell->get_level_grid_position());
    
    // start directly with neighbor's children
    std::vector<Cell*> stack(neighbor->get_children());
    
    while(!stack.empty())
    {
        Cell* cur_cell = stack.back();
        stack.pop_back();
        int cur_depth = cur_cell->get_level() - cell->get_level();
        
        // reached deeper level on which cell's grid_positions need to be 
        // computed
        if (cur_depth == cell_min_pos.size())
        {
            Point next_min_pos(cell_min_pos.at(cur_depth-1));
            Point next_max_pos(cell_max_pos.at(cur_depth-1));
            int cell_size = 1 << cur_depth;
            for(int i = 0; i<m_dim; i++)
            {
                next_min_pos[i] = (int)(2*next_min_pos[i]);
                next_max_pos[i] = next_min_pos[i] + cell_size - 1;
            }
            cell_min_pos.push_back(next_min_pos);
            cell_max_pos.push_back(next_max_pos);
        }
        
        // get coordinates of cur_cell and find out whether it is a direct
        // neighbor or not
        Point cur_cell_pos = lazy_get_set_cell_grid_pos(cur_cell,
                                                        cur_cell->get_level(),
                                                        m_root);
        
        // we have direct neighbor if in all dimensions we are within range
        // of minimum position lowered by one everywhere and max position 
        // heightened by one everywhere
        Point const & min_cell_pos = cell_min_pos[cur_depth];
        Point const & max_cell_pos = cell_max_pos[cur_depth];
        bool direct = true;
        for(int i = 0; (i<m_dim) && direct; i++)
        {
            direct = (0 <= (cur_cell_pos[i] - min_cell_pos[i] + 1));
            direct = direct && (0<= (max_cell_pos[i] - cur_cell_pos[i])+1);
        }
        if(direct)
        {
            // if we have a leaf that we are fine and add it to the direct list1
            if (cur_cell->is_leaf())
            {
                //cell is direct neighbor of cur_cell
                cur_cell->add_to_list(cell,0);
                //and vice versa
                cell->add_to_list(cur_cell,0);
            }
            else
            {
                //since it is direct but a leaf we need to examine its children
                std::vector<Cell*> const & children = cur_cell->get_children();
                stack.insert(stack.end(),children.begin(),children.end());
            }
        } else
        {
            // the cur_cell's father was a direct neighbor to cell, otherwise
            // it wouldn't be in the stack. but cur_cell is not a direct 
            // neighbor to cell --> cur_cell is in list 3 of cell and we don't
            // care for its children
            cell->add_to_list(cur_cell,2);
            cur_cell->add_to_list(cell,3);
        }
        
    }
    
}