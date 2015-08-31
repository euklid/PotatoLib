/* 
 * File:   tree2d.h
 * Author: euklid
 *
 * Created on June 13, 2015, 2:55 PM
 */

#ifndef TREE2D_H
#define	TREE2D_H
#include "tree.h"
#include <queue>

class Tree2D : public Tree 
{
public:
    Tree2D();
    virtual ~Tree2D();

    /**
     *  Builds 2-dimensional tree from given root cell. Can only be
     *  called once, further calls do nothing
     *  @param max_elements maximal number of elements within a cell
     */
    virtual void build_tree(int max_elements, int min_level = 2);
    virtual Tree_Iterator* upward_iterator();
    virtual Tree_Iterator* downward_iterator();
    
    static Point lazy_get_set_cell_grid_pos(Cell* cell,
                                            unsigned int lvl,
                                            Cell const * const root);
    
protected:
    virtual Tree_Iterator* bfs_iterator();
    void generate_cells(int max_elements, int min_level);
    virtual void generate_interaction_lists();
    virtual void init_first_level_lists();
    virtual void make_other_levels_lists();
    bool m_built;
    int m_min_level;
};

class Tree2D_BFS_Iterator : public Tree_Iterator
{
public:
    Tree2D_BFS_Iterator(Tree2D* tree);
    virtual Cell* next();
    virtual bool has_next();
private:
    Tree2D* m_tree;
    Cell* m_last;
    std::queue<Cell*> m_cell_queue;
};

class Tree2D_Upward_Iterator : public Tree_Iterator
{
public:
    Tree2D_Upward_Iterator(Tree2D* tree);
    virtual Cell* next();
    virtual bool has_next();
private:
    Tree2D* m_tree;
    Cell* m_last;
    std::queue<Cell*> m_leaf_queue;
    std::vector<bool> m_used_ids;
};

class Tree2D_Upward_Code_Iterator : public Tree_Iterator
{
public:
    Tree2D_Upward_Code_Iterator(std::vector<std::vector<unsigned int> > const & lvl_ids,
                                std::vector<Cell*> const & cells,
                                int min_level);
    virtual Cell* next();
    virtual bool has_next();
private:
    std::vector<std::vector<unsigned int> > const & m_lvl_ids;
    long m_outer_size;
    long int m_outer_idx;
    long int m_inner_size;
    long int m_inner_idx;
    int m_min_level;
    std::vector<Cell*> const & m_cells;
};

class Tree2D_Downward_Code_Iterator : public Tree_Iterator
{
public:
    Tree2D_Downward_Code_Iterator(std::vector<std::vector<unsigned int> > const & lvl_ids,
                                  std::vector<Cell*> const & cells,
                                  int min_level);
    virtual Cell* next();
    virtual bool has_next();
private:
    std::vector<std::vector<unsigned int> > const & m_lvl_ids;
    long int m_outer_size;
    long int m_outer_idx;
    long int m_inner_size;
    long int m_inner_idx;
    int m_min_level;
    std::vector<Cell*> const & m_cells;
};

#endif	/* TREE2D_H */

