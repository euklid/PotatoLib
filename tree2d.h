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
    virtual void build_tree(int max_elements);
    virtual Tree_Iterator* upward_iterator();
    virtual Tree_Iterator* downward_iterator();
    
private:
    virtual Tree_Iterator* bfs_iterator();
    void generate_cells(int max_elements);
    void generate_interaction_lists();
    bool m_built;
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


#endif	/* TREE2D_H */

