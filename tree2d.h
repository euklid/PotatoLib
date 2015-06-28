/* 
 * File:   tree2d.h
 * Author: euklid
 *
 * Created on June 13, 2015, 2:55 PM
 */

#ifndef TREE2D_H
#define	TREE2D_H
#include "tree.h"


class Tree2D : public Tree 
{
public:
    Tree2D();
    virtual ~Tree2D();

    virtual void build_tree(int max_elements);

private:

};

#endif	/* TREE2D_H */

