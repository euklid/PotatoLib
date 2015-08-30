/* 
 * File:   tree2d_ada.h
 * Author: euklid
 *
 * Created on August 28, 2015, 12:57 AM
 */

#ifndef TREE2D_ADA_H
#define	TREE2D_ADA_H

#include "tree2d.h"

class Tree2D_Ada : public Tree2D
{
public:
    Tree2D_Ada();
    virtual ~Tree2D_Ada();
        
protected:
    virtual void generate_interaction_lists();
    
private:
    void generate_lists134(Cell * const cell, Cell * const neighbor);
        

};

#endif	/* TREE2D_ADA_H */

