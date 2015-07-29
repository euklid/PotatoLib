/* 
 * File:   cell2d.h
 * Author: euklid
 *
 * Created on June 9, 2015, 4:35 PM
 */

#ifndef CELL2D_H
#define	CELL2D_H

#include "cell.h"

class Cell2D : public Cell   {
public:
    Cell2D(unsigned int terms, double size, Point const & center);

    /**
     *  divides the cell into 4 cells, sets them as children and cell is
     *  their father cell. Children cells contain elements that lie within them,
     *  such that every element lies in exact one cell
     *  @return children
     */
    virtual std::vector<Cell*> divide();
    virtual bool contains_point(Point const & pt) const;
private:

};

#endif	/* CELL2D_H */

