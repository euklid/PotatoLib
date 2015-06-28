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
    virtual void calc_moment();
    virtual std::vector<Cell*> divide();
    virtual bool contains_point(Point const & pt) const;
private:

};

#endif	/* CELL2D_H */

