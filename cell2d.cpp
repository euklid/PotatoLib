/* 
 * File:   cell2d.cpp
 * Author: euklid
 * 
 * Created on June 9, 2015, 4:35 PM
 */

#include "cell2d.h"


Cell2D::Cell2D(unsigned int terms, double size, Point const & center)
    : Cell(2, terms, size, center) {}

void Cell2D::calc_moment() {
    
}

std::vector<Cell*> Cell2D::divide() {
    
    
}
