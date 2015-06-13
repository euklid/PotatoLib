/* 
 * File:   cell2d.h
 * Author: euklid
 *
 * Created on June 9, 2015, 4:35 PM
 */

#ifndef CELL2D_H
#define	CELL2D_H

#include "cell.h"

namespace FMM2D {

class Cell2D : Cell   {
public:
	Cell2D();
	Cell2D(const Cell2D& orig);
	virtual ~Cell2D();
private:

};

}

#endif	/* CELL2D_H */

