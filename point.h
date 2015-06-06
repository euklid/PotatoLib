#ifndef POINT_H
#define POINT_H
#include <vector>
#include <assert.h>
#include <math.h>

class Point
{
public:

    Point();

    Point(unsigned int dimension);
    static double dist(Point const & p1, Point const & p2);
    inline const double operator[](unsigned int coord) const;
    inline double operator[](unsigned int coord);
    unsigned int get_dimension() const;

private:
    std::vector<double> m_coord;
    unsigned int m_dim;
};

#endif
