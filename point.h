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
    double & operator[](unsigned int coord);
    const double & operator[](unsigned int coord) const;
    unsigned int get_dimension() const;
    Point operator-(Point const & p1) const;
    Point operator+(Point const & p1) const;
    

private:
    std::vector<double> m_coord;
    unsigned int m_dim;
};

#endif
